module sirGillespie

include("./generateAdj.jl")
using .AdjGenerator
using Graphs, Distributions, SparseArrays
using Profile, PProf, BenchmarkTools, InvertedIndices

function SIRGillespie(edgeArray::Array, graph::SimpleGraph, modelParams::Array, maxTime::Float64, timeRes::Float64, initInfNodes::Array, numRuns::Int64)

    # Extract time components as time array
    maxTime = timeRes*ceil(maxTime/timeRes);
    t = Array(collect(Float64, 0:timeRes:maxTime)');
    numTimes = length(t);

    # Model parameters for infection and recovery
    lambda, gamma = modelParams;

    # Number of nodes and edges
    numNodes = nv(graph);
    numEdges = length(edgeArray);

    # Create initial conditions vector from initially infected nodes
    initConds = zeros(Int64, numNodes, 1);
    initConds[initInfNodes] .= 1;

    # Arrays of individual and pair states that will be used and reset at every run 
    storedStates = zeros(Int64, numNodes, numTimes);
    storedPairStates = zeros(Int64, numEdges, numTimes);

    # Individual and pair state arrays that will store cumulative runs
    # In all runs, indexing for individual states is: S=0 I=1 R=2
    storedS, storedI, storedR = [zeros(Int64, numNodes, numTimes) for _ = 1:3];
    
    # In all runs, indexing for pair states: SS=0 SI=1 SR=2  IS=3 II=4 IR=5  RS=6 RI=7 RR=8
    storedSS, storedSI, storedSR, storedIS, storedII, storedIR, storedRS, storedRI, storedRR = [zeros(Int64, numEdges, numTimes) for _ = 1:9];

    # Set initial pair-wise conditions, using pair indexing S-S=0, S-I=1 and I-S=3
    initPairConds = vec(zeros(Int64, numEdges, 1));
    initPairConds[AdjGenerator.edgeArrayGenList(edgeArray, initInfNodes, "RHS")] .= 1;
    initPairConds[AdjGenerator.edgeArrayGenList(edgeArray, initInfNodes, "LHS")] .= 3;

    # Ensure that edges connecting two initially infected nodes are not accidentally marked as S-I or I-S, but instead I-I
    for i in AdjGenerator.edgeArrayGenList(edgeArray, initInfNodes, "both")
        if in.(first.(edgeArray)[i], Ref(initInfNodes)) & in.(last.(edgeArray)[i], Ref(initInfNodes))
            initPairConds[i] = 4;
        end
    end

    for kRuns = 1:numRuns

        # Current time and position in time array
        currTime = 0.0;
        tPosition = 1;

        # Set/Reset current state as the initial individual and pair-based states
        currState = copy(initConds);
        currPairState = copy(initPairConds);

        # Set/Reset stored states arrays to initial conditions
        storedStates[:, tPosition] = currState;
        storedPairStates[:, tPosition] = currPairState;

        # Initialise the number of infected nodes, which is used to test whether there any infected nodes remaining
        numInfNodes = sum(vec(initConds.==1));

        # Iterate Gillespie algorithm till the max time is reached
        while currTime < maxTime

            # Define current state before any events take place, to copy over Tau-leaps that might be skipped
            oldState = copy(currState);
            oldPairState = copy(currPairState);

            # Define the vulnerable edges as pairs of S-I or I-S nodes
            newVulEdges = sort(reduce(vcat, [findall(x->x==i, vec(currPairState)) for i in [1, 3]]));

            # Probability of an event happening, infection or recovery, based on the number of susceptible and infected nodes
            infRate = lambda*size(newVulEdges)[1];
            recRate = gamma*sum(vec(currState.==1));
            eventRate = infRate + recRate;

            # Next timestep at which an event happens
            timeStep = rand(Exponential(1/eventRate));

            # Random value from 0 to 1
            randEventProb = rand();

            # if the random value is within the probability that the event occurring is an infection, then a node is infected
            if randEventProb < infRate/eventRate

                # Choose a random infected edge, then choose extract newly infected node based on S-I or I-S status
                newlyInfEdge = rand(newVulEdges);
                newlyInfNode = currPairState[newlyInfEdge, 1] == 1 ? edgeArray[newlyInfEdge][1] : edgeArray[newlyInfEdge][2]
                
                # Update current state vector and number of infected nodes
                currState[newlyInfNode] = 1;
                numInfNodes += 1;

                # Associated pair states are updated to take into account node switching from S to I
                currPairState[AdjGenerator.edgeArrayGenList(edgeArray, newlyInfNode, "LHS")] .+= 3;
                currPairState[AdjGenerator.edgeArrayGenList(edgeArray, newlyInfNode, "RHS")] .+= 1;
            
            # Otherwise a node recovers
            else
                
                # Choose recovered node
                newlyRecNode = rand(findall(x->x==1, currState))[1]

                # Update current state vector and number of infected nodes
                currState[newlyRecNode] = 2;
                numInfNodes -= 1;

                # Associated pair states are updated to take into account node switching from I to R
                currPairState[AdjGenerator.edgeArrayGenList(edgeArray, newlyRecNode, "LHS")] .+= 3;
                currPairState[AdjGenerator.edgeArrayGenList(edgeArray, newlyRecNode, "RHS")] .+= 1;
            
            end

            # Update timestep
            currTime += timeStep;

            # Update position within time array, depending on if the max time has been reached or exceeded 
            if currTime > maxTime
                currTime = maxTime;
                tPositionNew = numTimes;
            else
                tPositionNew = floor(Int64, currTime/timeRes);
            end

            # If the event timestep is larger than designated timestep, save the state to stored array
            # If the event timestep is multiple times larger than designated timestep, save replicated states to stored array
            # If the event timestep is smaller than designated timestep, simulate more events until the designated timestep is reached
            if tPositionNew > tPosition
                storedStates[:, (tPosition+1):tPositionNew] .= oldState;
                storedPairStates[:, (tPosition+1):tPositionNew] .= oldPairState;

                tPosition = tPositionNew;
            end

            # Stop simulation when there are no more infected nodes left
            if numInfNodes < 1
                currTime = maxTime;
            end

        end

        # Repeat final state until the end of the time array if the simulation is stopped early
        if tPosition < numTimes
            storedStates[:, (tPosition+1):numTimes] .= currState;
            storedPairStates[:, (tPosition+1):numTimes] .= currPairState;
        end
        
        storedS += storedStates.==0;
        storedI += storedStates.==1;
        storedR += storedStates.==2;

        storedSS += storedPairStates.==0;
        storedSI += storedPairStates.==1;
        storedSR += storedPairStates.==2;
        storedIS += storedPairStates.==3;
        storedII += storedPairStates.==4;
        storedIR += storedPairStates.==5;
        storedRS += storedPairStates.==6;
        storedRI += storedPairStates.==7;
        storedRR += storedPairStates.==8;

    end

    # Calculate individual state probabilities averaged over the number of runs
    probS, probI, probR = [storedState/numRuns for storedState in [storedS, storedI, storedR]];

    # Calculate pair state probabilities averaged over the number of runs
    probSS, probSI, probSR, probIS, probII, probIR, probRS, probRI, probRR = [storedPairState/numRuns for storedPairState in [storedSS, storedSI, storedSR, storedIS, storedII, storedIR, storedRS, storedRI, storedRR]];

    return probS, probI, probR, probSS, probSI, probSR, probIS, probII, probIR, probRS, probRI, probRR, t

end

end