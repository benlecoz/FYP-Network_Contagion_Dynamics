module GillespieSIR

include("./generateAdj.jl")
using .AdjGenerator, Graphs, Distributions, SparseArrays
using Profile, PProf, BenchmarkTools

function SIRGillespie(numNodes::Int64, graphType::String, graphParams::Array, modelParams::Array, maxTime::Float64, 
                        timeRes::Float64, initConds::Array, numRuns::Int64)

    # Graph, adjacency matrix and edge array
    graph = AdjGenerator.generateAdj(numNodes, graphType, graphParams);
    Adj = adjacency_matrix(graph);
    edgeArray = Tuple.(edges(graph));

    # Time components
    maxTime = timeRes*ceil(maxTime/timeRes);
    t = Array(collect(Float64, 0:timeRes:maxTime)');
    numTimes = length(t);

    # Model parameters
    lambda = modelParams[1];
    gamma = modelParams[2];

    # Number of nodes
    numEdges = length(edgeArray);

    # Arrays of individual and pair states that will be used and reset at every run 
    storedStates = zeros(Int64, numNodes, numTimes);
    storedPairStates = zeros(Int64, numEdges, numTimes);

    # Individual and pair state arrays that will store cumulative runs

    # In all runs, indexing for individual states is: S=0 I=1 R=2
    storedS, storedI, storedR = [zeros(Int64, numNodes, numTimes) for _ = 1:3]
    # In all runs, indexing for pair states: SS=0 SI=1 SR=2  IS=3 II=4 IR=5  RS=6 RI=7 RR=8
    storedSS, storedSI, storedSR, storedIS, storedII, storedIR, storedRS, storedRI, storedRR = [zeros(Int64, numEdges, numTimes) for _ = 1:9]

    # Find edge number of specific node
    # Function allows input on whether the node in question is on the LHS, RHS or both sides of the pair making up the edge 
    function edgeArrayGen(node, side)
        if side == "both"
            edgesNode = (findall(x->x==node, first.(edgeArray)))..., (findall(x->x==node, last.(edgeArray))...)
        elseif side == "LHS"
            edgesNode = findall(x->x==node, first.(edgeArray))
        elseif side == "RHS"
            edgesNode = findall(x->x==node, last.(edgeArray))
        end

        return edgesNode
    end

    # Set all edges as S-S, except the edges connecting edge 1 to other nodes
    storedSS[:, 1].=1;
    storedSS[edgeArrayGen(1, "LHS"), 1].=0;

    # Set all edges coming from node 1 to S-I
    storedSI[edgeArrayGen(1, "LHS"), 1].=1;

    # Set initial pair-wise conditions, using pair indexing S-S=0 and S-I=1
    initPairConds = storedSI[:, 1];

    @time for kRuns = 1:numRuns 

        # Current time and position in time array
        currTime = 0.0;
        tPosition = 1;

        # Set/Reset current state as the initial individual and pair-based states
        currState = copy(initConds);
        currPairState = copy(initPairConds)

        # Set/Reset stored states arrays to initial conditions
        storedStates[:, tPosition] = currState;
        storedPairStates[:, tPosition] = currPairState;

        # List of infected nodes
        infNodes = vec(initConds.==1);
        numInfNodes = sum(infNodes);
        infNodesList = zeros(Int, numNodes, 1);
        infNodesList[1:numInfNodes] = findall(infNodes);

        # Links between infected and susceptible nodes in adjacency matrix
        possLinks = copy(Adj);
        possLinks[:, infNodesList[1:numInfNodes]].=0;
        possLinks = dropzeros(possLinks);

        # Index of S-I pairings
        vulNodes = possLinks[infNodesList[1:numInfNodes], :];
        vulNodesPairs = Tuple.(findall(x->x==1, vulNodes));

        # Iterate Gillespie algorithm till the max time is reached
        while currTime < maxTime

            # Define current state before any events take place, to copy over Tau-leaps that were skipped
            oldState = copy(currState);
            oldPairState = copy(currPairState);

            # Probability of an event happening, infection or recovery, based on the number of susceptible and infected nodes
            infRate = lambda*size(vulNodesPairs)[1];
            recRate = gamma*sum(infNodes);
            eventRate = infRate + recRate;

            # Next timestep at which an event happens
            timeStep = rand(Exponential(1/eventRate));

            # Random value from 0 to 1
            randEventProb = rand();

            # if the random value is within the probability that the event occurring is an infection, then a node is infected
            if randEventProb < infRate/eventRate

                # Choose a random infected node, update current state and infected nodes list to infected
                newlyInfNode = rand(vulNodesPairs)[2];
                currState[newlyInfNode] = 1;
                infNodes[newlyInfNode] = true;
                numInfNodes += 1;
                infNodesList[numInfNodes] = newlyInfNode;

                # Possible infection links to the newly infected node are removed from the adjacency matrix
                possLinks[:, newlyInfNode] .= 0;
                possLinks = dropzeros(possLinks);

                # Associated pair states are updated to take into account node switching from S to I
                currPairState[edgeArrayGen(newlyInfNode, "LHS")] .+= 3;
                currPairState[edgeArrayGen(newlyInfNode, "RHS")] .+= 1;
            
            # Otherwise a node recovers
            else
                
                # Choose recovered node
                chosenRecNodeIndex = rand(1:numInfNodes);
                newlyRecNode = infNodesList[chosenRecNodeIndex];

                # Update list and number of infected nodes
                infNodesList[chosenRecNodeIndex] = infNodesList[numInfNodes];
                infNodesList[numInfNodes] = 0;
                currState[newlyRecNode] = 2;
                infNodes[newlyRecNode] = false;
                numInfNodes -= 1;

                # remove possible infection links coming out of the newly recovered node
                possLinks[newlyRecNode, :] .= 0;
                possLinks = dropzeros(possLinks);

                # Associated pair states are updated to take into account node switching from I to R
                currPairState[edgeArrayGen(newlyRecNode, "LHS")] .+= 3;
                currPairState[edgeArrayGen(newlyRecNode, "RHS")] .+= 1;
            
            end
            
            # Update vulnerable node list to either add nodes connected to newly infected node, or remove infection links from newly recovered node  
            vulNodes = possLinks[infNodesList[1:numInfNodes], :];
            vulNodesPairs = Tuple.(findall(x->x==1, vulNodes));

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
                oldStateRepeated = oldState*ones(Int64, 1, tPositionNew - tPosition);
                storedStates[:, (tPosition+1):tPositionNew] = oldStateRepeated;

                oldPairStateRepeated = oldPairState*ones(Int64, 1, tPositionNew - tPosition);
                storedPairStates[:, (tPosition+1):tPositionNew] = oldPairStateRepeated;

                tPosition = tPositionNew;
            end

            # Stop simulation when there are no more infected nodes left
            if sum(infNodes) < 1
                currTime = maxTime;
            end

        end

        # Repeat final state until the end of the time array if the simulation was stopped early
        if tPosition < numTimes
            endStateRepeated = currState*ones(Int64, 1, numTimes - tPosition);
            endPairStateRepeated = currPairState*ones(Int64, 1, numTimes - tPosition);

            storedStates[:, (tPosition+1):numTimes] = endStateRepeated;
            storedPairStates[:, (tPosition+1):numTimes] = endPairStateRepeated;
        end

        # Save first simulated run to boolean state arrays based on indexing 
        if kRuns == 1
            storedS, storedI, storedR = [storedStates.==i for i=0:2]
            storedSS, storedSI, storedSR, storedIS, storedII, storedIR, storedRS, storedRI, storedRR  = [storedPairStates.==i for i=0:8]
        
        # Add subsequent runs to stored boolean state arrays
        else

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
    end

    # Calculate individual state probabilities averaged over the number of runs
    storedAll = [storedS, storedI, storedR];
    probS, probI, probR = [storedState/numRuns for storedState in storedAll];

    # Calculate pair state probabilities averaged over the number of runs
    storedPairAll = [storedSS, storedSI, storedSR, storedIS, storedII, storedIR, storedRS, storedRI, storedRR];
    probSS, probSI, probSR, probIS, probII, probIR, probRS, probRI, probRR = [storedPairState/numRuns for storedPairState in storedPairAll];

    return probS, probI, probR, probSS, probSI, probSR, probIS, probII, probIR, probRS, probRI, probRR, t

end

end