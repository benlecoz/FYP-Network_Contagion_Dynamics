module GillespieSIR

include("./generateAdj.jl")
using .AdjGenerator, Graphs, Distributions, SparseArrays
using Profile, PProf, BenchmarkTools

function SIRGillespie(numNodes::Int64, graphType::String, graphParams::Array, modelParams::Array, maxTime::Float64, 
                        timeRes::Float64, initConds::Array, numRuns::Int64)

    graph = AdjGenerator.generateAdj(numNodes, graphType, graphParams);
    Adj = adjacency_matrix(graph);

    maxTime = timeRes*ceil(maxTime/timeRes);
    t = Array(collect(Float64, 0:timeRes:maxTime)');
    numTimes = length(t);

    lambda = modelParams[1];
    gamma = modelParams[2];

    storedStates = zeros(Int64, numNodes, numTimes);

    storedS = zeros(Int64, numNodes, numTimes);
    storedI = zeros(Int64, numNodes, numTimes);
    storedR = zeros(Int64, numNodes, numTimes);

    @time for kRuns = 1:numRuns

        currTime = 0.0;
        currState = copy(initConds);
        tPosition = 1;
        storedStates[:, tPosition] = currState;

        infNodes = vec(initConds.==1);
        numInfNodes = sum(infNodes);
        infNodesList = zeros(Int, numNodes, 1);
        infNodesList[1:numInfNodes] = findall(infNodes);

        possLinks = copy(Adj);
        possLinks[:, infNodesList[1:numInfNodes]].=0;
        possLinks = dropzeros(possLinks);

        vulNodes = possLinks[infNodesList[1:numInfNodes], :];
        vulNodesPairs = Tuple.(findall(x->x==1, vulNodes));

        while currTime < maxTime

            oldState = copy(currState);

            infRate = lambda*size(vulNodesPairs)[1];
            recRate = gamma*sum(infNodes);
            eventRate = infRate + recRate;

            timeStep = rand(Exponential(1/eventRate));

            randEventProb = rand();

            if randEventProb < infRate/eventRate

                newlyInfNode = rand(vulNodesPairs)[2];
                currState[newlyInfNode] = 1;
                infNodes[newlyInfNode] = true;

                numInfNodes += 1;
                infNodesList[numInfNodes] = newlyInfNode;

                possLinks[:, newlyInfNode] .= 0;
                possLinks = dropzeros(possLinks);
            
            else

                chosenRecNodeIndex = rand(1:numInfNodes);
                newlyRecNode = infNodesList[chosenRecNodeIndex];

                infNodesList[chosenRecNodeIndex] = infNodesList[numInfNodes];
                infNodesList[numInfNodes] = 0;

                currState[newlyRecNode] = 2;
                infNodes[newlyRecNode] = false;
                numInfNodes -= 1;

                possLinks[newlyRecNode, :] .= 0;
                possLinks = dropzeros(possLinks);
            
            end

            vulNodes = possLinks[infNodesList[1:numInfNodes], :];
            vulNodesPairs = Tuple.(findall(x->x==1, vulNodes));

            currTime += timeStep;

            if currTime > maxTime
                currTime = maxTime;
                tPositionNew = numTimes;
            else
                tPositionNew = floor(Int64, currTime/timeRes);
            end

            if tPositionNew > tPosition
                oldStateRepeated = oldState*ones(Int64, 1, tPositionNew - tPosition);
                storedStates[:, (tPosition+1):tPositionNew] = oldStateRepeated;
                tPosition = tPositionNew;
            end

            if sum(infNodes) < 1
                currTime = maxTime;
            end

        end

        if tPosition < numTimes
            endStateRepeated = currState*ones(Int64, 1, numTimes - tPosition);
            storedStates[:, (tPosition+1):numTimes] = endStateRepeated;
        end

        if kRuns == 1

            storedS = storedStates.==0;
            storedI = storedStates.==1;
            storedR = storedStates.==2;         
        
        else

            storedS += storedStates.==0;
            storedI += storedStates.==1;
            storedR += storedStates.==2;
        
        end

    end

    probS = storedS/numRuns;
    probI = storedI/numRuns;
    probR = storedR/numRuns;

    return probS, probI, probR, t, graph

end

end