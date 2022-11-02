include("./generateAdj.jl")
using .AdjGenerator, Graphs, Distributions, SparseArrays

# generate Adjacency matrix
numNodes = 50;

graphType = "complete";
graphParams = [];

graph = AdjGenerator.generateAdj(numNodes, graphType, graphParams);
Adj = adjacency_matrix(graph);

# define contagion dynamics parameters
lambda = 1;
gamma = 0.1;

modelParams = [lambda, gamma];

# time parameters
maxTime = 12;
timeRes = 0.01;

maxTime = timeRes*ceil(maxTime/timeRes);
t = Array(collect(Float64, 0:timeRes:maxTime)');
numTimes = length(t);

# initial conditions with one infected node
initConds = zeros(numNodes, 1);
initConds[1] = 1;

# Gillespie model parameters
numRuns = 10^3;

storedStates = zeros(numNodes, numTimes, numRuns);

for kRuns = 1:numRuns

    currTime = 0;
    currState = copy(initConds);
    tPosition = 1;
    storedStates[:, tPosition, kRuns] = copy(currState);

    infNodes = vec(initConds.==1);
    numInfNodes = Int(sum(infNodes)); 
    infNodesList = zeros(Int, numNodes, 1);
    infNodesList[1:numInfNodes] = findall(infNodes);

    possLinks = copy(Adj);
    possLinks[:, infNodesList[1:numInfNodes]].=0;
    possLinks = dropzeros(possLinks);

    vulNodes = possLinks[infNodesList[1:numInfNodes], :];
    vulNodesPairs = Tuple.(findall(x->x==1, vulNodes));

    vulNodesList = [];

    for i = 1:length(vulNodesPairs)
        push!(vulNodesList, vulNodesPairs[i][2]);
    end

    while currTime < maxTime

        oldState = copy(currState);

        infRate = lambda*length(vulNodesList);
        recRate = gamma*sum(infNodes);
        eventRate = infRate + recRate;

        timeStep = rand(Exponential(1/eventRate));

        randEventProb = rand();

        if randEventProb < infRate/eventRate

            newlyInfNode = rand(vulNodesList);
            currState[newlyInfNode] = 1;
            infNodes[newlyInfNode] = true;
            numInfNodes += 1;
            infNodesList[numInfNodes] = newlyInfNode;

            possLinks[:, newlyInfNode] .= 0;
            possLinks = dropzeros(possLinks);
                
            vulNodes = possLinks[infNodesList[1:numInfNodes], :];
            vulNodesPairs = Tuple.(findall(x->x==1, vulNodes));

            vulNodesList = [];

            for i = 1:length(vulNodesPairs)
                push!(vulNodesList, vulNodesPairs[i][2]);
            end
        
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

            vulNodes = possLinks[infNodesList[1:numInfNodes], :];
            vulNodesPairs = Tuple.(findall(x->x==1, vulNodes));

            vulNodesList = [];

            for i = 1:length(vulNodesPairs)
                push!(vulNodesList, vulNodesPairs[i][2]);
            end
        
        end

        currTime += timeStep;

        if currTime > maxTime
            currTime = maxTime;
            tPositionNew = numTimes;
        else
            tPositionNew = Int(floor(currTime/timeRes));
        end

        if tPositionNew > tPosition
            oldStateRepeated = oldState*ones(1, tPositionNew - tPosition);
            storedStates[:, (tPosition+1):tPositionNew, kRuns] = copy(oldStateRepeated);
            tPosition = tPositionNew;
        end

        if sum(infNodes) < 1
            currTime = maxTime;
        end

    end

    if tPosition < numTimes
        endStateRepeated = currState*ones(1, numTimes - tPosition);
        storedStates[:, (tPosition+1):numTimes, kRuns] = copy(endStateRepeated);
    end

end

storedStates

probS = sum(storedStates .== 0, dims=3)/numRuns
probI = sum(storedStates .== 1, dims=3)/numRuns
probR = sum(storedStates .== 2, dims=3)/numRuns