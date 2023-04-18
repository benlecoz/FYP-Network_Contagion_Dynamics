module sirTreeApprox

function SIRTreeApprox(graph::SimpleGraph, depth::Array, modelParams::Array, maxTime::Float64, timeRes::Float64)

    # extract model parameters and the number of nodes from the inputs
    lambda, gamma = modelParams;
    numNodes = nv(graph);

    # extract time components
    maxTime = timeRes*ceil(maxTime/timeRes);
    t = Array(collect(Float64, 0:timeRes:maxTime)');
    numTimes = length(t);

    # set the size of the Susceptible, Infected and Recovered arrays
    sTree, iTree, rTree = [zeros(numNodes, numTimes) for _ = 1:3];

    # at every time point, calculate the three state values for each node 
    for (index, time) in enumerate(t)

        iTree[1, index] = exp(-gamma * time)
        rTree[1, index] = 1 - exp(-gamma * time)

        for j=2:numNodes
            if depth[j] == 1
                sTree[j, index] = 1 - lambda/(lambda + gamma) + lambda/(lambda + gamma) * exp(-(lambda + gamma) * time)
                iTree[j, index] = exp(-gamma * time) - exp(-(lambda + gamma) * time)
                rTree[j, index] = lambda/(lambda + gamma) - exp(-gamma * time) + exp(-(lambda + gamma) * time) * (1 - lambda/(lambda + gamma))
            else
                sTree[j, index] = 1 - lambda^depth[j]/(lambda + gamma)^depth[j] + lambda^depth[j]/(lambda + gamma)^depth[j] * exp(-(lambda + gamma) * time) * sum((lambda + gamma)^n * time^n / factorial(n) for n in 0:(depth[j] - 1))
                iTree[j, index] = exp(-gamma * time) - exp(-(lambda + gamma) * time) * sum(lambda^n * time^n / factorial(n) for n in 0:(depth[j] - 1))
                rTree[j, index] = lambda^depth[j]/(lambda + gamma)^depth[j] - exp(-gamma * time) + exp(-(lambda + gamma) * time) * sum((lambda^n - lambda^depth[j]/(lambda + gamma)^(depth[j]-n)) * time^n/factorial(n) for n in 0:(depth[j] - 1))
            end
        end
    end

    return sTree, iTree, rTree

end

end
