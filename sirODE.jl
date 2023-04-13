module sirODE

using Graphs, SparseArrays, DifferentialEquations
using Plots, InvertedIndices, BenchmarkTools

function SIRODE(edgeArray::Array, graph::SimpleGraph, modelParams::Array, maxTime::Float64, timeRes::Float64, initInfNodes::Array)

    lambda, gamma = modelParams;
    numNodes = nv(graph);

    initConds = zeros(Int64, numNodes, 2);
    initConds[initInfNodes, 2] .= 1;
    initConds[Not(initInfNodes), 1] .= 1;

    maxTime = timeRes*ceil(maxTime/timeRes);
    t = Array(collect(Float64, 0:timeRes:maxTime)');

    u0 = [initConds[:, 1]; initConds[:, 2]];

    function diffEq(du, u, p, t)

        lambda, gamma = p

        s = u[1 : numNodes];
        i = u[(numNodes + 1) : 2*numNodes];

        adjIndex = (first.(edgeArray) .- 1) * numNodes + last.(edgeArray);

        infRates = zeros(numNodes, numNodes);
        infRates[adjIndex] = lambda .* s[last.(edgeArray)] .* i[first.(edgeArray)];
        totalInfRate = sum(infRates, dims = 2);

        du[1 : numNodes] .= -totalInfRate;
        du[(numNodes + 1) : 2*numNodes] .= totalInfRate .- gamma*i

    end;

    p = [lambda, gamma];

    prob = ODEProblem(diffEq, u0, (0, maxTime), p, dt = 0.01, saveat=t);
    sol = solve(prob);

    sSol = sol[1:numNodes, :];
    iSol = sol[(numNodes + 1):2*numNodes, :];
    rSol = 1 .- sSol .- iSol;

    return sSol, iSol, rSol

end

end