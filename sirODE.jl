include("./generateAdj.jl")
using .AdjGenerator, Graphs, SparseArrays
using DifferentialEquations, BenchmarkTools

numNodes = 100;

graphType = "ErdosRenyi";
graphParams = [0.05];

lambda = 1;
gamma = 0.1;
modelParams = [lambda, gamma];

maxTime = 12.0;
timeRes = 0.01;

initConds = zeros(Int64, numNodes, 2);
initConds[1, 2] = 1;
initConds[2:end, 1] .= 1;

numRuns = 10^3;

# function sirODE(numNodes::Int64, graphType::String, graphParams::Array, modelParams::Array, maxTime::Float64, 
#     timeRes::Float64, initConds::Array)

graph = AdjGenerator.generateAdj(numNodes, graphType, graphParams);
edgeArray = Tuple.(edges(graph));

maxTime = timeRes*ceil(maxTime/timeRes);
t = Array(collect(Float64, 0:timeRes:maxTime)');

function diffEq(du, u, p, t)

    modelRates = zeros(numNodes, numNodes);
    s = u[:, 1];
    i = u[:, 2];

    adjIndex = (first.(edgeArray) .- 1) * numNodes + last.(edgeArray);

    modelRates[adjIndex] = lambda .* s[first.(edgeArray)] .* i[last.(edgeArray)];

    totalInfRate = sum(modelRates);

    for k=1:numNodes
        du[k, 1] = -totalInfRate;
        du[k, 2] = totalInfRate - gamma*i[k]
    end

end;

prob = ODEProblem(diffEq, initConds, (0, maxTime));
sol = solve(prob)

# sSol = sol[1:numNodes,:]
# iSol = sol[(numNodes+1):2*numNodes,:]
# rSol = 1 .- sSol .- iSol

sSol = sol[:, 1, :]
iSol = sol[:, 2, :]
rSol = 1 .- sSol .- iSol

plot(iSol', label="")

