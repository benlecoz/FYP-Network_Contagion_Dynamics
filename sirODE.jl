include("./generateAdj.jl")
using .AdjGenerator, Graphs, SparseArrays
using DifferentialEquations, BenchmarkTools
using Plots
using InvertedIndices

numNodes = 100;

graphType = "ErdosRenyi";
graphParams = [0.05];

lambda = 1;
gamma = 0.1;
modelParams = [lambda, gamma];

maxTime = 50.0;
timeRes = 0.01;

# function sirODE(numNodes::Int64, graphType::String, graphParams::Array, modelParams::Array, maxTime::Float64, 
#     timeRes::Float64, initConds::Array)

graph = AdjGenerator.generateAdj(numNodes, graphType, graphParams);
edgeArray = Tuple.(edges(graph));

initInfNodes = [1, 3];

initConds = zeros(Int64, nv(graph), 2);
initConds[initInfNodes, 2] .= 1;
initConds[Not(initInfNodes), 1] .= 1;

maxTime = timeRes*ceil(maxTime/timeRes);
t = Array(collect(Float64, 0:timeRes:maxTime)');

function edgeArrayGenList(nodeList, side)
    if side == "both"
        edgesNodeList = unique([reduce(vcat, [findall(x->x==i, first.(edgeArray)) for i in nodeList])..., (reduce(vcat, [findall(x->x==i, last.(edgeArray)) for i in nodeList])...)])
    elseif side == "LHS"
        edgesNodeList = reduce(vcat, [findall(x->x==i, first.(edgeArray)) for i in nodeList])
    elseif side == "RHS"
        edgesNodeList = reduce(vcat, [findall(x->x==i, last.(edgeArray)) for i in nodeList])
    end

    return edgesNodeList
end;

u0 = [initConds[:, 1]; initConds[:, 2]];


function diffEq(du, u, p, t)

    s = u[1 : p[1]];
    i = u[(p[1] + 1) : 2*p[1]];

    adjIndex = (first.(edgeArray) .- 1) * p[1] + last.(edgeArray);

    infRates = zeros(p[1], p[1]);
    infRates[adjIndex] = lambda .* s[last.(edgeArray)] .* i[first.(edgeArray)];
    totalInfRate = sum(infRates, dims = 2);

    du[1 : p[1]] .= -totalInfRate;
    du[(p[1] + 1) : 2*p[1]] .= totalInfRate .- gamma*i

end;

p = [nv(graph), numEdges];

prob = ODEProblem(diffEq, u0, (0, maxTime), p, dt = 0.01, saveat=t);
sol = solve(prob);

sSol = sol[1:p[1], :];
iSol = sol[(p[1] + 1):2*p[1], :];
rSol = 1 .- sSol .- iSol;

plot(sSol', label="")
plot(iSol', label="")
plot(rSol', label="")

