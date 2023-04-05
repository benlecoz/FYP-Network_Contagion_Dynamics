include("./generateAdj.jl")
using .AdjGenerator, Graphs, SparseArrays
using DifferentialEquations, BenchmarkTools
using Plots
using InvertedIndices

# number of nodes
numNodes = 100;

# graph type and associated parameters
graphType = "WattsStrogatz";
graphParams = [4, 0.5];

# model parameters
lambda = 1;
gamma = 0.1;
modelParams = [lambda, gamma];

# time parameters
maxTime = 50.0;
timeRes = 0.01;

# generate and its array of edges
graph = AdjGenerator.generateAdj(numNodes, graphType, graphParams);
edgeArray = Tuple.(edges(graph));

# edgeArray only has directed edges, so generate oppositely directed edges for every node pair
for (i, t) in pairs(edgeArray)
    push!(edgeArray, reverse(t))
end;
edgeArray = sort(edgeArray);

# set initially infected node(s)
initInfNodes = [1, 3];

# set the initial conditions of each nodes
initConds = zeros(Int64, nv(graph), 2);
initConds[initInfNodes, 2] .= 1;
initConds[Not(initInfNodes), 1] .= 1;

# extract time components
maxTime = timeRes*ceil(maxTime/timeRes);
t = Array(collect(Float64, 0:timeRes:maxTime)');

# generate list of neighbour vertices for each node, based on which side of the edgeArray node pair it is on 
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

# calculate number of directed edges 
numEdges = length(edgeArray);

# based on initially infected nodes and edgeArray, set the initial pair conditions of the S-S and S-I states
initPairConds = zeros(Int64, numEdges, 2);
initPairConds[edgeArrayGenList(initInfNodes, "RHS"), 2] .= 1;
initPairConds[Not(edgeArrayGenList(initInfNodes, "both")), 1] .= 1;

# set the parameters associated with the ODEProblem
p = [nv(graph), numEdges];

# for each node, find the edge pairs for which the node is on the LHS (will be used to find which nodes are in the S-I state)  
edgeArrayIndexLHS = Array[];
for j=1:p[1]
    append!(edgeArrayIndexLHS, [edgeArrayGenList(j, "LHS")])
end;

# set the initial guess
u0 = [initConds[:, 1]; initConds[:, 2]; initPairConds[:, 1]; initPairConds[:, 2]];

# differential equation generator
function diffEq(du, u, p, t)

    # extract each of the states from the initial guess
    s = u[1 : p[1]];
    i = u[(p[1] + 1) : 2*p[1]];
    ss = u[(2*p[1] + 1) : (2*p[1] + p[2])];
    si = u[(2*p[1] + p[2] + 1) : end];

    # for each node, calculate the number of times the node is on the LHS of a S-I pairing
    sumSI = zeros(p[1], 1);
    for k=1:p[1]
        sumSI[k] = sum(si[edgeArrayIndexLHS[k]])
    end

    # calculate the sums of S-I states, depending on whether the S node is the i or j node in the SiIj expression
    neighbourInfRHS = (sumSI[last.(edgeArray)] - si) ./ s[last.(edgeArray)];
    neighbourInfLHS = (sumSI[first.(edgeArray)] - si) ./ s[first.(edgeArray)]; 

    # test if the S value in the above sums are different than 0 by testing whether the result is finite
    for i=1:p[2]
        neighbourInfRHS[i] = isfinite(neighbourInfRHS[i]) ? neighbourInfRHS[i] : 0
        neighbourInfLHS[i] = isfinite(neighbourInfLHS[i]) ? neighbourInfLHS[i] : 0
    end

    # set the differential equations associated with the S, I, S-S and S-I states
    du[1 : p[1]] .= -lambda * sumSI;
    du[(p[1] + 1) : 2*p[1]] .= lambda * sumSI .- gamma*i

    du[(2*p[1] + 1) : (2*p[1] + p[2])] = - ss.*neighbourInfLHS - ss.*neighbourInfRHS
    du[(2*p[1] + p[2] + 1) : end] = ss.*neighbourInfRHS - si.*neighbourInfLHS - lambda*si - gamma*si

end;

# set the ODE problem and solve it
prob = ODEProblem(diffEq, u0, (0, maxTime), p, dt = 0.01, saveat=t);
sol = solve(prob);

# extract and plot the solutions of the individual states
sSol = sol[1:p[1], :]
iSol = sol[(p[1] + 1):2*p[1], :]
rSol = 1 .- sSol .- iSol

plot(sSol', label="")
plot(iSol', label="")
plot(rSol', label="")

# extract and plot the solutions of the pair states
ssSol = sol[(2*p[1] + 1) : (2*p[1] + p[2]), :]
siSol = sol[(2*p[1] + p[2] + 1) : end, :]

plot(ssSol', label="")
plot(siSol', label="")



# s = u0[1 : p[1]]
# i = u0[(p[1] + 1) : 2*p[1]]
# ss = u0[(2*p[1] + 1) : (2*p[1] + p[2])]
# si = u0[(2*p[1] + p[2] + 1) : end]

# neighbourInfRHS = ss.*(sumSI[last.(edgeArray)] - si) ./ s[last.(edgeArray)]
# neighbourInfLHS = (sumSI[first.(edgeArray)] - si) ./ s[first.(edgeArray)]


