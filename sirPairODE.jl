module sirPairODE

include("./generateAdj.jl")
using .AdjGenerator: edgeArrayGenList 
using Graphs, SparseArrays, DifferentialEquations
using Plots, InvertedIndices, BenchmarkTools

function SIRPairODE(edgeArray::Array, graph::SimpleGraph, modelParams::Array, maxTime::Float64, timeRes::Float64, initInfNodes::Array)

    # extract model parameters and the number of nodes from the inputs
    lambda, gamma = modelParams;
    numNodes = nv(graph);
    
    # edgeArray only has directed edges, so create new array of edges with oppositely directed edges for every node pair
    edgeArrayNew = copy(edgeArray);
    for (i, t) in pairs(edgeArray)
        push!(edgeArrayNew, reverse(t))
    end;

    # calculate the number of edges in the initial and the extended arrays of edges
    numEdges = length(edgeArrayNew);
    numEdgesTrue = length(edgeArray);
    
    # set the initial conditions of each nodes
    initConds = zeros(Int64, numNodes, 2);
    initConds[initInfNodes, 2] .= 1;
    initConds[Not(initInfNodes), 1] .= 1;

    # extract time components
    maxTime = timeRes*ceil(maxTime/timeRes);
    t = Array(collect(Float64, 0:timeRes:maxTime)');

    # ensure that edges containing two initially infected nodes are not accidentally marked as S-I 
    initIIEdges = zeros(Int64, 0);
    for i in AdjGenerator.edgeArrayGenList(edgeArrayNew, initInfNodes, "RHS")
        if in.(first.(edgeArrayNew)[i], Ref(initInfNodes))
            append!(initIIEdges, i)
        end
    end
    
    # filter the I-I edges from the S-I and I-S edges
    initSIEdges = filter!(x->!(x in initIIEdges), AdjGenerator.edgeArrayGenList(edgeArrayNew, initInfNodes, "RHS"));

    # based on initially infected nodes and edgeArray, set the initial pair conditions of the S-S and S-I states
    initPairConds = zeros(Int64, numEdges, 2);
    initPairConds[initSIEdges, 2] .= 1;
    initPairConds[Not(edgeArrayGenList(edgeArrayNew, initInfNodes, "both")), 1] .= 1;

    # set the parameters associated with the ODEProblem
    p = [lambda, gamma];

    # for each node, find the edge pairs for which the node is on the LHS (will be used to find which nodes are in the S-I state)  
    edgeArrayIndexLHS = Array[];
    for j=1:numNodes
        append!(edgeArrayIndexLHS, [edgeArrayGenList(edgeArrayNew, j, "LHS")])
    end;

    # set the initial guess
    u0 = [initConds[:, 1]; initConds[:, 2]; initPairConds[:, 1]; initPairConds[:, 2]];

    # differential equation generator
    function diffEq(du, u, p, t)

        lambda, gamma = p

        # extract each of the states from the initial guess
        s = u[1 : numNodes];
        i = u[(numNodes + 1) : 2*numNodes];
        ss = u[(2*numNodes + 1) : (2*numNodes + numEdges)];
        si = u[(2*numNodes + numEdges + 1) : end];

        # since the edgeArray is mirrored at the halfway point, the I-S array is a pivot of the S-I array at the halfway point
        is = [si[numEdgesTrue+1:end]; si[1:numEdgesTrue]];

        # for each node, calculate the number of times the node is on the LHS of a S-I pairing
        sumSI = zeros(numNodes, 1);
        for j=1:numNodes
            sumSI[j] = sum(lambda*si[edgeArrayIndexLHS[j]])
        end

        # calculate the sums of S-I states, depending on whether the susceptible node is the i or j node in the SiIj expression
        neighbourInfRHS = (sumSI[last.(edgeArrayNew)] - lambda*is) ./ s[last.(edgeArrayNew)];
        neighbourInfLHS = (sumSI[first.(edgeArrayNew)] - lambda*si) ./ s[first.(edgeArrayNew)]; 

        # test if the S value in the above sums are non-zero by testing whether the result is finite
        for k=1:numEdges
            neighbourInfRHS[k] = isfinite(neighbourInfRHS[k]) ? neighbourInfRHS[k] : 0
            neighbourInfLHS[k] = isfinite(neighbourInfLHS[k]) ? neighbourInfLHS[k] : 0
        end

        # set the differential equations associated with the S, I, S-S and S-I states
        du[1 : numNodes] .= -sumSI;
        du[(numNodes + 1) : 2*numNodes] .= sumSI .- gamma*i

        du[(2*numNodes + 1) : (2*numNodes + numEdges)] = - ss.*neighbourInfLHS - ss.*neighbourInfRHS
        du[(2*numNodes + numEdges + 1) : end] = ss.*neighbourInfRHS - si.*neighbourInfLHS - lambda*si - gamma*si

    end;

    # set the ODE problem and solve it
    prob = ODEProblem(diffEq, u0, (0, maxTime), p, dt = timeRes, saveat=t);
    sol = solve(prob);

    # extract the solutions of the individual states
    sSol = sol[1:numNodes, :];
    iSol = sol[(numNodes + 1):2*numNodes, :];
    rSol = 1 .- sSol .- iSol;

    # extract the solutions of the pair states
    ssSol = sol[(2*numNodes + 1) : (2*numNodes + numEdges), :];
    siSol = sol[(2*numNodes + numEdges + 1) : end, :];

    return sSol, iSol, rSol, ssSol, siSol

end

end