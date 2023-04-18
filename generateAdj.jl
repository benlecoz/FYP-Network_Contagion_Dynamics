module AdjGenerator

using Graphs

function generateAdj(numNodes, graphType, graphParams)

    plot = "normal"

    if graphType == "complete"
        graph = complete_graph(numNodes)
        plot = "circular plotting"

    elseif graphType == "WattsStrogatz"
        graph = watts_strogatz(numNodes, trunc(Int, graphParams[1]), graphParams[2])
        plot = "circular plotting"

    elseif graphType == "BarabasiAlbert"
        graph = barabasi_albert(numNodes, graphParams[1])
    
    elseif graphType == "barbell"
        n = trunc(Int, floor(numNodes/2));
        graph = barbell_graph(n, n)
    
    elseif graphType == "ErdosRenyi"
        graph = erdos_renyi(numNodes, graphParams[1])

    elseif graphType == "ring"
        graph = cycle_graph(numNodes)
        plot = "circular plotting"
    
    elseif graphType == "tree"
        graph = SimpleGraph{Int64}(numNodes)
        depth = zeros(Int64, numNodes)

        for i in 2:numNodes
            parent = rand(1:i-1)
            add_edge!(graph, parent, i)

            if i==2
                depth[i] = 1
            else
                depth[i] = depth[parent] + 1
            end
        end

        return graph, depth

    elseif graphType == "chain"
        graph = SimpleGraph{Int64}(numNodes)

        for i in 2:numNodes
            parent = i-1
            add_edge!(graph, parent, i)
        end

    end

    return graph
   
end

# Find edgeArray index of edges for node(s) in nodeList
# Function allows input on whether the node in question is on the LHS, RHS or both sides of the edge pair 
function edgeArrayGenList(edgeArray, nodeList, side)
    if side == "both"
        edgesNodeList = unique([reduce(vcat, [findall(x->x==i, first.(edgeArray)) for i in nodeList])..., (reduce(vcat, [findall(x->x==i, last.(edgeArray)) for i in nodeList])...)])
    elseif side == "LHS"
        edgesNodeList = reduce(vcat, [findall(x->x==i, first.(edgeArray)) for i in nodeList])
    elseif side == "RHS"
        edgesNodeList = reduce(vcat, [findall(x->x==i, last.(edgeArray)) for i in nodeList])
    end

    return edgesNodeList
end

end