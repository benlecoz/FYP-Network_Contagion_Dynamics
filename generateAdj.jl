module AdjGenerator

using GraphPlot, Graphs

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
        n = trunc(Int, floor(numNodes/2))
        graph = barbell_graph(n, n)
    
    elseif graphType == "binaryTree"
        graph = binary_tree(numNodes)
    
    elseif graphType == "ErdosRenyi"
        graph = erdos_renyi(numNodes, graphParams[1])

    elseif graphType == "ring"
        graph = cycle_graph(numNodes)
        plot = "circular plotting"
    
    end

    if plot == "circular plotting"
        evenlySpaced = transpose(LinRange(0, 1, numNodes + 1))

        locs_x = vec(cos.(2*pi .* evenlySpaced))
        locs_y = vec(sin.(2*pi .* evenlySpaced))
        graph_plot = gplot(graph, locs_x, locs_y)
    
    else
        graph_plot = gplot(graph)
    
    end
    
    # display(graph_plot)

    return graph
   
end

end