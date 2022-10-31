using GraphPlot, Graphs

function generateAdj(numNodes, graphType, graphParams, graphOptions)

    plot = "normal"

    if graphType == "complete"
        graph = complete_graph(numNodes)
        plot = "circular plotting"

    elseif graphType == "WattsStrogatz"
        graph = watts_strogatz(numNodes, trunc(Int, graphParams[1]), graphParams[2])

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

    if isnothing(graphOptions) 

    elseif graphOptions == "plot"

        if plot == "circular plotting"
            evenlySpaced = transpose(LinRange(0, 1, numNodes + 1))

            locs_x = vec(cos.(2*pi .* evenlySpaced))
            locs_y = vec(sin.(2*pi .* evenlySpaced))
            graph_plot = gplot(graph, locs_x, locs_y)
        
        else
            graph_plot = gplot(graph)
        
        end

    else
        savefig(pwd() * "\\" * graphOptions * ".png")

    end
    
    return graph, graph_plot
   
end

generateAdj(12, "ErdosRenyi", [0.15], "plot")[2]

Graphs.cycle_graph()