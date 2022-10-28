using GraphPlot, Graphs

function generateAdj(numNodes, graphType, graphParams, graphOptions)
    if graphType == "complete"
        graph = complete_graph(numNodes)

        evenlySpaced = transpose(LinRange(0, 1, numNodes + 1))

        locs_x = vec(cos.(2*pi .* evenlySpaced))
        locs_y = vec(sin.(2*pi .* evenlySpaced))

    elseif graphType == "WattsStrogatz"
        graph = watts_strogatz(numNodes, trunc(Int, graphParams[1]), graphParams[2])
    elseif graphType == "BarabasiAlbert"
        graph = barabasi_albert(numNodes, graphParams[1])
    elseif graphType == "barbell"
        n = trunc(Int, floor(numNodes/2))
        graph = barbell_graph(n, n)
    elseif graphType == "binaryTree"
        graph = binary_tree(numNodes)

    end

    if isnothing(graphOptions) 
    elseif graphOptions == "plot1"
        graph_plot = gplot(graph, locs_x, locs_y)
        graph_plot
    elseif graphOptions == "plot2"
        graph_plot = gplot(graph)
        graph_plot
    else
        savefig(pwd() * "\\" * graphOptions * ".png")
    end
    
    return graph, graph_plot
   
end

generateAdj(4, "binaryTree", [], "plot2")[2]

