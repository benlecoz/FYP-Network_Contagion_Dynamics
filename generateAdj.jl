using GraphPlot, Graphs

function generateAdj(numNodes, graphType, graphParams, graphOptions)
    if graphType == "complete"
        graph = complete_graph(numNodes)
    end

    evenlySpaced = transpose(LinRange(0, 1, numNodes + 1))

    locs_x = vec(cos.(2*pi .* evenlySpaced))
    locs_y = vec(sin.(2*pi .* evenlySpaced))

    if graphOptions == "plot"
        graph_plot = gplot(graph, locs_x, locs_y)
    else
        savefig(pwd() * "\\" * graphOptions * ".png")
    end
    
    return graph, graph_plot
   
end

graph, graph_plot = generateAdj(10, "complete", [], "plot")

graph_plot

