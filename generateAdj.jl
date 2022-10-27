using GraphPlot, Graphs, LightGraphs


function generateAdj(numNodes, graphType, graphParams)
    if graphType == "complete"

        Adj = ones(numNodes, numNodes) - I(numNodes)
        
        evenlySpaced = transpose(LinRange(0, 1, numNodes + 1))
        
        G_com = LightGraphs.complete_graph(numNodes)

        locs_x = vec(cos.(2*pi .* evenlySpaced))
        locs_y = vec(sin.(2*pi .* evenlySpaced))
        
        GraphPlot.gplot(G_com)

    end

    return G, G_com

end

G, G_com = generateAdj(8, "complete", []);