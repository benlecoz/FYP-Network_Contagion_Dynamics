using GraphPlot, Graphs

function generateAdj(numNodes, graphType, graphParams)
    if graphType == "complete"

        Adj = ones(numNodes, numNodes) - I(numNodes)
        
        evenlySpaced = transpose(LinRange(0, 1, numNodes + 1))
              
        G = Graphs.Graph(Adj)

        locs_x = vec(cos.(2*pi .* evenlySpaced))
        locs_y = vec(sin.(2*pi .* evenlySpaced))
        
        GraphPlot.gplot(G, locs_x, locs_y)

    end

    return G, G_com

end

G, G_com = generateAdj(8, "complete", []);