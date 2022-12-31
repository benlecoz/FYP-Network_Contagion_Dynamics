include("./runSirGillespie.jl")
using .runGillespieSIR
using Statistics 
using Plots

probS, probI, probR, probSS, probSI, probSR, probIS, probII, probIR, probRS, probRI, probRR, t = runGillespieSIR.RunSIRGillespie();

function plotAvgIndivStates()

    # Plot averaged evolution of individual states throughout the simulations

    meanS = mean(probS, dims=1)
    meanI = mean(probI, dims=1)
    meanR = mean(probR, dims=1)

    plot(vec(t).*100, vec(meanS).*100, label="Susceptible", xlabel="τ-leaps", ylabel="Percentage of individuals in each state (%)", legend=:right)
    plot!(vec(t).*100, vec(meanI).*100, label="Infected")
    plot!(vec(t).*100, vec(meanR).*100, label="Recovered")
end

function plotAvgPairStates()

    # Plot 3x3 grid of the averaged evolution of pairs of states during the simulations

    meanSS = mean(probSS, dims=1);
    meanSI = mean(probSI, dims=1);
    meanSR = mean(probSR, dims=1);
    meanIS = mean(probIS, dims=1);
    meanII = mean(probII, dims=1);
    meanIR = mean(probIR, dims=1);
    meanRS = mean(probRS, dims=1);
    meanRI = mean(probRI, dims=1);
    meanRR = mean(probRR, dims=1);

    l = @layout [a b c; d e f; g h i];
    plt1 = plot(vec(t).*100, vec(meanSS).*100, label="S-S", legend=:topright, ylims=(0,100), linecolor="brown");
    plt2 = plot(vec(t).*100, vec(meanSI).*100, label="S-I", ylims=(0,100), linecolor="red");
    plt3 = plot(vec(t).*100, vec(meanSR).*100, label="S-R", ylims=(0,100), linecolor="orange");
    plt4 = plot(vec(t).*100, vec(meanIS).*100, label="I-S", ylims=(0,100), linecolor="green", ylabel="Averaged evolution of pair states (%)");
    plt5 = plot(vec(t).*100, vec(meanII).*100, label="I-I", ylims=(0,100), linecolor="cyan");
    plt6 = plot(vec(t).*100, vec(meanIR).*100, label="I-R", ylims=(0,100), linecolor="blue");
    plt7 = plot(vec(t).*100, vec(meanRS).*100, label="R-S", ylims=(0,100), linecolor="magenta");
    plt8 = plot(vec(t).*100, vec(meanRI).*100, label="R-I", ylims=(0,100), linecolor="grey", xlabel="τ-leaps");
    plt9 = plot(vec(t).*100, vec(meanRR).*100, label="R-R", ylims=(0,100), linecolor="black");

    pairPlot = plot(plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8, plt9, layout = l)

    display(pairPlot)
end

function plotAvgNodeStates(state)

    # Plot each nodes average probability evolution for a given state

    stateTranslation = Dict("S" => probS, "I" => probI, "R" => probR, "SS" => probSS, "SI" => probSI, "SR" => probSR, "IS" => probIS, "II" => probII, "IR" => probIR, "RS" => probRS, "RI" => probRI, "RR" => probRR)
    plot(stateTranslation[state]'.*100, label="", xlabel="τ-leaps", ylabel="Probability of each node (%)")
end

function graphEvo()

    # Plot four networds showing a color-coded epidemic spread from start to end
    # Requires sirGillespie code to be changed to return matrices of an individual run (called indivRunS, indivRunI, indivRunR), and the graph variable

    function colour_nodes(state)

        membershipStart = zeros(Int64, numNodes);
        membershipMid1 = zeros(Int64, numNodes);
        membershipMid2 = zeros(Int64, numNodes);
        membershipEnd = zeros(Int64, numNodes);

        if state == "start"
            vector = membershipStart
            time_step = 1
        elseif state == "mid1"
            vector = membershipMid1
            time_step = 150
        elseif state == "mid2"
            vector = membershipMid2
            time_step = 300
        elseif state == "end"
            vector = membershipEnd
            time_step = 3000
        end

        for i = 1:numNodes
            if indivRunS[i, time_step] == 1
                vector[i, 1] = 1
            elseif indivRunI[i, time_step] == 1
                vector[i, 1] = 2
            elseif indivRunR[i, time_step] == 1
                vector[i, 1] = 3
            end
        end
        
        return vector
    end

    start = colour_nodes("start")
    mid1 = colour_nodes("mid1")
    mid2 = colour_nodes("mid2")
    ending = colour_nodes("end")

    nodecolour = [colorant"lightblue", colorant"red", colorant"green"];
    edgecolour = [colorant"gray"];

    start_plot = gplot(graph, layout=circular_layout, nodefillc = nodecolour[start], edgestrokec=[colorant"gray"])
    gplot(graph, layout=circular_layout, nodefillc = nodecolour[mid1], edgestrokec=[colorant"gray"])
    gplot(graph, layout=circular_layout, nodefillc = nodecolour[mid2], edgestrokec=[colorant"gray"])
    gplot(graph, layout=circular_layout, nodefillc = nodecolour[ending], edgestrokec=[colorant"gray"])

    # savefig("C:\\Users\\benle\\GitHub\\FYP-Network_Contagion_Dynamics\\epi_start.png")

end