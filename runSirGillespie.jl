module runGillespieSIR

include("./sirGillespie.jl")
using .GillespieSIR
using Statistics 
using Plots

function RunSIRGillespie()

    # generate Adjacency matrix
    numNodes = 50;

    graphType = "WattsStrogatz";
    graphParams = [4, 0.5];

    # define contagion dynamics parameters
    lambda = 1;
    gamma = 0.1;
    modelParams = [lambda, gamma];

    # time parameters
    maxTime = 12.0;
    timeRes = 0.01;

    # initial conditions with one infected node
    initConds = zeros(Int64, numNodes, 1);
    initConds[1] = 1;

    # Gillespie model parameters
    numRuns = 10^3;

    probS, probI, probR, probSS, probSI, probSR, probIS, probII, probIR, probRS, probRI, probRR, t, edgeArray = GillespieSIR.SIRGillespie(numNodes, graphType, graphParams, modelParams, maxTime, timeRes, initConds, numRuns)

    return probS, probI, probR, probSS, probSI, probSR, probIS, probII, probIR, probRS, probRI, probRR, t, edgeArray

end

end