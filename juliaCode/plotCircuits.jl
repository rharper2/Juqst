using PyCall
using PyPlot
using Juqst
pushfirst!(PyVector(pyimport("sys")."path"), "./PlotQCircuit/")
pushfirst!(PyVector(pyimport("sys")."path"), "../PlotQCircuit/")
pqc = pyimport("plot_quantum_circuit")


function scheduleIt(t::Tableau)
    (li,la)=drawCircuit(t)
    pqc.plot_quantum_schedule(pqc.make_schedule(li),labels=la)
end

function plotIt(t::Tableau)
    (li,la)=drawCircuit(t)
    pqc.plot_quantum_circuit(li,labels=la)
end
    
