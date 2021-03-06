using PyCall
using PyPlot
unshift!(PyVector(pyimport("sys")["path"]), "./PlotQCircuit/")
unshift!(PyVector(pyimport("sys")["path"]), "../PlotQCircuit/")
@pyimport plot_quantum_circuit as pqc


function scheduleIt()
    (li,la)=drawCircuit()
    pqc.plot_quantum_schedule(pqc.make_schedule(li),labels=la)
end

function plotIt()
    (li,la)=drawCircuit()
    pqc.plot_quantum_circuit(li,labels=la)
end
    
