# getSomeCliffords.jl

#Creates a bunch of 1 and 2 qubit Cliffords (complex and real) for ease of use.
#Defines a number of obvious operators
try 
    display("text/markdown","Creating a few things, this may take a short (but noticeable) time\n ")
catch
    print("Creating a few things, this may take a short (but noticeable) time\n")
end


#Pkg.update()
include("Initial.jl")
include("Symplectic.jl")
include("open-systems.jl")
include("PTSM.jl")
include("loadHelpers.jl")

using PyCall
using ProgressMeter
using DelimitedFiles

# Generate the full set of 2 qubit cliffords for comparison purposes
FullCliffords= Array{Complex{Float64},2}[]
FullCliffordCommands = Array{String,1}[]
FullCliffordsD=Array{Int64,1}[]
n=11520
for i=0:719 # number of symplectics
    for j=0:15 # number of phase variations.
        global state
        state=setup(2)
        decompose(i,j,2,true,false) # symplectic i, bits j, 2 qubits, no output and rationalise
        push!(FullCliffordsD,[i,j])
        push!(FullCliffordCommands,commands)
        push!(FullCliffords,makeFromCommand(commands))
    end
end
SuperCliffs = [makeSuper(s) for s in FullCliffords];

try 
    display("text/markdown","1 of 3, all two qubit cliffords and super cliffords generated.")
catch 
    print("1 of 3, all two qubit cliffords and super cliffords generated.\n")
end

#Some basic gates

pI=[1 0;0 1]
pX=[0 1;1 0]
pZ=[1 0;0 -1]
pP=[1 0;0 im]
pH=[1 1;1 -1]/sqrt(2)
cnot12 = kron([1 0]'*[1 0],pI)+kron([0 1]'*[0 1],pX)


xby2 = exp(-im*π*pX/4) # same as sqrt(X)
cnot21 = kron(pI,[1 0]'*[1 0])+kron(pX,[0 1]'*[0 1])

pXpI = kron(pX,pI)
pIpX = kron(pI,pX)
pZpI = kron(pZ,pI)
pIpZ = kron(pI,pZ)
pHpI = kron(pH,pI)
pIpH = kron(pI,pH)
pPpI = kron(pP,pI)
pIpP = kron(pI,pP)

twoQubitGens = [pXpI,pIpX,pZpI,pIpZ,pHpI,pIpH,cnot12,cnot21];

#Index to the cliffords
singleReals=[]
doubleReals=[]

# Just for the sake of it find the single qubit Real clifford group
# The generators we will use here are pZ and PH
push!(singleReals,findClifford(makeSuper(pZ)))
push!(singleReals,findClifford(makeSuper(pH)))

# Note don't use the isequal predicate in findfirst as we need 0.0 to equal -0.0

let doneOne = true
    while doneOne
        doneOne=false
        toIterate = copy(singleReals)
        for i in toIterate
            n1 = findClifford(makeSuper(pZ)*superCliffords[i])
            if (findfirst(x->x==n1,singleReals) == nothing)
                push!(singleReals,n1)
                doneOne=true
            end
            n2 = findClifford(makeSuper(pH)*superCliffords[i])
            if (findfirst(x->x==n2,singleReals) == nothing)
                push!(singleReals,n2)
                doneOne=true
            end
        end
    end
end    



#Check the frame reference, this is an orthogonal 2-design, should be 3

@assert round(checkFrame([operatorCliffords[x] for x in singleReals]),digits = 8)==3

# Then generate using the generators for the real cliffords
# Note here we have both X and Z generators, they are not a minimal set
superReal2Gens = [makeSuper(x) for x in twoQubitGens];

doubleReals = []
for i in superReal2Gens
    push!(doubleReals,findfirst(x->x==i,SuperCliffs))
end

let doneOne = true
    while doneOne
        doneOne=false
        toIterate = copy(doubleReals)
        for i in toIterate
            for j in superReal2Gens
                lookFor = j*SuperCliffs[i]
                n1 = findfirst(x->x == lookFor,SuperCliffs)
                if (findfirst(x->x == n1,doubleReals) == nothing)
                    push!(doubleReals,n1)
                    doneOne=true
                end
            end
        end
    end
end


#Check the frame reference, this is an orthogonal 2-design, should be 3

@assert round(checkFrame([FullCliffords[x] for x in doubleReals]),digits=8)==3

try 
    display("text/markdown","2 of 3, all two qubit real cliffords indexed in.")
catch
    print("2 of 3, all two qubit real cliffords indexed in.\n")
end

# This bit is directory dependent - should think about it

pushfirst!(PyVector(pyimport("sys")["path"]), "./finite-fields/")
@pyimport finitefield as ff

F22 = ff.FiniteField(2,2)

@pyimport SLFunctions as sl


function getStabilisers()
    jstabs=Array{Int64,2}[]
    stabs = sl.getF2Stabilisers()
    for i in 1:60
        push!(jstabs,reshape(stabs[i,:,:],4,4))
    end
    return jstabs
end

stabs=getStabilisers();


#since we saved them previously we can load them
gates=Array{Float64}[]
for i=1:60
    push!(gates,readdlm("./compiledGates/gateList$i.csv"))
end

gates[29]=[7.0] # This is the identity gate, it was saved wrong


gatesCommand = (Expr)[]
	svec = []
	for i = 1:2
	  for j = 1:2
	  	if i!=j
            push!(gatesCommand,Expr(:call,:cnot,:svec,i,j,false))
	  	end
	  end
    push!(gatesCommand,Expr(:call,:hadamard,:svec,i,false))
    push!(gatesCommand,Expr(:call,:phase,:svec,i,false))
	end
push!(gatesCommand,Expr(:call,:noop))


#Make them, not sure why I don't just save them!

sl2Cliffords= Array{Complex{Float64},2}[]
sl2CliffordCommands = Array{String,1}[]
n=60#*4*4
for i=1:60
    for j=1:4
        for z=1:4
            global state
            state=setup(2)
            state[:,1:4]=stabs[i]
            decomposeState(state,true)
            #Apply the PAULIS!
            if (j==2)
                push!(commands,"hadamard(1)")
                push!(commands,"phase(1)")
                push!(commands,"phase(1)")
                push!(commands,"hadamard(1)")
            elseif (j==3)
                push!(commands,"phase(1)")
                push!(commands,"phase(1)")
            elseif (j==4)
                push!(commands,"phase(1)")
                push!(commands,"phase(1)")
                push!(commands,"hadamard(1)")
                push!(commands,"phase(1)")
                push!(commands,"phase(1)")
                push!(commands,"hadamard(1)")
            end
            if (z==2)
                push!(commands,"hadamard(2)")
                push!(commands,"phase(2)")
                push!(commands,"phase(2)")
                push!(commands,"hadamard(2)")
            elseif (z==3)
                push!(commands,"phase(2)")
                push!(commands,"phase(2)")
            elseif (z==4)
                push!(commands,"phase(2)")
                push!(commands,"phase(2)")
                push!(commands,"hadamard(2)")
                push!(commands,"phase(2)")
                push!(commands,"phase(2)")
                push!(commands,"hadamard(2)")
            end
            push!(sl2CliffordCommands,commands)
            push!(sl2Cliffords,makeFromCommand(commands))
       end
    end
end

try
    display("text/markdown","3 of 3, all Cleve cliffords found.")
catch
    print("3 of 3, all Cleve cliffords found.\n")
end


sl60_2Cliffords= Array{Complex{Float64},2}[]
sl60_2CliffordCommands = Array{String,1}[]
n=60#*4*4
for i=1:60
        global state
       state=setup(2)
       state[:,1:4]=stabs[i]
       decomposeState(state,true)
       push!(sl60_2CliffordCommands,commands)
       push!(sl60_2Cliffords,makeFromCommand(commands))
end

super60_sl2=[makeSuper(i) for i in sl60_2Cliffords];


supersl2=[makeSuper(i) for i in sl2Cliffords];

try
    display("text/markdown","Created:\n   ")
    display("text/markdown","- **FullCliffords** - all two qubit cliffords")
    display("text/markdown","- **SuperCliffs** - two qubit cliffords as superOperators in the Pauli basis")
    display("text/markdown","- **sl60_2Cliffords** - the cleve 60 2 qubits and **super60_sl2**")
    display("text/markdown","- **sl2Cliffords** - the cleve 60 twirled by 16 paulis and **supersl2**")
    display("text/markdown","- **singleReals** and **doubleReals** the index of the real cliffords into the full clifford arrays")
catch 
    print("Created:\n")
    print("\t- **FullCliffords** - all two qubit cliffords\n")
    print("\t- **SuperCliffs** - two qubit cliffords as superOperators in the Pauli basis\n")
    print("\t- **sl60_2Cliffords** - the cleve 60 2 qubits and **super60_sl2**\n")
    print("\t- **sl2Cliffords** - the cleve 60 twirled by 16 paulis and **supersl2**\n")
    print("\t- **singleReals** and **doubleReals** the index of the real cliffords into the full clifford arrays\n")
end
