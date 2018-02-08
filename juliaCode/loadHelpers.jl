
#generateClifford actually allows gate dep errors to be applied, but not doing just now.

superCliffords = Array{Float64,2}[round.(generateClifford(i,0,0,0),5) for i=1:24];
operatorCliffords = generateRawCliffords();

function makeSuper(u)
    return round.(real(liou2pauliliou(liou(u))),10)
end

# Just some convenience functions
operatorPaulis= Array{Complex{Float64},2}[[1 0;0 1], [0 1;1 0],[0 -im;im 0], [1 0;0 -1]];

superPaulis = Array{Float64,2}[makeSuper(i) for i in operatorPaulis];
pI=operatorPaulis[1]
pX=operatorPaulis[2]
pY=operatorPaulis[3]
pZ=operatorPaulis[4];

piBy8 = exp(im*π/8)*[exp(-im*π/8) 0;0 exp(im*π/8)]
superPiBy8 = makeSuper(piBy8)

display("text/markdown","Created:\n-  superCliffords\n-  operatorCliffords\n-  superPaulis\n-  operator Paulis\n-  pI,pX,pY,pZ\n- T gate as piBy8 and superPiBy8");
display("text/markdown","- **makeSuper(operator)**");

# Quick function, given x, checks if its in the cliffords.
# Uses the SuperClifford to ignore phase
function findClifford(x::Array{Complex{Float64},2})
    test = makeSuper(x)
    return findfirst(x->x==test,superCliffords)
end

function findClifford(test::Array{Float64,2})
        return findfirst(x->x==test,superCliffords)
end


display("text/markdown","- **findClifford(x)** (super and operator");

function checkFrame(x)
    sum = 0
    for i=1:length(x)
        for j=1:length(x)
            sum += abs(trace(x[i]'*x[j]))^4
        end
    end
    sum/(length(x)^2)
end

display("text/markdown","- **checkFrame(x)** - operator")

function findInThis(needle,haystack)
    tofind = round.(needle,10)
    return findfirst(x->round.(x,10)==tofind,haystack)
end


display("text/markdown","- **findInThis(needle,haystack)**");

gx=expm(im*π*operatorPaulis[2]/4)
Gx=makeSuper(gx)
gy=expm(im*π*operatorPaulis[3]/4)
Gy=makeSuper(gy)
mgx=expm(-im*π*operatorPaulis[2]/4)
mGx=makeSuper(mgx)
mgy=expm(-im*π*operatorPaulis[3]/4)
mGy=makeSuper(mgy)


# My supeoperators are not in the same order as the 'operator' cliffords 
# reorder them.
reorderedSCliffords = copy(operatorCliffords)
for i=1:24
    reorderedSCliffords[i]= operatorCliffords[findClifford(makeSuper(operatorCliffords[i]))]
end

for i=1:24
    assert(makeSuper(reorderedSCliffords[i])==superCliffords[i])
end
#print("Sanity check - Ok\n")

operatorCliffords=reorderedSCliffords;

#This is our set of Generators (note I included the I)
twoGenSet=[]
push!(twoGenSet,superCliffords[24])
push!(twoGenSet,Gx)
push!(twoGenSet,Gy)

# Not particularly inspiring way to find the minimum generators.

function findAGeneratorFor(clifford)
    for fifth=1:2
        for i=1:3
            for j=1:3
                for k=1:3
                    for z=2:3 #stops I being done with just I generators.
                        if (twoGenSet[fifth]*twoGenSet[i]*twoGenSet[j]*twoGenSet[k]*twoGenSet[z] == superCliffords[clifford])
                            return [fifth,i,j,k,z]
                        end
                    end
                end
            end
        end
    end
    return []
end

# Find the minimum length generators for the Cliffords (1=I 2=Gx 3=Gy)
minGensWithI=[findAGeneratorFor(i) for i=1:24]

# Might as well get rid of the 1s
minGens = [i[find(x->x!=1,i)] for i in minGensWithI]
# Note the empty vector in position 24 - that is the identity.


fourGen=[]
push!(fourGen,superCliffords[24])
push!(fourGen,Gx)
push!(fourGen,Gy)
push!(fourGen,mGx)
push!(fourGen,mGy)
operatorFourGen=[]
push!(operatorFourGen,operatorCliffords[24])
push!(operatorFourGen,gx)
push!(operatorFourGen,gy)
push!(operatorFourGen,mgx)
push!(operatorFourGen,mgy)


function findAGeneratorFor(clifford,gens)
    for fifth=1:length(gens)
        for i=1:length(gens)
            for j=1:length(gens)
                for k=1:length(gens)
                    for z=2:length(gens)
                        if (gens[fifth]*gens[i]*gens[j]*gens[k]*gens[z] == superCliffords[clifford])
                            return [fifth,i,j,k,z]
                        end
                    end
                end
            end
        end
    end
    return []
end
min4Gens = [i[find(x->x!=1,i)] for i in [findAGeneratorFor(clif,fourGen) for clif=1:24]];


function findAllGeneratorsFor(clifford)
    gens=[]
    for fifth=1:3
    for i=1:3
        for j=1:3
            for k=1:3
                for z=2:3
                        if (twoGenSet[fifth]*twoGenSet[i]*twoGenSet[j]*twoGenSet[k]*twoGenSet[z] == superCliffords[clifford])
                            push!(gens,[fifth,i,j,k,z])
                    end
                end
            end
        end
    end
    end
    return gens
end

# all Gens contains all possible generators (maximum of 5 for the cliffords)
# this is how they are distributed.
allGens=[length(findAllGeneratorsFor(i)) for i=1:24]

display("text/markdown","- minGens, min4Gens, allGens *twoGenSet, fourGen")

function genA2RotateClifford(x)
    generator = minGens[x]
    start = superCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*twoGenSet[i]
        end
    end
    return start
end
operatorTestSet=[]
push!(operatorTestSet,operatorCliffords[24])
push!(operatorTestSet,gx)
push!(operatorTestSet,gy)
function genA2RotateOperatorClifford(x)
    generator = minGens[x]
    start = operatorCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*operatorTestSet[i]
        end
    end
    return start
end

display("text/markdown","- **genA2RotateOperatorClifford(x)**")

function genAFaulty2RotateClifford(x,noiseRotation)
    generator = minGens[x]
    start = superCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*twoGenSet[i]*noiseRotation
        end
    end
    return start
end

display("text/markdown","- **genAFaulty2RotateClifford(x,noiseRotation-matrix)**")
display("text/markdown","- **genAFaulty2RotateOperatorClifford(x,noiseRotation-matrix)**")


#Note the type of the noise here, so I don't mess up and pass in the PTSM noise
function genAFaulty2RotateOperatorClifford(x,noiseRotation::Array{Complex{Float64},2})
    generator = minGens[x]
    start = operatorCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*noiseRotation*operatorTestSet[i]
        end
    end
    return start
end

function genAFaulty4RotateClifford(x,noiseRotation)
    generator = min4Gens[x]
    start = superCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*fourGen[i]*noiseRotation
        end
    end
    return start
end 

#Note the type of the noise here, so I don't mess up and pass in the PTSM noise

function genAFaulty4RotateOperatorClifford(x,noiseRotation::Array{Complex{Float64},2})
    generator = min4Gens[x]
    start = operatorCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*noiseRotation*operatorFourGen[i]*noiseRotation'
        end
    end
    return start
end

display("text/markdown","- **genAFaulty4RotateClifford(x,noiseRotation-matrix)**")
display("text/markdown","- **genAFaulty4RotateOperatorClifford(x,noiseRotation-matrix)**")


#This final one just adds the noise to the cliffords - i.e as per protocol assumptions.
function genAFaultyClifford(x,noiseRotation)
    return superCliffords[x]*noiseRotation
    
end


display("text/markdown","- **genAFaultyClifford(x,noiseRotation-matrix)**")

#Check this out with some example noise
θ=0.01
znoise=expm(-im*π*θ*[1 0;0 -1]/2)
Znoise = makeSuper(znoise)

function getFidelity(gate,clifford)
    actualχ = pauliliou2chi(gate*superCliffords[clifford]')[1,1]/4.0
    actualFidelity = (2*actualχ+1)/3
end

function twirl(sequence,victim)
    return sum([i*victim*i' for i in sequence])/length(sequence)
end

function extractFidelity(sc)
    return (trace(sc)-1)/3
end

display("text/markdown","- **getFidelity(gate,clifford)**\n-  **twirl(sequence,victim)**\n-  **extractFidelity(sc)**")


display("text/markdown","- **makeRuns(getNumbers,number,gateGenerator,noise,initialState,endState)**");

"""
    makeRuns(function: (number of gates)->seqence of gates, number of gates,
             function: (gate number)->gate to test,
             initial state (vectorized),
             endMeasurement (vectorized)')

Conducts an RB run starting in initial state, measuring in endMeasurement and using the gates supplied by the gateGenerator function.
Returns the probability of a +1 measurement.

"""
function makeRuns(getNumbers,number,gateGenerator,initialState,endState)
    gates = getNumbers(number)
    # state keeps track of the theoretical state so we know what the inversion clifford is
    state = superCliffords[gates[1]]
    # state E is the noisy state.
    stateE = gateGenerator(gates[1])
    for gate = 2:length(gates)
        a = superCliffords[gates[gate]]
        ae = gateGenerator(gates[gate])
        state = a*state
        stateE = ae*stateE
    end
    # state' is the inverse SuperOperator of the state 
    # findClifford(state') returns the clifford number of the inverse so we can create a noisy version.
    return ((endState*gateGenerator(findClifford(state'))*stateE*initialState)[1]+1)/2
end

function getACliffordSequence(number)
    return rand(1:24,1,number)
end

initialZState=1/sqrt(2)*[1;0;0;1]
measurementZ =sqrt(2)*[0 0 0 1]

#Benchmark model
model(x, p) = p[1]*(p[2].^x)+p[3]

display("text/markdown","- **doRuns(theGeneratorFunction,)**")

# Single clifford seqence - the GeneratorFunction is a function that takes 
# a clifford number with noise and returns the noisy gate.
"""
    doRuns(function GateGeneratorFunction(x)->Gate)

# Arguments
* 'theGeneratorFunction::(x)->Gate
* 'sequences::Array{Int64}'=1:10:300
* 'repetitions::Int'=300

returns the calculated p values(collection)

"""
function doRuns(theGeneratorFunction,sequences=1:10:300,repetitions=300)
    measurements = repetitions
    results= Float64[mean(
        [makeRuns(
                    getACliffordSequence,j,theGeneratorFunction,initialZState,measurementZ) 
                    for i =1:measurements]
                 )[1] for j in sequences]
    fit = curve_fit(model, sequences,results, [0.5,0.9, 0.5]);
   # print("P is $(fit.param[2])\n")
    return fit.param[2]
end

