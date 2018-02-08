


# My supeoperators are not in the same order as the 'straight' cliffords 
# reorder them.
straightCliffords = generateRawCliffords();
reorderedSCliffords = copy(straightCliffords)
for i=1:24
    reorderedSCliffords[i]= straightCliffords[findInThis(makeSuper(straightCliffords[i]),superCliffords)]
end
straightCliffords=reorderedSCliffords;

print("Defined straightCliffords (2x2 cliffords)\n")

gx=expm(im*π*straightPaulis[2]/4)
Gx=makeSuper(gx)
gy=expm(im*π*straightPaulis[3]/4)
Gy=makeSuper(gy)
mgx=expm(-im*π*straightPaulis[2]/4)
mGx=makeSuper(mgx)
mgy=expm(-im*π*straightPaulis[3]/4)
mGy=makeSuper(mgy)
#This is our set of Generators (note I included the I)
testSet=[]
push!(testSet,superCliffords[24])
push!(testSet,Gx)
push!(testSet,Gy)

function findAGeneratorFor(clifford)
    for fifth=1:3
        for i=1:3
            for j=1:3
                for k=1:3
                    for z=1:3
                        if (testSet[fifth]*testSet[i]*testSet[j]*testSet[k]*testSet[z] == superCliffords[clifford])
                            return [fifth,i,j,k,z]
                        end
                    end
                end
            end
        end
    end
    return []
end

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
straightFourGen=[]
push!(straightFourGen,straightCliffords[24])
push!(straightFourGen,gx)
push!(straightFourGen,gy)
push!(straightFourGen,mgx)
push!(straightFourGen,mgy)


function findAGeneratorFor(clifford,gens)
    for fifth=1:length(gens)
        for i=1:length(gens)
            for j=1:length(gens)
                for k=1:length(gens)
                    for z=1:length(gens)
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

#Check it works

function genA2RotateClifford(x)
    generator = minGens[x]
    start = superCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*testSet[i]
        end
    end
    return start
end
straightTestSet=[]
push!(straightTestSet,straightCliffords[24])
push!(straightTestSet,gx)
push!(straightTestSet,gy)
function genA2RotateStraightClifford(x)
    generator = minGens[x]
    start = straightCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*straightTestSet[i]
        end
    end
    return start
end

#Check it works
for i=1:24
    assert(genA2RotateClifford(i)==superCliffords[i])
end
print("genA2RotateClifford(Clifford_number) defined and working!\n")



function genAFaulty2RotateClifford(x,noiseRotation)
    generator = minGens[x]
    start = superCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*testSet[i]*noiseRotation
        end
    end
    return start
end

#Note the type of the noise here, so I don't mess up and pass in the PTSM noise
function genAFaulty2RotateStraightClifford(x,noiseRotation::Array{Complex{Float64},2})
    generator = minGens[x]
    start = straightCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*straightTestSet[i]*noiseRotation
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

function genAFaulty4RotateStraightClifford(x,noiseRotation::Array{Complex{Float64},2})
    generator = min4Gens[x]
    start = straightCliffords[24] # Identity
    for i in generator
        if i!=1 # Identity
            start = start*straightFourGen[i]*noiseRotation
        end
    end
    return start
end

#This final one just adds the noise to the cliffords - i.e as per protocol assumptions.
function genAFaultyClifford(x,noiseRotation)
    return superCliffords[x]*noiseRotation
    
end


#Check it works
for i=1:24
    assert(genAFaulty4RotateClifford(i,superCliffords[24])==superCliffords[i])
end
print("Defined: genAFaulty2RotateClifford, genAFaulty4RotateClifford, genAFaultyClifford")





