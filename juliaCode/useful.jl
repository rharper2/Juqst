
superCliffords = Array{Float64,2}[round.(generateClifford(i,0,0,0),5) for i=1:24];
straightCliffords = generateRawCliffords();

function makeSuper(u)
    return round.(real(liou2pauliliou(liou(u))),10)
end

# Just some convenience functions
straightPaulis= Array{Complex{Float64},2}[[1 0;0 1], [0 1;1 0],[0 -im;im 0], [1 0;0 -1]];
superPaulis = Array{Float64,2}[makeSuper(i) for i in straightPaulis];
pI=straightPaulis[1]
pX=straightPaulis[2]
pY=straightPaulis[3]
pZ=straightPaulis[4];
piBy8 = exp(im*π/8)*[exp(-im*π/8) 0;0 exp(im*π/8)]
superPiBy8 = makeSuper(piBy8)

# Quick function, given x, checks if its in the cliffords.
# Uses the SuperClifford to ignore phase
function findClifford(x::Array{Complex{Float64},2})
    test = makeSuper(x)
    return findfirst(x->x==test,superCliffords)
end

function findClifford(test::Array{Float64,2})
        return findfirst(x->x==test,superCliffords)
end

function checkFrame(x)
    sum = 0
    for i=1:length(x)
        for j=1:length(x)
            sum += abs(trace(x[i]'*x[j]))^4
        end
    end
    sum/(length(x)^2)
end

function findInThis(needle,haystack)
    tofind = round.(needle,10)
    return findfirst(x->round.(x,10)==tofind,haystack)
end

#Super operator version
function sfindInThis(needle,haystack)
    check = makeSuper(needle)
    return findfirst(x->round.(x,10)==tofind,haystack)

end



gx=expm(im*π*straightPaulis[2]/4)
Gx=makeSuper(gx)
gy=expm(im*π*straightPaulis[3]/4)
Gy=makeSuper(gy)
mgx=expm(-im*π*straightPaulis[2]/4)
mGx=makeSuper(mgx)
mgy=expm(-im*π*straightPaulis[3]/4)
mGy=makeSuper(mgy)


gx=expm(im*π*straightPaulis[2]/4)
Gx=makeSuper(gx)
gy=expm(im*π*straightPaulis[3]/4)
Gy=makeSuper(gy)
mgx=expm(-im*π*straightPaulis[2]/4)
mGx=makeSuper(mgx)
mgy=expm(-im*π*straightPaulis[3]/4)
mGy=makeSuper(mgy)

# My supeoperators are not in the same order as the 'straight' cliffords 
# reorder them.
reorderedSCliffords = copy(straightCliffords)

for i=1:24
    reorderedSCliffords[i]= straightCliffords[findClifford(makeSuper(straightCliffords[i]))]
end


for i=1:24
    assert(makeSuper(reorderedSCliffords[i])==superCliffords[i])
end

straightCliffords=reorderedSCliffords;


#This is our set of Generators (note I included the I)
testSet=[]
push!(testSet,superCliffords[24])
push!(testSet,Gx)
push!(testSet,Gy)


# Not particularly inspiring way to find the minimum generators.

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


# Find the minimum length generators for the Cliffords (1=I 2=Gx 3=Gy)
minGensWithI=[findAGeneratorFor(i) for i=1:24];


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


function findAllGeneratorsFor(clifford)
    gens=[]
    for silly=1:3
    for i=1:3
        for j=1:3
            for k=1:3
                for z=1:3
                        if (testSet[silly]*testSet[i]*testSet[j]*testSet[k]*testSet[z] == superCliffords[clifford])
                            push!(gens,[silly,i,j,k,z])
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

# And the generated 'straight' cliffords should be a unitary-2 design
assert(round.(checkFrame([genA2RotateStraightClifford(i) for i=1:24]),5)==2)


