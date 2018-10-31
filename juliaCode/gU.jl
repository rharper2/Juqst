using HDF5
# Its been saved - to re-run this worksheet, just load it.

include("Initial.jl")
include("Symplectic.jl")
include("open-systems.jl")
include("PTSM.jl")
include("loadHelpers.jl")
FullCliffords= Array{Complex{Float64},2}[]
FullCliffordCommands = Array{String,1}[]
FullCliffordsD=Array{Int64,1}[]
n=11520
for i=0:719 # number of symplectics
    for j=0:15 # number of phase variations.
        state=setup(2)
        decompose(i,j,2,true,false) # symplectic i, bits j, 2 qubits, no output and rationalise
        push!(FullCliffordsD,[i,j])
        push!(FullCliffordCommands,commands)
        push!(FullCliffords,makeFromCommand(commands))
    end
end
SuperCliffs = [makeSuper(s) for s in FullCliffords];
outFile = "test.csv"


first = 0
second = 100

start = first*second+1
toDoTo = (first+1)*second-1
    
# We need to convert the way the noise maps are stored
function convertBack(maps)
    backed=[]
    no_of_maps =size(maps,1)
    for i in 1:no_of_maps
        
        push!(backed,maps[i,:,:])
    end
    return backed
end

FileName = "../storedMaps/highFidelity2QUnitarityv2.hdf"
maps=convertBack(h5read(FileName,"Map"));
higherFidelityFiltered=copy(maps);


FileName = "../TwoQubitSPAMMaps.hdf"
SPAM1=convertBack(h5read(FileName,"Map1"));
SPAM2=convertBack(h5read(FileName,"Map2"));

pauliModifiers = [makeSuper(pI⊗pI),makeSuper(pX⊗pX),makeSuper(pX⊗pI),makeSuper(pI⊗pX)]
rX = pX^0.5
rY = pY^0.5
rZ = pZ^0.5
rZZ = makeSuper(rZ⊗rZ)
rXX = makeSuper(rX⊗rX)
rYY = makeSuper(rY⊗rY)
rXY = makeSuper(rX⊗rY)
rYX = makeSuper(rY⊗rX)
rXZ = makeSuper(rX⊗rZ)
rZX = makeSuper(rZ⊗rX)
rYZ = makeSuper(rY⊗rZ)
rZY = makeSuper(rZ⊗rY);


function extract(v)
    (zz,zo,oz,oo) = v
    return (0.5*(oo+zz-oz-zo),0.5*(oz+zz-zo-oo),0.5*(zo+zz-oz-oo))
    
end

function getA2QCliffordSequence(number)
    return rand(1:11520,1,number)
end

function getA2PauliSequence(number)
    return rand(1:16,1,number)
end

function allocatE(x)
    shots=zeros(4,1)
    for i = 1:1024
        shots[genN(x)] += 1
    end
    return shots./1024
end

function genBinomial(n,p)
    d = Binomial(n,p)
    return rand(d)
end

function genN(r)
    rn = rand()
    for (idx,i) in enumerate(r)
        if rn<i return idx
        else rn=rn-i
        end
    end
end

function getResults(x)
    return real.((initial2ZZ'*x,initial2ZO'*x,initial2OZ'*x,initial2OO'*x))
end
ψz = [1 0]
ψo = [0 1]
ψzz = kron(ψz,ψz)
ψzo = kron(ψz,ψo)
ψoz = kron(ψo,ψz)
ψoo = kron(ψo,ψo)
dzz = ψzz'*ψzz
dzo = ψzo'*ψzo
doz = ψoz'*ψoz
doo = ψoo'*ψoo
initial2ZZ= vec([(trace(_toPauli(_num2quat(i,2))*dzz)) for i = 0:15])
initial2ZZ = 1/norm(initial2ZZ)*initial2ZZ
#print(initial2ZZ)
initial2ZO= vec([(trace(_toPauli(_num2quat(i,2))*dzo)) for i = 0:15])
initial2ZO = 1/norm(initial2ZO)*initial2ZO
initial2OZ = vec([(trace(_toPauli(_num2quat(i,2))*doz)) for i = 0:15])
initial2OZ = 1/norm(initial2OZ)*initial2OZ
initial2OO = vec([(trace(_toPauli(_num2quat(i,2))*doo)) for i = 0:15])
initial2OO = 1/norm(initial2OO)*initial2OO;

shotsToDo = 1024
function allocatU(x)
    shots=zeros(4,1)
    for i = 1:shotsToDo
        shots[genN(x)] += 1
    end
    return shots./shotsToDo
end

positions = [16, 13, 4, 11, 9, 3, 6, 5, 2, 10, 9, 2, 7, 5, 3, 15, 13, 3, 12, 9, 4, 8, 5, 4, 14, 13, 2];

function doAUnitaritySequence(number,noise,noise1,noise2,invert=[SuperCliffs[1],rXX',rYY',rXY',rYX',rZX',rXZ',rYZ',rZY'])
    gates = getA2QCliffordSequence(number)
    results = []
    seqE =  noise*SuperCliffs[gates[1]];
    for g in gates[2:end]
            seqE = noise*SuperCliffs[g]*seqE
    end
    
    # for each sequence we need to take the four computational bases.
    for i in invert
        
        
        p1 = i'*seqE*i*initial2ZZ
        p2 = i'*seqE*i*pauliModifiers[2]*initial2ZZ
        p3 = i'*seqE*i*pauliModifiers[3]*initial2ZZ
        p4 = i'*seqE*i*pauliModifiers[4]*initial2ZZ
        # Simulate results of z measurements on each qubit.
        r1 = allocatU(getResults(p1))
        r2 = allocatU(getResults(p2))
        r3 = allocatU(getResults(p3))
        r4 = allocatU(getResults(p4))
        # Convert this into the results of the pauli state vector
        (p1zz,p1iz,p1zi) = extract(r1)
        (p2zz,p2iz,p2zi) = extract(r2)
        (p3zz,p3iz,p3zi) = extract(r3)
        (p4zz,p4iz,p4zi) = extract(r4)
        # Finally use the four results to get the individual entries in the vector
        
        # This is iz = P3, ZI = P2 and ZZ = p1
        P1 = 0.5*(p1zz+p2zz-(p3zz+p4zz))
        P2 = 0.5*(p1zi+p4zi-(p3zi+p2zi))
        P3 = 0.5*(p1iz+p3iz-(p2iz+p4iz))
        push!(results,[P1,P2,P3])
    end
    return (results,seqE)
    
end

newr = []
for co in start:(toDoTo-1)
    y= mean([map(x->x.^2,doAUnitaritySequence(4,higherFidelityFiltered[co],SPAM1[co],SPAM2[co])[1]) for i=1:5000])
    x2= mean([map(x->x.^2,doAUnitaritySequence(200,higherFidelityFiltered[co],SPAM1[co],SPAM2[co])[1]) for i=1:5000])
    push!(newr,mean(([b for b in Base.Iterators.flatten(x2)]./[a for a in Base.Iterators.flatten(y)]).^(1/195)))
end
outFile = "unis_$(start)_$(toDoTo).csv"

writecsv(outFile,newr)
