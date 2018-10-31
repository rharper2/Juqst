
include("Initial.jl")
include("Symplectic.jl")


# Its been saved - to re-run this worksheet, just load it.
    
# We need to convert the way the noise maps are stored




	if length(ARGS) != 2
		println("Expected 2 arg got $(length(ARGS))\n")
		exit()
	end
	
	first = parse(Int64,ARGS[1])
	if (first < 0)
		println("Unexpected to start with $first\n")
		exit()
	end

	second = parse(Int64,ARGS[2])
	if (second < 1 || second > 100000)
		println("Unexpected batch size $second\n")
		exit()
	end



FileName = "4QReals.csv"
reals_totest = readcsv("4QReals.csv",Int64)
combos = readcsv("combos.csv",Int64)
pI = [1 0;0 1]
pX = [0 1;1 0]
pY = [0 -im;im 0]
pZ = [1 0;0 -1]
paulis = [pI,pX,pY,pZ]

function gen256(no,combos,paulis)
	svec = stabiliseSymp(symplectic(no,4),0)
	decomposeState(svec,true)
	x = makeFromCommand(commands)
	cliffords = []
	for j in 1:size(combos)[1]
        i = combos[j,:]
		todo = copy(x)
		todo = todo * kron(paulis[i[1]],kron(paulis[i[2]],kron(paulis[i[3]],paulis[i[4]])))
		push!(cliffords,todo)
	end
	return cliffords
end


function rawCliffords(no)
	svec = stabiliseSymp(symplectic(no,4),0)
	decomposeState(svec,true)
	x = makeFromCommand(commands)
	return x
end

function ptwirl(x)
	cliffords = []
	for j in 1:size(combos)[1]
        i = combos[j,:]
		todo = x * kron(paulis[i[1]],kron(paulis[i[2]],kron(paulis[i[3]],paulis[i[4]])))
		push!(cliffords,todo)
	end
	return cliffords
end

first = first+288
println("Doing $(first*second+1) to $((first+1)*second)\n")

# Actually just create the first set of cliffords
# Note here we are already conjugating them.
toTestCliffords = []
for i = (first*second+1):((first+1)*second)
	got = gen256(reals_totest[i],combos,paulis)
	for c in got
		push!(toTestCliffords,c')
	end
end

baseCliffords = [rawCliffords(x) for x in reals_totest]

print("Created twirlers")

sum = 0
for j in baseCliffords
	for i in toTestCliffords
		interim = ptwirl(j)
		for c in interim
			sum = sum + (abs(trace(i*c)))^4
		end
		
	end
end

print("Sum: $(sum)\n")
outFile = "result_$(first)_$(second).csv"

writecsv(outFile,sum)



