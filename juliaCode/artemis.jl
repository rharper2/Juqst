
using HDF5
# Its been saved - to re-run this worksheet, just load it.
    
# We need to convert the way the noise maps are stored
function convertBack(maps)
    backed=[]
    no_of_maps =size(maps,1)
    for i in 1:no_of_maps
        push!(backed,maps[i,:,:])
    end
    return backed
end


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



FileName = "toTestReals.hdf"
testCliffordsR=convertBack(h5read(FileName,"RealCliffords_R"));
testCliffordsI=convertBack(h5read(FileName,"RealCliffords_I"));
testCliffords = testCliffordsR + im*testCliffordsI
noToDo = length(testCliffords)
println("Doing $(first*second+1) to $((first+1)*second)\n")
sum = 0
for i = first*second+1:(first+1)*second
	for j = 1:noToDo
		sum = sum + (abs(trace(testCliffords[i]*testCliffords[j]')))^4
	end
end

print("Sum: $(sum)\n")
outFile = "result_$(first)_$(second).csv"

writecsv(outFile,sum)



