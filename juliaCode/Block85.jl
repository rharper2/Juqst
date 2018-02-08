
module B85


export complexFNT85

DigitizingFactor = 100 # factor to take into account we have 24 fp bits in the digitizer.
WeightsFactor = 128 # controls the weights.


export do85Block

include("RF5.jl")
include("B17.jl")
include("B5.jl")
using B17
using B5

N1=17
N2=5
K1=N2*(inv(ModInt{N1}(N2))).k
K2=N1*(inv(ModInt{N2}(N1))).k

index1=[(n1*N2+n2*N1)%(17*5) for n2=0:4,n1=0:16]

index2=[(k1*K1+k2*K2)%(17*5) for k2=0:4,k1=0:16]

# This function just pulls out the sample pointed to by the index
function extract(samples,index,no)
    return [samples[j+1] for j in index[no,:]]
end

function do85Block(samplesR,samplesI)

	inputR1 = digitiseAndReduceSamples(extract(samplesR,index1,1),DigitizingFactor)
	inputR2 = digitiseAndReduceSamples(extract(samplesR,index1,2),DigitizingFactor)
	inputR3 = digitiseAndReduceSamples(extract(samplesR,index1,3),DigitizingFactor)
	inputR4 = digitiseAndReduceSamples(extract(samplesR,index1,4),DigitizingFactor)
	inputR5 = digitiseAndReduceSamples(extract(samplesR,index1,5),DigitizingFactor)
	inputI1 = digitiseAndReduceSamples(extract(samplesI,index1,1),DigitizingFactor)
	inputI2 = digitiseAndReduceSamples(extract(samplesI,index1,2),DigitizingFactor)
	inputI3 = digitiseAndReduceSamples(extract(samplesI,index1,3),DigitizingFactor)
	inputI4 = digitiseAndReduceSamples(extract(samplesI,index1,4),DigitizingFactor)
	inputI5 = digitiseAndReduceSamples(extract(samplesI,index1,5),DigitizingFactor)

	FNTFT1=complexFNT17(vcat(inputR1,inputI1))
	FNTFT2=complexFNT17(vcat(inputR2,inputI2))
	FNTFT3=complexFNT17(vcat(inputR3,inputI3))
	FNTFT4=complexFNT17(vcat(inputR4,inputI4))
	FNTFT5=complexFNT17(vcat(inputR5,inputI5));

	
	# reorganise so thats 17 by (5 real x 5 imaginary) (i.e. 5 part input)
	firstOutputs=[]
	for i=1:17
 		push!(firstOutputs,vcat(FNTFT1[i],FNTFT2[i],FNTFT3[i],FNTFT4[i],FNTFT5[i],FNTFT1[i+17],FNTFT2[i+17],FNTFT3[i+17],FNTFT4[i+17],FNTFT5[i+17]))
	end


	# Then do the 5 part on each of these 17 groupings.
	results = []
	for i=1:17
    	results = vcat(results,[complexFNT5(firstOutputs[i])])
	end

	# reorder to correct output (keeping real and imaginary seperate just now.)
	outputReal = zeros(Int64,85,1)
	outputImaginary = zeros(Int64,85,1)
	
	for no=1:17
        for j=1:5
            #print("Setting X[$(index[j,no]+1)] to equal samples[$(4*(no-1)+j)]\n")
            outputReal[index[j,no]+1] = reducedToInt(samples[no][j])
            outputImaginary = reducedToInt(samples[no][j+5])
        end
    end
    

end

function finito(samples,index)
    # Use the index (Good-Thomas unscramble index) to re-order the output signals.
    results = zeros(Complex{Float64},85,1)
    for no=1:17
        for j=1:5
            #print("Setting X[$(index[j,no]+1)] to equal samples[$(4*(no-1)+j)]\n")
            results[index[j,no]+1] = reducedToInt(samples[no][j])+im*reducedToInt(samples[no][j+5])
        end
    end
    return results
end


end




