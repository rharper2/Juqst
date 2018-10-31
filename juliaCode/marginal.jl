module Marginal

using LsqFit, Hadamard, DelimitedFiles, LinearAlgebra
using PyPlot
import Base.show

export readInCSVFile,transformToFidelity,fitTheFidelities,convertAndProject,gibbsRandomField
export getGrainedP,mutualInformation,hinton,relativeEntropy,conditionalMutualInfo,covarianceMatrix

⊗ = kron

""" 
readInCSVFile(filename)

## Arguments 
	Filename: full name of a csv delimited Filename

## File Format
	Assumes that each sequence appears on a new line and each line has the same number of entries.
	Assumes that it is a raw count for each of the possible measurement outcomes.
	Note it is important to know the order of the outcomes (this is significant).
	Although this can be overriden, it is a assumed that they are stored as:

```
			00000000
			00000001
			00000010
			00000011
			
			...
			11111110
			11111111

		i.e. lsb to the right.
```
Returns a full rows x col matrix depending on the number of entries of the file.
for a 16 qubit machine, where the data was taken over 6 sequences we will return a 6x655536 array of Ints.

"""
function readInCSVFile(filename)
	return readdlm(filename,',',Int64)
end

"""
	transformToFidelity(data)

##	Arguments:
	Data, as returned by readInCSVFile - see help for that function.
		
## Returns:
	The Hadamard transformed fidelities for each sequence as a list of lists.

"""
function transformToFidelity(data)
	splitMatrix = [data[i,:]/sum(data[i,:]) for i in 1:size(data,1)]
	return [ifwht_natural(x) for x in splitMatrix]
end


modelF(x, p) = p[1]*(p[2].^x)

"""
	fitTheFidelties(transformedData,lengths)

## Arguments
-	`transformedData: Array{Array{Float64,1},1}`, think of as a List of List of floating point observations (see transformToFidelity)

-	`lengths: Array{Int64,1}`, e.g. collect(1:2:24) - the length at which each of the above sets of data was gathered.

Take the transformed data and see if we can fit it to an exponential decay.
This proves problematic where the data is a bit sparse, as the decays are all over the place.
	
Here we check for convergence, if that fails we arbitrarily set the A of Af^m to 1.
If we have to do the above we @warn.

## Returns
	three things:
	-  the fitting parameters, 
	-  the indices of those that we fit with A=1 
	-  the indices of those that still failed to converge.

## Typical usage
```
	params,_ = fithTheFidelities(transformToFidelity(data),collect(1:2:22))
```
"""
function fitTheFidelities(transformedData,lengths)
	# Fit the Fs. Pretty naive fitting here.
	# This default weight vector that I found to help with a specific dataset.
	
	params= Array{Float64,1}[]
	push!(params,[1.0,1.0])
	secondfit=[]
	failedToC = []
	for i = 2:length(transformedData[1])
		extract = map(x->x[i],transformedData)
    	if (transformedData[1][i] < 0 ) 
    		@warn "Index $i started with a negative value, which might be problematic\n"
		end
    	fit = curve_fit(modelF, lengths,extract, [1,0.6],lower=[.4,0.03],upper=[1.0,1.0])
	    if !fit.converged
            # Just force A to be 1.
                fit = curve_fit(modelF,lengths,extract,[1,0.3],lower=[1.0,0.15],upper=[1.0,1.0])
                if fit.converged
                    push!(secondfit,i)
                else
                    push!(failedToC,i)
                end
       
	    end
    	push!(params,fit.param)
    end
	if length(secondfit) + length(failedToC) > 0
		@warn "Forced: $(length(secondfit)) entries to have A=1, and of these $(length(failedToC)) failed to converge.\n"
	end

	return params, secondfit, failedToC
end



"""
	Projects the P's onto the simplex where they are all >= 0
"""
function ProjectSimplex(ps)
    sortedSeq = reverse(sort(ps))
    # This is a bit crytic.
    #   a) Sort (greatest first)
    #   b) create a culmlative sum (cumsum)
    #   c) we want the last one where (Cum[k]-1)/k < sorted[k]
    #   d) we then reduce all probs by the last value we had in c, setting them to zero if they are negative.
    τ=filter(x->x[1]<x[2],collect(zip((cumsum(sortedSeq).-1)./(1:length(ps)),sortedSeq)))[end][1]
    projected = map(x->max(x-τ,0),ps)
    return projected
end




"""
	Convert and project the fitted fidelities, to rerieve the 'p's representing the joint probabilities.

	Input: The parameters that have been fit.
	Output: The projected ps.
"""
function convertAndProject(params)
	ps = fwht_natural(map(x->x[2],params))
	pps = ProjectSimplex(ps)
	if !(isapprox(sum(pps),1))
		@warn "There may be a problem, the probabilities do not add to 1."
	end
	return pps
end


function getIndices(qs;dimension=16)
    sqs = sort(qs)
    indices = [2^(sqs[1]-1),2]
    for i in 2:length(sqs)
            push!(indices,2^(sqs[i]-sqs[i-1]-1))
            push!(indices,2)
    end
    push!(indices,2^(dimension-sqs[end]))
    return indices
end


""" 
    Returns the marginalised probability, given the (probjected) probabilities in pps.
    Now with optional second parameter - default is pps i.e. the project (global) ps
"""
function marginalise(q,pps)
    # Get indices sorts the entries
    dimension = 0
    try
    	dimension = Integer(log(2,length(pps)))
    catch
    	@warn "Error, the size of pps needs to be an integer power of 2"
    	return 
    end
    x = reshape(pps,tuple((getIndices(q,dimension=dimension))...))
    # Below the permutation that will get us back.
    sorted = sort(q)
    # This is the reverse of sortperm is from sorted -> original
    permute = [findfirst(isequal(x),sorted) for x in q]
  
    # print(permute)
    # Sum over all the odd ones.
    marginalised = sum(x,dims=1:2:length(size(x)))
    # Work out how many 2 variables we have
    indices = length(q)
    # Shove them in a tuple (note splat so it works well with reshape)
    fullIndices = tuple([2 for i in 1:indices]...)
    # Might as well get it in a 2 dim array to return
    rows = round.(Int,floor(sqrt(2^indices)))
    cols = round.(Int,2^indices/rows)
    return reshape(permutedims(reshape(marginalised,fullIndices),permute),rows,cols)
end


"""
gibbsRandomField(pps,constraints)

## Arguments
-	`pps: Array{Float64,1}` For example the output of convertAndProject
-	`constraints: Array{Int64,1}`: The division of the gibbs field. See below for example.

	Takes a joint probability and the gibbs variable constraints and returns a vector of reduced probability distributions.
	For example the contrains might be [[1,2,3,4],[3,4,5,6],[5,6,7,8]] over a field of 8 qubits
	This would return an array of 3 ϕ's each 16 long (2^4), from which the joint probability can be extracted on the assumption
	That, say, qubits 1,2 are independent of qubits 5-8.
"""
function gibbsRandomField(pps,constraints)
	pxx = [marginalise(x,pps) for x in constraints]
	px =  [vec(marginalise(x,pps))' for x in [i[3:end] for i in constraints[1:end-1]]]
	# Map any zero entries in px to something non-zero or division will fail.
	# Not an issue as if px is zero, the corresponding pxx has to be zero. (and in this case by definition the probability is zero)
	pxNz = map(x-> x == 0 ? 1e-8 : x,px)
	ϕ = [vec(pxx[i]./px[i]) for i in 1:(length(constraints)-1)]
	push!(ϕ, vec(marginalise(constraints[end],pps)))
	return ϕ
end


# The function below is probably worth a bit of comment. 
# uses the information in 'graining' to tells us which index the pattern we want resides in
# For instance say we want all bits zero apart from qubit 3 (which = 1)
# Then for ϕ[1], which represents bits [2,1,3,6] - we need the entry corresponding to decimal equiv of reverse(0,0,1,0)
# plus one, because 0 0 0 0 = first entry i.e. index 1.
# We will have length(ϕ) of these entries, this gives us the perVecIndex.
# Then we just read it out for ϕ, with a map, using foldl to multiply them for us.


"""
	For a given ϕ, the list of qubits we used to construct ϕ (graining) and the 
    Bit pattern, the probability of which we want to extract
"""
function getGrainedP(ϕ,tomatch,graining)
    vectorIndex = [[x&tomatch>0 ? 1 : 0 for x in (map(x->2^(x-1),i))] for i in graining]
    perVecIndex = map(x->parse(Int64,join(reverse(x)),base=2)+1,vectorIndex)
    return foldl(*,map(x->ϕ[x][perVecIndex[x]],1:length(ϕ)))
end

"""
Some gotchas:
    I return -1 if p1==p2 , this allows an identification of the qubit where I loop.
    if p12 is 0, then we get Nan error, replace with zero on the basis that although the logs
    are undefined, we are multiplying them with zero.
"""
function mutualInformation(p1,p2,p)
    # Allow this for easy looping and identificaiton of controlling qubit.
    if p1 == p2 
        return -1
    end
    p12 = marginalise([p1,p2],p)
    p1 = sum(p12,dims=2)
    p2 = sum(p12,dims=1)
    logPart = (log.(vec(p12)) - log.(vec(p1⊗p2)))
    # Zap the ones where p12 is zero (gets rid of NAN)
    
   toReturn =  vec(p12)'*(log.(vec(p12)) - log.(vec(p1⊗p2)))
   return map(x->isnan(x) ? 0 : x,toReturn)
end


"""
Adapted from demo algorithm found at 
https://matplotlib.org/gallery/specialty_plots/hinton_demo.html

No default on ax or max_weight - pass it in (allows different plots to keep it constant)

Modified to accept keyword, which allows for a fixed highlight of negative values.

## Arguments:
-	`matrix: Array{Float64,1}` of values to be plotted
-	`max_weight: Float` The value required to 'fill a box'
-	`ax`, the graphical axis (get from gca())
-	`highlightNegative`: If true will print negative values as a qubit, value given by qubit and font by fontsize.
-	`addScale`: If true will draw the scale of the box (from max_weight) to the rhs of the plot.
    
"""
function hinton(matrix, max_weight, ax;highlightNegative = false,fontsize = 12,qubit=0,addScale=false)
    #"""Draw Hinton diagram for visualizing a weight matrix."""
    matrix=transpose(matrix)
    ax[:patch][:set_facecolor]("gray")
    ax[:set_aspect]("equal", "box")
    ax[:xaxis][:set_major_locator](plt[:NullLocator]())
    ax[:yaxis][:set_major_locator](plt[:NullLocator]())
    for x in 1:size(matrix)[1]
        for y in 1:size(matrix)[2]
            color="lightgray"
            size = 1.0
            rect = plt[:Rectangle]([x - size / 2, y - size / 2], size, size,
                             facecolor="gray", edgecolor="lightgray")
            ax[:add_patch](rect)
        end
    end
    if addScale
        ax[:add_patch](rect)
        rect = plt[:Rectangle]([size(matrix)[1]+xoffset+0.5,yoffset+1.3],0.1,0.1,
                facecolor = "gray")
    
        ax[:add_patch](rect)
    
        plt[:annotate](s="", xy=(size(matrix)[1]+0.7,1.5), xytext=(size(matrix)[1]+0.7,0.5), arrowprops=Dict([(:arrowstyle,"<->")]))
        plt[:annotate]( string(max_weight), 
                    xy=(size(matrix)[1]+0.8,0.8),
                    xytext=(size(matrix)[1]+0.8,0.8), 
                    va = "top", ha="left")

    end
    for x in 1:size(matrix)[1]
        for y in 1:size(matrix)[2]
            w= matrix[x,y]
            size = min(sqrt(abs(w) / max_weight),1)
            if w >= 0 
                color = "white" 
                rect = plt[:Rectangle]([x - size / 2, y - size / 2], size, size,
                             facecolor=color, edgecolor= color)
                ax[:add_patch](rect)
            else 
                if highlightNegative
                    color = "blue"
                    size=0.7
                    circ = plt[:Circle]([x - 0.01  , y - 0.03  ], size/2,
                             facecolor=color, edgecolor= "black")
                    ax[:add_patch](circ)
                    plt[:text](x,y,"Q"*string(qubit),color="white",fontsize=fontsize,
                                horizontalalignment="center",
                                verticalalignment="center"
                                )
                else
                    color = "black"
                    rect = plt[:Rectangle]([x - size / 2, y - size / 2], size, size,
                             facecolor=color, edgecolor= color)
                    ax[:add_patch](rect)
                end
            end
            
        end
    end
    ax[:autoscale_view]()
    ax[:invert_yaxis]()

end


"""
	relativeEntropy(jointProbabilityDistribution1, jpS2)

	Calculates the relative entropy between two joint probability distributions.
	D(P||Q) = Sum p_j*log(p_j/q_j)

	Undefined if any of the q_js are zero, unless the corresponding p_j is zero (incomplete statement of the rule, but that's the relevant part.)
"""
function relativeEntropy(pps,p̃)
	# First of all assert that where p̃ is zero pps is also zero
	@assert(length(filter(x->x!=0,pps.*p̃)) == length(filter(x->x!=0,pps)))

	# We only need the values of p and p̃ where pps is non zero
	xp = pps[map(x->x!=0,pps)]
	xp̃ = p̃[map(x->x!=0,pps)]
	pDp̃ = xp./xp̃
	return xp'*log.(pDp̃)
end





""" 
    Helper function that given an xyz correctly indexes and pulls out the
    xyz from the joint probability distribution, jpXYZ
    and the xz and yz from jpXZ and jpYZ 
    then does the calculation, inside the sum for that x,y,z
"""
function getSummand(x,y,z,jpXYZ,jpXZ,jpYZ,jpZ)
    multipland = size(jpYZ)[1]
    pxyz = jpXYZ[multipland*(y-1)+x,z]
    if pxyz == 0 
        return 0
    end
    pz = jpZ[z]
    pxz = jpXZ[x,z]
    pyz = jpYZ[y,z]
    return pxyz*log((pz*pxyz)/(pxz*pyz))
end




""" 
    Gives the conditional mutual information for the qubits 
    supplied in X,Y and Z where we want I(XY|Z)

	I(XY|Z) = Sum( P(x,y,z)log( P(Z)P(X,Y,Z)/P(X,Z)P(Y,Z) )


    All the heavy lifting is in getSummand.
"""
function conditionalMutualInfo(X,Y,Z,p)
    XYZ = vcat(X,Y,Z)
    XZ = vcat(X,Z)
    YZ = vcat(Y,Z)
    xSize = 2^length(X)
    ySize = 2^length(Y)
    zSize = 2^length(Z)
    jpXYZ = reshape(newMarginalise(XYZ,p),:,zSize)
    jpXZ = reshape(newMarginalise(XZ,p),:,zSize)
    jpYZ = reshape(newMarginalise(YZ,p),:,zSize)
    jpZ = reshape(newMarginalise(Z,p),1,zSize)
    return sum([getSummand(x,y,z,jpXYZ,jpXZ,jpYZ,jpZ) for x=1:xSize,y=1:ySize,z=1:zSize])
end

"""
	translateLocationd(data,qubit)

Rather specific IBMQX16 function. Basically remaps the qubits into the same order as they appear in the IBM 16 Qubit machine.
"""
function translateLocation(data,i)
    a = Array{Float64}(undef,2,8)
    a[2,1]=data[1]
    a[1,1]=data[2]
    a[1,2:8]= data[3:9]
    a[2,2:8]=reverse(data[10:16])
    
    return a
end



""" 
	drawIBMHinton(titleText)

Another specific function (specific to IBMQX16). Basically draws hinton plots for each of the qubits, given the mutualInformation Ps.
If you specify a save file it will save it.

"""
function drawIBMHinton(mutualPs,titleText;saveFile = "",titleFontSize=24)
	addScale = false
	fig = figure("Slightly larger",figsize=(18,6))
	scale = round(maximum(maximum.(mutualPs)),RoundUp,digits=3)
	fig[:suptitle](titleText+" Full box = $scale",fontsize=titleFontSize)
	ax = gca()
	ax[:set_facecolor]("gray")
	for (ix,i) in enumerate(2:4)
  	  subplot(4,4,ix)
  	  ax=gca()
   	 hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1)
	end


	for (ix,i) in enumerate(5:5)
   	 subplot(4,4,4)
  	  ax=gca()
  	  hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1,addScale=addScale)
	end
	for (ix,i) in enumerate(6:9)
    	subplot(4,4,4+ix)
    	ax=gca()
    	hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1)
	end
	for (ix,i) in enumerate(vcat([1],16:-1:14))
	    subplot(4,4,8+ix)
	    ax=gca()
	    hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1,
        	fontsize=i > 10 ? 9 : 12)
	end
	for (ix,i) in enumerate(13:-1:10)
	    subplot(4,4,12+ix)
	    ax=gca()
	    hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1,
        	fontsize=i >10 ? 9 : 12)
	end
	fig[:subplots_adjust](top=1.3)
	plt[:tight_layout]()
	if saveFile != ""
		savefig(saveFile)
	end
end



"""
	covarianceMatrix(p)

## Covariance
	We can compute the covariance matrix between the 0/1 random variables representing no error / error. 
	If x is a column vector of bits representing an error pattern, then we can compute the matrix
	Expect_p[(x-μ) (x-μ)^T]
	where μ = Expect_p[x]. 
"""
function covarianceMatrix(p)
# Take a probability distribution on n qubits and compute the covariance matrix as above.
    n = convert(Int,log2(length(p)));
    v = reverse.(digits.(0:(2^n-1),base=2, pad=n)); # the 0/1 random variable of errors
    v = convert(Array{Array{Int,1},1}, v);
    mu= sum(p.*v); # compute the mean
    v = hcat(v...).-mu; # center the random variable v
    covmat = p[1]*v[:,1]*v[:,1]';
    for i in 2:2^n
        covmat+=p[i]*v[:,i]*v[:,i]';
    end
    return covmat
end


end
