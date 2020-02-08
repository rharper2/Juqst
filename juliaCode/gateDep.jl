 """
# extractLRD(perfectGates,noisyGates)

Code written to implement the ideas presented in: Joel Wallman - Randomized benchmarking with gate-dependent noise Arxiv:1793.09835


   (Left,Right,Delta) = extractLRD(perfectGates,noisyGates)
   
   Given a set of noisy gates (that form a unitary-2 design) extract
   the averaged twirled noise channel between pairs of the gates (see https://arxiv.org/pdf/1703.09835.pdf)

 
## Equations 13(a)--(c)

``\\mathbb{E}_G(\\mathcal{\\tilde{G}}\\mathcal{L}\\mathcal{G}^{\\dagger})=\\mathcal{L}\\mathcal{D}_{p,t}``

``\\mathbb{E}_G(\\mathcal{G}^{\\dagger}\\mathcal{L}\\mathcal{\\tilde{G}})=\\mathcal{D}_{p,t}\\mathcal{R}``

``\\mathbb{E}_G(\\mathcal{G}^{\\dagger}\\mathcal{R}\\mathcal{L}\\mathcal{G})=\\mathcal{D}_{p,t}``

## Become: 17(a)--(d)

``\\mathbb{E}_G(\\mathcal{\\tilde{G}}) | L \\rangle\\!\\rangle = t|L\\rangle\\!\\rangle``

``\\mathbb{E}(\\mathcal{\\tilde{G}}\\mathcal{L}^\\prime{}\\mathcal{G}^{\\dagger}_{u})=p\\mathcal{L}^{\\prime{}}``

``\\langle\\!\\langle R | \\mathbb{E}_G(\\mathcal{\\tilde{G}})=t\\langle\\!\\langle R|``

``\\mathbb{E}_G(\\mathcal{G}_U^{\\dagger}\\mathcal{R}^{\\prime{}}\\mathcal{\\tilde{G}})=p\\mathcal{R}^{\\prime}``


## Which we use as (20)

``\\mathbb{E}_G(\\mathcal{G}_u \\otimes \\mathcal{\\tilde{G}})\\text{ vec}(\\mathcal{L}^{\\prime})=p\\text{ vec}(\\mathcal{L}^{\\prime})``

``\\mathbb{E}_G(\\mathcal{\\tilde{G}} \\otimes \\mathcal{G}_u)^{T}\\text{ vec}(\\mathcal{R}^{\\prime})=p\\text{ vec}(\\mathcal{R}^{\\prime})``

  
	


## Arguments

-----

	perfectGates: a vector of each of the "perfect" gates (in the PauliLiouville superoperator basis)

	noisyGates: a vector of each of the gates (in the PauliLiouville superoperator basis)

------	


"""
function extractLRD(expectedGates,noisyGates) # the documentation is perfectGates, the code expectedGates, sorry
	dsq = size(noisyGates[1])[1] # dimension squared
	B1 = zeros(Int64,dsq,1)
	B1[1]=1
	B11B=B1*B1'
	numberOfGates = length(expectedGates)
	Gu=[expectedGates[i]-B11B for i=1:numberOfGates];
	Expect1 = mean([kron(Gu[i],noisyGates[i]) for i = 1:numberOfGates])
	Expect2 = mean([transpose(kron(noisyGates[i],Gu[i])) for i=1:numberOfGates])
	ExpectG = mean(noisyGates)
	# We need the left hand eigenvector
	ExpectGT = collect(ExpectG')

	eig1=eigen(Expect1)
	eig2=eigen(Expect2)
	eig3=eigen(ExpectG)
	eig4=eigen(ExpectGT)

	# So repectively these are the L', R', L and R values (17b,17d,17a,17c of Joel's paper.)
	Lp_position = findmax(abs.(eig1.values))[2]
	Rp_position = findmax(abs.(eig2.values))[2]
	M1_position = findmax(abs.(eig3.values))[2]
	M2_position = findmax(abs.(eig4.values))[2]

	Lp=eig1.vectors[:,Lp_position]
	Rp=eig2.vectors[:,Rp_position]
	M1=eig3.vectors[:,M1_position]
	M2=transpose(eig4.vectors[:,M2_position])

	# The p's - before we rescale to make them equal.
	LV=eig1.values[Lp_position]
	RV=eig2.values[Rp_position]

	MV1=eig3.values[M1_position]
	MV2=eig4.values[M2_position]


	Lprime = reshape(Lp,dsq,dsq)
	Rprime = reshape(Rp,dsq,dsq)
    @assert(isapprox(LV,RV))

	Rscaled=(dsq-1)*RV/tr(Lprime*Rprime)*Rprime
	η1=MV1/(M1[1])
	η2=MV2/(M2[1])

	finalR = η2*B1*M2+Rscaled
	finalL = η1*M1*B1'+Lprime
	deltaG=[expectedGates[i]'*finalR*(noisyGates[i]-(finalL*expectedGates[i]*finalR)) for i=1:numberOfGates];
	if !(isapprox(round(sum(abs.(mean(deltaG))),digits=6),0))
		@warn "Warning, check deltaG values\n"
	end
	return (finalL,finalR,deltaG)
end


"""
# extractLRD(noisyGates)

##	Form:

### (finalL,finalR,deltaG) = extractLRD(noisyGates)

Everything in the PauliLiouville superoperator basis.

Kept for compatibility reasons - prefer *extractLRD(expectedGates,noisyGates)* 

"""

function extractLRD(noisyGates)
	
	return extractLRD(superCliffords,noisyGates)
end

