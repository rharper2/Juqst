require("LinearEstimate.jl")

#copied this from Steve's source

#It seems that we take a nxr matrix
#and expand it to one that is nxn but will 
#decompose to nxr (degenerate?)
function randomRankReduced(n,r)
   dim=2^n # n qubits n^2 vectors
   A = randn(dim,r)+im*randn(dim,r)
   (Q,R)=qr(A) # slim is the default
   S=R*R'; S=S/trace(S)
   # diagonalise S
   # note this is different from matlab, Z here is a vector of the eignvalues
   (Z,U)=eig(S)
   # we need to "flesh" it out
   Z = eye(r)
   U=Q*U # redefine U
   return U*Z*U'
 end


function ketPrint(x,num)
  sofar=1;
  x=x-1
  outString=""
  while (x > 0) 
    if (x & 1 == 1)
      outString = string("1",outString);
    else 
      outString = string("0",outString);
    end
    sofar = sofar+1
    x >>= 1
   end
   for i=sofar:num
   	 outString=string("0",outString)
   end
   #outString = reverse(outString)
   print("|",outString,">");
  end

  function showState(a)
  	# shows the state it is in e.g. if the sixth row is 1, it would return the correct ket
  	found = false
  	toReturn = 1
  	for i=1:size(a,1)
  		if a[i] != 0
  			if (found) 
  				println("More than one element of array set, problem")
  				return
  			else 
  				found = true
  				toReturn = i
  			end
  		end
  	end
  	if (!found) 
  		println("Didn't find a set value which is a problem")
  		return
  	else
  		if a[toReturn] < 0 print("-")
  		end
  		ketPrint(toReturn,log2(size(a,1)))
  	end
  end

  function deKet(a)
  	qzero = [1,0]
  	qone = [0,1]
  	daQubit=Array(Int64,2,length(a))
  	for i=1:length(a)
  		if (a[i] == '0')
  			daQubit[:,i]= qzero
  		else 
  			daQubit[:,i]= qone
  		end
  	end
    sofar=daQubit[:,length(a)]
    for loop=length(a)-1:-1:1
      sofar=kron(daQubit[:,loop],sofar)
    end
    return sofar
   end


  function gen7(q1,q2,q3,q4,q5,q6,q7)
  	kron(q1,kron(q2,kron(q3,kron(q4,kron(q5,kron(q6,q7))))))
  end

  function setupSteane()
  	global Px,Py,Pz,Pi,zeroL,oneL
  	Pi = getFullPauli(1)
  	Px = getFullPauli(2)
  	Py = getFullPauli(3)
  	Pz = getFullPauli(4)
  	global g1,g2,g3,g4,g5,g6
    global x1,x2,x3,x4,x5,x6,x7
    global z1,z2,z3,z4,z5,z6,z7
    global xL,zL,yL,iL
	x1 = gen7(Px,Pi,Pi,Pi,Pi,Pi,Pi)
	x2 = gen7(Pi,Px,Pi,Pi,Pi,Pi,Pi)
	x3 = gen7(Pi,Pi,Px,Pi,Pi,Pi,Pi)
	x4 = gen7(Pi,Pi,Pi,Px,Pi,Pi,Pi)
	x5 = gen7(Pi,Pi,Pi,Pi,Px,Pi,Pi)
	x6 = gen7(Pi,Pi,Pi,Pi,Pi,Px,Pi)
	x7 = gen7(Pi,Pi,Pi,Pi,Pi,Pi,Px)

	z1 = gen7(Pz,Pi,Pi,Pi,Pi,Pi,Pi)
	z2 = gen7(Pi,Pz,Pi,Pi,Pi,Pi,Pi)
	z3 = gen7(Pi,Pi,Pz,Pi,Pi,Pi,Pi)
	z4 = gen7(Pi,Pi,Pi,Pz,Pi,Pi,Pi)
	z5 = gen7(Pi,Pi,Pi,Pi,Pz,Pi,Pi)
	z6 = gen7(Pi,Pi,Pi,Pi,Pi,Pz,Pi)
	z7 = gen7(Pi,Pi,Pi,Pi,Pi,Pi,Pz)
	 
  xL =gen7(Px,Px,Px,Px,Px,Px,Px)
  zL = gen7(Pz,Pz,Pz,Pz,Pz,Pz,Pz)
  yL = -gen7(Py,Py,Py,Py,Py,Py,Py)
  iL = gen7(Pi,Pi,Pi,Pi,Pi,Pi,Pi)
  	g1 = kron(Pi,kron(Pi,kron(Pi,kron(Px,kron(Px,kron(Px,Px))))))
	g2 = kron(Pi,kron(Px,kron(Px,kron(Pi,kron(Pi,kron(Px,Px))))))
  	g3 = kron(Px,kron(Pi,kron(Px,kron(Pi,kron(Px,kron(Pi,Px))))))
  	g4 = kron(Pi,kron(Pi,kron(Pi,kron(Pz,kron(Pz,kron(Pz,Pz))))))
  	g5 = kron(Pi,kron(Pz,kron(Pz,kron(Pi,kron(Pi,kron(Pz,Pz))))))
  	g6 = kron(Pz,kron(Pi,kron(Pz,kron(Pi,kron(Pz,kron(Pi,Pz))))))
  	 zeroL = 1/sqrt(8).*(deKet("0000000") + deKet("1010101") + deKet("0110011") + deKet("1100110") + deKet("0001111") + deKet("1011010") + deKet("0111100") + deKet("1101001"))
  	 oneL = 1/sqrt(8).*(deKet("1111111") + deKet("0101010") + deKet("1001100") + deKet("0011001") + deKet("1110000") + deKet("0100101") + deKet("1000011") + deKet("0010110"))
  	 return "Steane setup"
   end
#setupSteane()

  function areTheyTheSame(a,b) 
    return real(trace(a*b)) != 0
  end

  function showMeGs(psi)
  	@printf("Generator 1 returns %f\n",real((psi'*g1*psi)[1]))
  	@printf("Generator 2 returns %f\n",real((psi'*g2*psi)[1]))
  	@printf("Generator 3 returns %f\n",real((psi'*g3*psi)[1]))
  	@printf("Generator 4 returns %f\n",real((psi'*g4*psi)[1]))
  	@printf("Generator 5 returns %f\n",real((psi'*g5*psi)[1]))
  	@printf("Generator 6 returns %f\n",real((psi'*g6*psi)[1]))
end

function testMyHyp()
  for i=1:10
    test1=randomNormalisedVector(1)
    testDM = test1*test1'
    testSteane = test1[1]*zeroL + test1[2]*oneL
    testStDM= testSteane*testSteane'
    @printf("On the I %f then %f\n",real(trace(Pi*testDM)),real(trace(iL*testStDM)))
    @printf("On the X %f then %f\n",real(trace(Px*testDM)),real(trace(xL*testStDM)))
    @printf("On the Y %f then %f\n",real(trace(Py*testDM)),real(trace(yL*testStDM)))
    @printf("On the Z %f then %f\n",real(trace(Pz*testDM)),real(trace(zL*testStDM)))
  end
end


function generateTheStabilisers()
   global stabilisers=[]
   gens = [g1,g2,g3,g4,g5,g6]
   push!(stabilisers,g1)
   push!(stabilisers,g2)
   push!(stabilisers,g3)
   push!(stabilisers,g4)
   push!(stabilisers,g5)
   push!(stabilisers,g6)
   for a1=1:5
    for a2=a1+1:6
         push!(stabilisers,gens[a1]*gens[a2])
     end
   end
   for a1=1:4
      for a2=a1+1:5
        for a3=a2+1:6
                push!(stabilisers,gens[a1]*gens[a2]*gens[a3])
        end
      end
   end
for a1=1:3
  for a2=a1+1:4
    for a3=a2+1:5
      for a4=a3+1:6
                push!(stabilisers,gens[a1]*gens[a2]*gens[a3]*gens[a4])
        
      end
    end
  end
end

for a1=1:2
  for a2=a1+1:3
    for a3=a2+1:4
      for a4=a3+1:5
        for a5=a4+1:6
          push!(stabilisers,gens[a1]*gens[a2]*gens[a3]*gens[a4]*gens[a5])
        end
      end
    end
  end
end

push!(stabilisers,g1*g2*g3*g4*g5*g6)
end

function commute(a1,a2)
   if sum(abs(a1*a2-a2*a1)) == 0 return true
   else return false
   end
end

function commutes(a1)
   if !commute(a1,g1) return false
   end
   if !commute(a1,g2) return false
   end
   if !commute(a1,g3) return false
   end
   if !commute(a1,g4) return false
   end
   if !commute(a1,g5) return false
   end
   if !commute(a1,g6) return false
   end
   return true
end

function inStabilisers(a1) 
   for test in stabilisers
      if areTheyTheSame(a1,test) return true
      end
   end
   return false
end

function checkPaulis()
  lengthOfPaulis = size(ALLPAULI,2) #number of tensored matrixes.
   global nonCommuters=[]
   global pauliInStabilisers=[]
   global commutingNonStabilisers=[]
  dim=convert(Int64,sqrt(lengthOfPaulis))
  if (dim !=128 )
      print("Just now we are only checking the paulis against steane code generators, so you need to initialise(7) make a hot cup of tea and try again")
      return
   end
    count = 0
  
  for i=1:lengthOfPaulis
    test = reshape(ALLPAULI[:,i],dim,dim)
    if !commutes(test) 
      count = count+1
      #@printf("The following [%d]: %s did not commute\n",count,AllPauliNames[i])
      push!(nonCommuters,i)
      if i % 2000 == 0 
         @printf("Checked %d\n",i)
      end
   
   elseif inStabilisers(test) 
      push!(pauliInStabilisers,i)
   else 
      @printf("%s commuted but was not in stablilisers\n",AllPauliNames[i])
       push!(commutingNonStabilisers,i)
    end
      
  end
  @printf("Total non-commuting is %d\n",count)
 
end

function getThePauli(i)
    dim=convert(Int64,sqrt(lengthOfPaulis))
  if (dim !=128 )
      print("Just now we are only doing this paulis against steane code generators, so you need to initialise(7) make a hot cup of tea and try again")
      return
   end
   return reshape(ALLPAULI[:,i],dim,dim)
end

#Try this with 7 qubits, if you have a day to spare.
#
#function depolorise(p,j,n)
#   q=p/4
#   tX = kron(eye(2^(j-1)),Px)
#   tX = kron(tX,eye(2^(n-j)))
#   tY = kron(eye(2^(j-1)),Py)
#   tY = kron(tY,eye(2^(n-j)))
#   tZ = kron(eye(2^(j-1)),Pz)
#   tZ = kron(tZ,eye(2^(n-j)))
#   return (1-3*q)*eye(4^n) + q*(kron(tX,tX) + kron(tY,conj(tY)) + kron(tZ,tZ))
#end

function depol(p,n,rho)
   srho = sparsevec(vec(rho))
   for j=1:n
      srho = spDepolorise(p,j,n)*srho
   end
   return dense(reshape(srho,2^n,2^n))
end


 function s7labStats(s1,n,noise)
  # Using the matlab interface at s1, calls the script
  # specified in sendItToMatlab to use Matlab
  # to work out the reconstructed density matrix.
  # Averages n times, using a qubits system.
  # I am going to assume that initialise(qubits) has been called
  global densityMatrix
  qubits=7
  multiple=1000
  indexC = 11
  a=zeros(Float64,indexC)
  println(a)
  
  for j=1:n
    count =1
    i=5
    while (count <= indexC)
     testVector=randomNormalisedVector(1)
     testSteane = testVector[1]*zeroL + testVector[2]*oneL
     testStDM= testSteane*testSteane'
     # here we might want to restrict the paulis we use or add depolorising noise
      @printf("Construction density Matrix\n")
     densityMatrix=depol(noise,7,testStDM)
     fnameDR = string("DensityReal3_",count)
     fnameDI = string("DensityImag3_",count)
     
     writedlm(fnameDR,real(densityMatrix),',')
     writedlm(fnameDI,imag(densityMatrix),',')
     @printf("Now with added noise!\n")
     
         (R,U) = constructRj(i)
         @printf("Sending to matlab count %f\n",count);
         phHat = sendItToMatlab(s1,R,U)
         fnameR = string("ReturnedReal2_",count);
         fnameI = string("ReturnedImag2_",count);
         interimFname = string("Interim",count,".ssv")
         writedlm(fnameR,real(phHat),',')
         writedlm(fnameI,imag(phHat),',')
         @printf("Got back\n")
         diff = phHat-densityMatrix
         ret = sqrt(real(trace(diff'*diff)))
         a[count]=a[count]+ret
         count = count+1
         i = i+multiple
         writedlm(interimFname,a,' ')
      end
      println(a)
  end
  println(":",j)
  av = a/n
  writedlm("Average.ssv",av)
 return a/n
end 

function generateSparseBasePaulis()
  global sPi,sPx,sPz,sPy
  sPi = sparse([1,2],[1,2],[1.0+0*im,1.0+0*im])
  sPx = sparse([1,2],[2,1],[1.0+0*im,1.0+0*im]) 
  sPz = sparse([1,2],[1,2],[1.0+0*im,-1.0+0*im]) 
  sPy = sparse([1,2],[2,1],[-1.0*im,1.0*im]) 
  global sparsePaulis=[sPi,sPx,sPy,sPz]
end

function spDepolorise(p,j,n)
   q=p/4
   tX = kron(speye(Complex{Float64},2^(j-1)),sPx)
   tX = kron(tX,speye(Complex{Float64},2^(n-j)))
   tY = kron(speye(Complex{Float64},2^(j-1)),sPy)
   tY = kron(tY,speye(Complex{Float64},2^(n-j)))
   tZ = kron(speye(Complex{Float64},2^(j-1)),sPz)
   tZ = kron(tZ,speye(Complex{Float64}, 2^(n-j)))
   return (1-3*q)*speye(4^n) + q*(kron(tX,tX) + kron(tY,conj(tY)) + kron(tZ,tZ))
end





   