#Basic use here is intialise(number_of_qubits)
#Alternatively we can use funStats(number_of_iterations_to_average_over,qubits)
#this provides an array of 1..2^numberofqubits showing
#the error between the starting initial random vector
#and the rebuilt one, using the index number of the constructed
#Pauli tensor.

include("densityMatrixCreation.jl")
# require("sparseDensityMatrixCreation.jl")
include("densePauliFunctions.jl")

# Initialise the system (setting up globals etc) for n quibits
function initialise(n::Int64)
    global vectorSize = 2^n
    global pauliSize = n
    global qubits = n
    global AllPauliNames
    global FullBasePaulis=generateFullBasePaulis()
    
    initialiseVandDM(n)
    
  # when we were generating the Paulis, I omitted the 1/sqrt(2) factor
  # so introduce it here.
    global ALLPAULI
    global AllPauliNames 
    (ALLPAULI,AllPauliNames) = getAllPaulis(pauliSize) 
    ALLPAULI = ALLPAULI * *(1/sqrt(2))^n
  # once I know I have it right, I should save these somewhere.
 
 end

function initialiseVandDM(n)
  global testVector=randomNormalisedVector(n)
    global densityMatrix = testVector*testVector'
  end

function constructRj(m) 
# m is the number of Pauli matrices we want to use
# obviously this can't be larger than size(ALLPAULI,2) - i.e the number 
# number of created Paulis.

  length = size(ALLPAULI,2) #number of tensored matrixes.
  toUse = randperm(length)[1:m]
  Rjs = constructRjWith(m,toUse)
  #return the results and the Paulis used as a tuple.
  return (Rjs,toUse)
end

function constructRjWith(m,toUse) 
  length = size(ALLPAULI,2) #number of tensored matrixes.
  dim=convert(Int64,sqrt(length))
  
  Rjs=Array(Complex{Float64},m)

  for i=1:m
    Rjs[i]=trace((reshape(ALLPAULI[:,toUse[i]],dim,dim))*densityMatrix)
  end
  # return the results 

  Rjs=real(Rjs)
  return Rjs
end

function rebuildIt(Rjs,toUse) 
  # calculate the number of PauliMatrices used.
  length = size(toUse,1)
  # the dimensions of the nxn Pj
  dim=convert(Int64,sqrt(size(ALLPAULI,2)))
  # create the empty container
  pHat=zeros(Complex{Float64},dim,dim)
  # reconstruct it.
  for i=1:length
    pHat = pHat + Rjs[i]*reshapedFAllPauli(toUse[i],dim)
  end
  

  # CHECK - BUT INITIALISE NOW SCALES THEM
  # The first pFactor is designed to deal with the fact that the Pauli Matrices weren't 
  # scaled by 1/sqrt(2)
  # pFactor = (1/sqrt(2))^pauliSize
  pFactor = 1

  #The second pfactor is to do with the number of matrices taken?
  # TO CHECK NOT SURE ABOUT THIS - might not be necessary?
  pFactor2 = 1 #sqrt(dim/length)
 
  return pHat * pFactor*pFactor2
end

function funStats(n::Int64,qubits)
  initialise(qubits)
  a=zeros(Float64,size(ALLPAULI,2))
  
  for j=1:n
  initialiseVandDM(qubits)
  for i = 1:size(ALLPAULI,2)
    (R,U) = constructRj(i)
    phHat = rebuildIt(R,U)
    diff = phHat-densityMatrix
    # ret = real(trace(sqrtm(diff'*diff)))
    # ret = .5*sum(svd(diff)[2])
    ret = sqrt(real(trace(diff'*diff)))
    print(".");
    a[i]=a[i]+ret
  end
  println(":",j)
  end
  return a/n
end

function funOldStats(n::Int64,qubits)
  initialise(qubits)
  a=zeros(Float64,size(ALLPAULI,2))
  
  for j=1:n
  initialiseVandDM(qubits)
  for i = 1:size(ALLPAULI,2)
    (R,U) = constructRj(i)
    phHat = rebuildIt(R,U)
    diff = phHat-densityMatrix
    #ret = real(trace(sqrtm(diff'*diff)))
    ret = .5*sum(svd(diff)[2])
    # ret = sqrt(real(trace(diff'*diff)))
    print(".");
    a[i]=a[i]+ret
  end
  println(":",j)
  end
  return a/n
end

function labStats(s1,n::Int64,qubits)
  # Using the matlab interface at s1, calls the script
  # specified in sendItToMatlab to use Matlab
  # to work out the reconstructed density matrix.
  # Averages n times, using a qubits system.
  initialise(qubits)
  a=zeros(Float64,size(ALLPAULI,2))
  println(a)

  for j=1:n
  initialiseVandDM(qubits) #note this randomises the vector.
  for i = 1:size(ALLPAULI,2)
    (R,U) = constructRj(i)
    phHat = sendItToMatlab(s1,R,U)
    diff = phHat-densityMatrix
    ret = sqrt(real(trace(diff'*diff)))
    a[i]=a[i]+ret
    println(a)
  end
  println(":",j)
  end
  return a/n
end




function fakeStats(n::Int64)
  a=zeros(Float64,size(ALLPAULI,2))
  
  fakeNorm = randn(vectorSize)+rand(vectorSize)*im
  fakeNorm = fakeNorm/norm(fakeNorm)
  fakeDensity = fakeNorm*fakeNorm'

  for j=1:n
  for i = 1:size(ALLPAULI,2)
    (R,U) = constructRj(i)
    phHat = rebuildIt(R,U)
    
    diff = phHat-fakeDensity
    ret = real(trace(sqrtm(diff'*diff)))

    print(".");
    a[i]=a[i]+ret
    
  end
  println(":",j)
  end
  return a/n
end

function makeMatlabPauli(toUse)
  # unfortunately the mxarray doesn't seem to send Complex
  # arrays properly, so sending it in two parts, then
  # we recombine in Matlab.

  length = size(toUse,1)

  par = zeros(Float64,vectorSize*vectorSize,length) 
  pai = zeros(Float64,vectorSize*vectorSize,length) 
  for i=1:length 
    par[:,i]=real(ALLPAULI[:,toUse[i]])
    pai[:,i]=imag(ALLPAULI[:,toUse[i]])
  end
  par = transpose(par)
  pai = transpose(pai)
  return(par,pai)
end

function sendItToMatlab(s1,rjs,toUse)
  (pr,pi)=makeMatlabPauli(toUse)
  put_variable(s1,:PR,mxarray(pr))
  put_variable(s1,:PI,mxarray(pi))
  put_variable(s1,:data,mxarray(rjs))
  eval_string(s1,"PAULI=PR+PI*i;")
  #we need the conjugate of it on matlab.
  eval_string(s1,"PAULI=conj(PAULI);")
  #passing dim in an mxarray didn't work , so get matlab to calculate it.
  eval_string(s1,"dim=sqrt(size(PAULI,2));")
  eval_string(s1,"CVXLeastSquares")
  eval_string(s1,"xr=full(real(x));")
  eval_string(s1,"xi=full(imag(x));")
  #again we need to bring it back from matlab in parts.
  xr_m=get_mvariable(s1,:xr)
  xi_m=get_mvariable(s1,:xi)
  xr=jarray(xr_m)
  xi=jarray(xi_m)
  x=xr+xi*im
  return x
end

# gets matlab to reconstruct using
# CVXLestSquares.
# sendItToMatlab does the heavy lifting.
function checkWithMatlab(session,m) 
    (R,U) = constructRj(m)
    phHat = sendItToMatlab(session,R,U)
    diff = phHat-densityMatrix
    return phHat,diff
end

function rollAndTest(session,qubits)
  # save each iteraton seperately incase we want to look at 
  # variance etc or to combine seperate runs.
  maxTests = 2000
  interval = 20
  # I have assumed maxTests%interval = 0.
  numberOfTests = maxTests/interval
  initialise(qubits)

    pSize = size(ALLPAULI,2)
    a = zeros(numberOfTests*pSize,3)
    
    # we need to randomise the vector and density matrix
    # for each iteration
    initialiseVandDM(qubits)
    
    # Iterate over each of our possible maxTests
    # i.e. where we try and work out the probabilities in the experimant
    # with one measurement, two measurements, three etc up to maxTests.

    # now I need to remember to store it in x y z format for easy
    # splot ting.

    for rolls = interval:interval:maxTests #obviously we don't start at 0 or 1.
      # Iterate over a possible selection of 1 to all the Paulis.
      for i = 1:pSize
          # so for the Paulis we choose...
          print(i,":")
          (R,U) = constructRj(i)
          
          # Pj will hold the probability of getting a +1 on a measurable p(+)
          Pj=zeros(size(R,1))
          
          # rR is the randomly constructed R based on rolls .v. the probability.
          (rR,Pj) = rollNewR(R,rolls) 
          # println(R-rR) - sanity check print, seems to work.
          # now that we know the estimated values, 
          # do something to construct our guess of the density matrix
          # based on pseudo measurements.

          phHat = sendItToMatlab(session,rR,U)
          print(".")

          diff = phHat-densityMatrix
          ret = sqrt(real(trace(diff'*diff)))
          # save the results in a.
          index = (pSize*(rolls-interval)/interval)+i
          a[index,1]= rolls
          a[index,2]= i
          a[index,3] = ret
      end # of for each of our pauliSize
      println("")
   end # of iterating over the number of experiments performed.
   
  return a
end # of function.

function doAfew(session,qubits)
  #writing them out both comma seperated and (for gnuplot) space delimited.
  println("Starting first")
  a = rollResourcesAndTest(session,qubits)
  writecsv("LR800x5Q4Runa1.csv",a)
  writedlm("LR800x5Q4Runa1.ssv",a,' ')
  a = rollResourcesAndTest(session,qubits)
  writecsv("LR800x5Q4Runa2.csv",a)
  writedlm("LR800x5Q4Runa2.ssv",a,' ')
  a = rollResourcesAndTest(session,qubits)
  writecsv("LR800x5Q4Runa3.csv",a)
  writedlm("LR800x5Q4Runa3.ssv",a,' ')
  a = rollResourcesAndTest(session,qubits)
  writecsv("LR800x5Q4Runa4.csv",a)
  writedlm("LR800x5Q4Runa4.ssv",a,' ')
end

function rollResourcesAndTest(session,qubits)
  # save each iteraton seperately incase we want to look at 
  # variance etc or to combine seperate runs.
  # I have assumed maxTests%interval = 0.
  # also saving in a format for gnuplot i.e. x, y, z
  numberOfTests=800
  initialise(qubits)
  pSize = size(ALLPAULI,2)
  a = zeros(numberOfTests*pSize,3)
    
  # we need to randomise the vector and density matrix
  # for each iteration
    
    # Iterate over each of our possible maxTests
    # i.e. where we try and work out the probabilities in the experimant
    # with one measurement, two measurements, three etc up to maxTests.

    # now I need to remember to store it in x y z format for easy
    # splot ting.

    for resources = 1:numberOfTests #obviously we don't start at 0 
      initialiseVandDM(qubits)
      # Iterate over a possible selection of 1 to all the Paulis.
      for i = 1:pSize
          # so for the Paulis we choose...
          print(i,":")
          (R,U) = constructRj(i)
          
          # Pj will hold the probability of getting a +1 on a measurable p(+)
          Pj=zeros(size(R,1))
          
          # need to calculate number of rolls given limited resources
          rollsToDistribute = resources * 5 * pSize
          baseRolls = div(rollsToDistribute,i)
          rolls=ones(Int64,i)*baseRolls
          for extra=1:(rollsToDistribute%i) 
            rolls[extra]+=1
          end

          # rR is the randomly constructed R based on rolls .v. the probability.
          (rR,Pj) = rollNewR(R,rolls) 
          # println(R-rR) - sanity check print, seems to work.
          # now that we know the estimated values, 
          # do something to construct our guess of the density matrix
          # based on pseudo measurements.

          phHat = sendItToMatlab(session,rR,U)
          print(".")

          diff = phHat-densityMatrix
          ret = sqrt(real(trace(diff'*diff)))
          # save the results in a.
          index = (pSize*(resources-1))+i
          a[index,1]= resources
          a[index,2]= i
          a[index,3] = ret
      end # of for each of our pauliSize
      println("")
   end # of iterating over the number of experiments performed.
   
  return a
end # of function.

function rollNewR(R,numberOfRolls)
      # Iterate over a possible selection of 1 to all the Paulis.
      # Pj will hold the probability of getting a +1 on a measurable p(+)
      Pj=zeros(size(R,1))
      
      # rR is the randomly constructed R based on rolls .v. the probability.
      rR = zeros(size(R,1))
      

      
          # R is the expectation value, use it to construct the probabilities i.e.  
          # p(+) = (expectation_value+1)/2
          # then for each of the measurements we are doing
          # roll randomly against this probability to simulate
          # the results of "rolls" number of tests.
          for j=1:size(R,1)
            Pj[j]=(R[j]+1)/2 #probability of +1

            # conduct the number of tests we are doing and remember if we measured
            # +1 or -1.
            for doingRoll=1:numberOfRolls[j]
              if rand() <= Pj[j] 
                rR[j]+=1
              else 
                rR[j]-=1
              end
            end
            # then average the result, to achieve our guess of the probability.
            rR[j]=rR[j]/numberOfRolls[j]
          end
         
  return rR,Pj
end # of function.
#crib sheet s1=MSession()
#close(s1)


