
#NOTE Y WAS WRONG ITS BEEN FIXED
#As per generateBasePaulis, but keeps them as dense matrices
#note we are missing the 1/sqrt(2) factor - I am keeping it
#as Complex{Int64} just now as its easier to check its working.
#and Ill multiply by (1/sqrt(2))^n later (where n = number of qubits)
# moved to a more natural way of instantiating them.

function generateFullBasePaulis()
  i = [1 0;0 1] #reshape([1,0,0,1],2,2)
  x = [0 1;1 0] #reshape([0,1,1,0],2,2)
  z = [1 0;0 -1] #reshape([1,0,0,-1],2,2)
  y = [0 -1*im;im 0] #reshape([0,-1*im,1*im,0],2,2)
  global FullBasePaulis=[i,x,y,z]
end

function getPauliName(i)
  if i==1 return "I"
  elseif (i==2) return "X"
  elseif (i==3) return "Y"
  else return "Z"
  end
end


#Decodes the vectorised Pauli to return it as a 2x2 dense matrix.
#Note this is quick its just playing with indices, the matrix
#we get back is the same matrix, just differently indexed.

function getFullPauli(n::Int64)
 if (n<1) || n>4
    return "Undefined"
 end
 index = n*2-1
 return FullBasePaulis[index:index+1,1:2]
end




# DENSE VERSION
# using a recursive function here to construct all possible {P1..P4} tensor {P1..P4} tensor {P1..P4}
# Pauli functions where d is the number of tensor functions to carry out.
#
# the returned matrix has each comination vectorised
# and returned as a column of the matrix. The matrix will have 4^d columns and 4^d rows.
# This is the dense version and is the quicker way of doing it.

# Note still missing factor of (1/sqrt(2))^d

function getAllPaulis(d)
  global ALLPAULI
  println("Get All Paulis ",d)
# d should never be less than one.
  if d<1 return 1 end

# if d is one, we just need to return each of the Paul (2x2) matrices vectorised.
  if d==1 
     return [(getFullPauli(1)[:]) (getFullPauli(2)[:]) (getFullPauli(3)[:]) (getFullPauli(4)[:])], ["I","X","Y","Z"]
  end

# otherwise get the result of all the paulis for one dimension down.
# we need to return each column of that matrix (resized to a square matrix)
# tensored with each of the four Pauli matrices.

  (temp,names) = getAllPaulis(d-1)
  length=size(temp,1)
  width=size(temp,2)
  dim=convert(Int64,sqrt(length))

# create an empty array to hold the tensored product.
  answerSoFar=Array(eltype(temp),4^d,4^d)
  namesSoFar = Array(ASCIIString,4^d)

#loop over each column and combine with each of the four paulis.

  for i=1:width
    for j=1:4
        # using a temporary variable here for readability
        # this takes the ith column and reshapes it to a dimxdim matrix
        # its just indexes I think, so its cheap.

        t1=reshape(temp[:,i],dim,dim)
        # we then want to take the jth Pauli, tensor multiply with t1
        # vectorise the resultant matrix and store it as a column.
        # println("Setting to",(kron(getFullPauli(j),t1))[:])

        answerSoFar[:,(i-1)*4+j] = (kron(getFullPauli(j),t1))[:]
        namesSoFar[(i-1)*4+j] = string(getPauliName(j),names[i])
    end
  end
  # stops me from thinking it crashed...
   println("Get All Paulis completed for ",d)
   ALLPAULI=(answerSoFar, namesSoFar)
  return answerSoFar, namesSoFar
end

# convience function.
# index of AllPauli and the dimensions we need - no error checking
function reshapedAllPauli(n) 
  global ALLPAULI
  dim=int(sqrt(size(ALLPAULI[1])[1]))
  reshape((ALLPAULI[1])[:,n],dim,dim)
  end


#return the reshaped Paul as a Complex{Float64}
# correctly normalised.
function reshapedNormalPauli(n) # index of AllPauli and the dimensions we need - no error checking
  global ALLPAULI
  dim=int(sqrt(size(ALLPAULI[1])[1]))
  (1/sqrt(dim))*reshape((ALLPAULI[1])[:,n],dim,dim) + zeros(Complex{Float64},dim,dim)
end
