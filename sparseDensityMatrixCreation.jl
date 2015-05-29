#Generates the Base Paulis in a 4x4 (each Pauli is vectorised)
#This generates them in sparse matrices.
#This is slower than the dense ones, but might be necessary for
#memory requirements?

function generateSparseBasePaulis()
  global sPi,sPx,sPz,sPy
  sPi = sparse([1,2],[1,2],[1,1])
  sPx = sparse([1,2],[2,1],[1,1]) 
  sPz = sparse([1,2],[1,2],[1,-1]) 
  sPy = sparse([1,2],[2,1],[-1*im,1*im]) 
  global sparsePaulis=[sPi,sPx,sPy,sPz]
end

#decodes the vectorised Pauli to return it as a 2x2 sparse matrix.

function getPauli(n::Int64)
 if (n<1) || n>4
    return "Undefined"
 end
 index = n*2-1
 return sparsePaulis[index:index+1,1:2]
end

# SPARSE version
# Using a recursive function here to construct all possible {P1..P4} tensor {P1..P4} tensor {P1..P4}
# Pauli functions where d is the number of tensor functions to carry out.
#
# The returned matrix has each combination vectorised
# and returned as a column of the matrix. The matrix will have 4^d columns and 4^d rows.
#
# This keeps everything in sparse arrays, its a lot slower than the dense version but (presumably) uses
# less memory.
# Note still missing the factor of (1/sqrt(2))^d

# Note this will throw errors unless we are using the "beta" version of julia.
# The latest stable version has a problem with converting to a sparse array of complex numbers.

function getSparsePaulis(d)
  println("Get Sparse Paulis ",d)
# d should never be less than one.
  if d<1 return 1 end
# if d is one, we just need to return each of the Paul (2x2) matrices vectorised.
  if d==1 
    return [sparse(getPauli(1)[:]) sparse(getPauli(2)[:]) sparse(getPauli(3)[:]) sparse(getPauli(4)[:])]
  end

#otherwise get the result of all the paulis for one dimension down.
#we need to return each column of that matrix (resized to a square matrix)
#tensored with each of the four Pauli matrices.
  
  temp = getSparsePaulis(d-1)
  length=size(temp,1)
  width=size(temp,2)
  println("Width",width)
  dim=convert(Int64,sqrt(length))
# directly creating the sparse array was giving us problems, so create a 
# dense array which is converted to sparse
# its not in the loop so the speed loss isn't important.
  answerSoFar=Array(eltype(temp),4^d,4^d)
  answerSoFar=sparse(answerSoFar)
  #loop over each column and combine with each of the four paulis.
  for i=1:width
    for j=1:4
        # using a temporary variable here for readability
        # this takes the ith column and reshapes it to a dimxdim matrix
        t1=reshape(temp[:,i],dim,dim)
        #we then want to take the jth Pauli, tensor multiply with t1
        #vectorise the resultant matrix and add it on as a new column.
        #println("Setting to",(kron(getPauli(j),t1))[:])
        answerSoFar[:,(i-1)*4+j] = sparse((kron(getPauli(j),t1))[:])
    end
  end
  
   println("Get All Paulis completed for ",d)
   return answerSoFar
end

