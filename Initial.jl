# This file contains all the basic functions needed to run a quantum simulator
# using the stabiliser formalism
# The format used is as described in 
# Improved Simulation of Stabilizer Circuits
# Scott Aaronson and Daniel Gottesman
# arXiv:quant-ph/0406196v5

# It is, effectively, a Julia Port of the CHP program by Scott Aaronson
# The original CHP program can be found at http://www.scottaaronson.com/chp

# I have tried to write for clarity rather than efficiency.
# (Although thanks to the algorithms by Aaronson/Gottesman its fast enough for hundreds of qubits)

function setup(n)
  # sets up and returns the initial state vector for n qubits.
  # Heisenberg representation
  # We are representing a |00...0> ket in the stabiliser state
  # This stabilisises with "Z" and anti-stabilises with "X"
  # So the tableau starts off with an Identity in the Anti commute top X section
  # and identity in the commuting - Z section.

  # For the purposes of this port, the tableau is exactly replicated as drawn
  # i.e. the "state" is an Int32 array (should really be just binary)
  # of the form 
  #
  #
  # x11   .....  x1n | z11   ...      z1n | r1
  #  .    \       .  |  .    \          . |  .   
  #  .     \      .  |  .     \         . |  .     Destabilisers
  #  .      \     .  |  .      \        . |  .
  # xn1      \   xnn | zn1      \      znn| rn 
  # ______________________________________________
  # x(n+1)1. x(n+1)n | z(n+1) ... z(n+1)n | r(n+1)
  #  .    \       .  |  .      \        . |  .   
  #  .     \      .  |  .       \       . |  .     Stabilisers
  #  .      \     .  |  .        \      . |  .
  # x(2n)1   \x(2n)n | z(2n)1     \ z(2n)n| r(2n)
  
  top = hcat(eye(Int32,n),zeros(Int32,n,n),zeros(Int32,n,1))
  bottom = hcat(zeros(Int32,n,n),eye(Int32,n),zeros(Int32,n,1))
  state = vcat(top,bottom)
end

function output(state)
	# prints the Heisnberg representation of the state array supplied
	n=div(size(state)[1],2)
	#  print("Size of system (N) is ",n,"qubits\n\n")
	#  1 in a and a+n = Y
  #  else 1 in a = X
	#  else 1 in a+n=Z
	#  otherwise I
	# last bit is sign.
	for row=1:n
	  if state[row,2*n+1] == 0 
		print("+")
	  else
		print("-")
	  end
	  for col=1:n
	     if state[row,col] == 1
		      if state[row,col+n] == 1
			       print("Y")
		      else
			       print("X")
	        end
	     elseif state[row,n+col]==1
		      print("Z")
	     else
		      print("I")
	     end
	  end
	  print("\n")
	end
	for dot=1:n+1
		print("-")
	end
	print("\n")
	for row=n+1:2*n
	  if state[row,2*n+1] == 0 
		  print("+")
	  else
		  print("-")
	  end
	  for col=1:n
	     if state[row,col] == 1
		      if state[row,col+n] == 1
			       print("Y")
		      else
			       print("X")
	        end
	     elseif state[row,n+col]==1
		      print("Z")
	     else
		      print("I")
	     end
	  end
	  print("\n")
	end
end

function xor(a,b)
  return mod(a+b,2)
end




function cnot(state,a,b,showOutput=true)
  # From control a to taget b
  # You can supress the output by calling with the fourth parameter false
  # otherwise it will "output" the stabiliser state.
  # note the state is changed, (in Julia this should really be cnot!)
  n=div(size(state)[1],2)
  endC=2*n+1

  if a < 0 || a > n
    print("Control qubit a(",a,") out of range of state which has ",n," qubits\n")
    throw(BoundsError())
    return
   end
  if b < 0 || b > n
    print("Target qubit b(",b,") out of range of state which has ",n," qubits\n")
    throw(BoundsError())
    return
   end
  for i = 1:(2*n)
      ri=state[i,endC]
      xia = state[i,a]
      xib = state[i,b]
      zia = state[i,a+n]
      zib = state[i,b+n]
      state[i,endC] = xor(ri,xia*zib*(xor(xib,xor(zia,1))))
      state[i,b]=xor(xib,xia)
      state[i,a+n]=xor(zia,zib)
   end
   if showOutput
      output(state)
    end
end



function hadamard(state,a,showOutput=true)
  # Apply a hadamard gate to quibit a in the state, state.
  # again note it changes the state, so it should be hadamard!

  n=div(size(state)[1],2)
  endC=2*n+1

  if a < 0 || a > n
    print("Qubit a(",a,") out of range of state which has ",n," qubits\n")
    return
   end
  for i = 1:(2*n)
      ri=state[i,endC]
      xia = state[i,a]
      zia = state[i,a+n]
      state[i,endC]=xor(ri,xia*zia)
      state[i,a]=zia
      state[i,a+n]=xia
   end
   if showOutput
    output(state)
  end
end

function phase(state,a,showOutput=true)
  # Apply the phase gate to qubit a
  # again changes the state (phase!)

  n=div(size(state)[1],2)
  endC=2*n+1

  if a < 0 || a > n
    print("Qubit a(",a,") out of range of state which has ",n," qubits\n")
    return
   end
  for i = 1:(2*n)
      ri=state[i,endC]
      xia = state[i,a]
      zia = state[i,a+n]
      state[i,endC]=xor(ri,xia*zia)
      state[i,a+n]=xor(xia,zia)
   end
   if showOutput
    output(state)
  end
end   


function GottesmanG(x1,z1,x2,z2)
# x1z1 and x2z2 represent pauli matrices by the usual manner 
# e.g. both are 1 = Y gate, both 0 = I, otherwise same as the one thats 1.
# this returns the exponent to which i is raised if we multiply x1z1 by x2z2
 
  if x1==0
    if z1 == 0 # x1=0 z1=0, so multiplying by I
       return 0
    else  # x1=0 z1=1, so multiplying by Z
       return x2*(1-2*z2)
    end
  else
    if z1 == 0 # x1=1, z1=0, so multiplying by X
       return z2*(2*x2-1)
    else # x1=1 z1=1
       return z2-x2
    end
  end
end


function rowsum(state,h,i)
  # "sums" two rows in the state matrix
  # pass in integers to the rows, not the rows themselves.
  n=div(size(state)[1],2)
  total = 0;
  for j=1:n
	xij=state[i,j]
        zij=state[i,n+j]
        xhj=state[h,j]
	zhj=state[h,n+j]
	total += GottesmanG(xij,zij,xhj,zhj)
  end
  rh=state[h,2*n+1]
  ri=state[i,2*n+1]
  grandTotal=2*rh+2*ri+total
  if mod(grandTotal,4) == 0
	state[h,2*n+1]=0
  elseif mod(grandTotal,4) == 2
	state[h,2*n+1]=1
  else
	print("Failed sanity check the grandTotal ",grandTotal, " which mods to ",mod(grandTotal,4),"\n")
  	return
  end
  for j=1:n
	state[h,j]=xor(state[i,j],state[h,j]) # xhj=xij xor xhj
	state[h,n+j]=xor(state[i,n+j],state[h,n+j])
  end
end
	
function measure(state,a)
  n=div(size(state)[1],2)
  endC=2*n+1
  print("Before\n")
  output(state)
  print("\n")
  if a < 0 || a > n
    print("Qubit a(",a,") out of range of state which has ",n," qubits\n")
    return
   end

  aSlice=state[(n+1):end,:] # so aSlice is a "slice" of the tableau, being the commuting part.
  if sum(aSlice[:,a]) > 0 
	   # we have a one on the xia somewhere. This makes it a random measuremnt. 
  	p = n+indmax(aSlice[:,a]); # first one will do
	  print("Random first x row is ",p,"\n\n")
  	 # First rowsum(state,i,p) for all i 1..2n where i!=p and xia =1
   	for h=1:2n # h is what we use in rowsum, i is p in rowsum 
	     if (h!=p) && (state[h,a]==1)
		      rowsum(state,h,p)
	     end
	  end
	  # Second step set p-nth to equal p
	  state[p-n,:]=state[p,:]
    # Third step set pth row to be 0, except rp = 0 or 1 and zpa = 1
	  value = int(rand()) # check distribution
	  state[p,endC]=value
	  state[p,:]=zeros(2*n+1)
	  state[p,endC]=value
	  state[p,n+a]=1
    output(state)
	  print("\n")
    return value
  else
       print("Deterministic\n\n")
       #first set the 2n+1 row i.e. scratch space to be zero
       scratch=zeros(Int32,2n+1)' # hey 2n works dont you know
       state = vcat(state,scratch) # add it on.
       #second call rowsum(2n+1,i+n) for all i in 1..n where xia = 1
       for i = 1:n
	       if state[i,a]==1
	         rowsum(state,2n+1,i+n)
         end
       end
       value = state[2n+1,endC]
       state=state[1:end-1,:]
       print("After\n")
       output(state)
       print("\n")
       return value
   end
end



