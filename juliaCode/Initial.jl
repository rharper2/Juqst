# This file contains all the basic functions needed to run a quantum simulator
# using the stabiliser formalism
# The format used is as described in 
# Improved Simulation of Stabilizer Circuits
# Scott Aaronson and Daniel Gottesman
# arXiv:quant-ph/0406196v5

# It is, effectively, a Julia Port of the CHP program by Scott Aaronson
# The original CHP program can be found at http://www.scottaaronson.com/chp

# Primarily, I have tried to write for clarity rather than efficiency.
# (Although thanks to the algorithms by Aaronson/Gottesman its fast enough for hundreds of qubits)

# I am partially through rationalising some things for ease of REPL
# If you use a ! command like setup!(4), hadamard!(2) then it acts on the global state
# Alternatively you can pass in a state like this
# mystate = setup(3)
# cnot(mystate,1,2)
# In which case it will act on the supplied state.
# Note that this isn't quite in accord with Julia syntax since we still alter the supplied state
# But there you go.



function setup!(n)
  global state
  state = setup(n)
end

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

function output()
  global state 
  output(state)
end

function output(state)
	# prints the Heisnberg representation of the state array supplied
	n=div(size(state,1),2)
	#  print("Size of system (N) is ",n,"qubits\n\n")
	#  1 in a and a+n = Y
  #  else 1 in a = X
	#  else 1 in a+n=Z
	#  otherwise I
	# last bit is sign.
	for row=1:n
	  if state[row,2*n+1] == 1 
		print("-")
	  else
		print("+")
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
	  if state[row,2*n+1] == 1 
		  print("-")
	  else
		  print("+")
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


function j_xor(a,b)
  return mod.(a+b,2)
end


function cnot!(a,b,showOutput=true)
   global state;
  cnot(state,a,b,showOutput)
end

function cnot(state,a,b,showOutput=true)
  # From control a to taget b
  # You can supress the output by calling with the fourth parameter false
  # otherwise it will "output" the stabiliser state.
  # note the state is changed, (in Julia this should really be cnot!)
  n=div(size(state,1),2)
  endC=2*n+1

  if a < 1 || a > n
    print("Control qubit a(",a,") out of range of state which has ",n," qubits\n")
    throw(BoundsError())
    return
   end
  if b < 1 || b > n
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
      state[i,endC] = j_xor(ri,xia*zib*(j_xor(xib,j_xor(zia,1))))
      state[i,b]=j_xor(xib,xia)
      state[i,a+n]=j_xor(zia,zib)
      #print("a is $(a) b is $(b), xia is $(xia) xib is $(xib), setting $([i,b]) to $(state[i,b])\n")
   end
   if showOutput
      output(state)
    end
end


function hadamard!(a,showOutput=true)
  global state;
  hadamard(state,a,showOutput)
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
      state[i,endC]=j_xor(ri,xia*zia)
      state[i,a]=zia
      state[i,a+n]=xia
   end
   if showOutput
    output(state)
  end
end


function phase!(a,showOutput=true)
  global state
  phase(state,a,showOutput)
end

function phase(state,a,showOutput=true)
  # Apply the phase gate to qubit a
  # again changes the state (phase!)

  n=div(size(state,1),2)
  endC=2*n+1

  if a < 0 || a > n
    print("Qubit a(",a,") out of range of state which has ",n," qubits\n")
    return
   end
  # so we need two loops!
  for i = 1:(2*n)
      ri=state[i,endC]
      xia = state[i,a]
      zia = state[i,a+n]
      state[i,endC]=j_xor(ri,xia*zia)
   end
   for i = 1:(2*n)
      ri=state[i,endC]
      xia = state[i,a]
      zia = state[i,a+n]
      state[i,a+n]=j_xor(xia,zia)
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
  n=div(size(state,1),2)
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
	state[h,j]=j_xor(state[i,j],state[h,j]) # xhj=xij j_xor xhj
	state[h,n+j]=j_xor(state[i,n+j],state[h,n+j])
  end
end
	
function measure!(a)
  global state 
  return measure(state,a)
end 


function measure(state,a)
  n=div(size(state,1),2)
  endC=2*n+1
  #print("Before\n")
  #output(state)
  #print("\n")
  if a < 0 || a > n
    print("Qubit a(",a,") out of range of state which has ",n," qubits\n")
    return
   end

  aSlice=state[(n+1):end,:] # so aSlice is a "slice" of the tableau, being the commuting part.
 # print("$(aSlice)\n$(aSlice[:,a])\n")
  if sum(aSlice[:,a]) > 0 
	   # we have a one on the xia somewhere. This makes it a random measuremnt. 
  	p = n+indmax(aSlice[:,a]); # first one will do
	  # print("Random first x row is ",p,"\n\n")
  	# First rowsum(state,i,p) for all i 1..2n where i!=p and xia =1
   	for h=1:2n # h is what we use in rowsum, i is p in rowsum 
	     if (h!=p) && (state[h,a]==1)
		      rowsum(state,h,p)
	     end
	  end
	  # Second step set p-nth to equal p
	  state[p-n,:]=state[p,:]
    # Third step set pth row to be 0, except rp = 0 or 1 and zpa = 1
	  value = round.(Integer,rand()) # check distribution
	  state[p,endC]=value
	  state[p,:]=zeros(2*n+1)
	  state[p,endC]=value
	  state[p,n+a]=1
   # output(state)
	 # print("\n")
    return value
  else
      # print("Deterministic\n\n")
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
      # print("After\n")
      # output(state)
      # print("\n")
       return value
   end
end

function rowswap!(alteredState,i,k)
      temp = alteredState[i,:]
      alteredState[i,:]=alteredState[k,:]
      alteredState[k,:]=temp
end

function cliffordPhase(i,k)
  # Returns the phase (0,2,3,4) when row i is left multiplied by row k
  # println("CliffordPhase $i,$k")
  e=0
  n=div(length(i),2)
  #println("n is $n")
  for j=1:n
    if k[j] == 1 && k[j+n]==0 # its an X
      if i[j] == 1 && i[j+n] == 1 # its a Y
        #print("XY\n")
        e+=1 # XY=iZ
      elseif   i[j] == 0 && i[j+n] == 1 # its a Z
        #print("XZ\n")
        e+=-1 # XY=-iZ
      end
    elseif  k[j] == 1 && k[j+n]==1 # its a Y
      if i[j] == 0 && i[j+n] == 1       # its a Z
        #print("YZ\n")
        e+=1 # YZ=iXZ
      elseif i[j] == 1 && i[j+n] == 0 # its a X
        #print("YX\n")
        e+=-1 # YX=-iZ
      end
    elseif  k[j] == 0 && k[j+n]==1 # its a Z
      if i[j] == 1 && i[j+n] == 0       # its a X
        #print("ZX\n")
        e+=1 # ZX=iY
      elseif i[j] == 1 && i[j+n] == 1 # its a Y
        #print("ZY\n")
        e+=-1 # ZY=-iX
      end
    end
  end

 # NOTE the 2* in front of k but not i
 # This is a TERRIBLE hack
 # See rowmult! for the full grimy disclosure.
  e = e + i[2*n+1] + 2*k[2*n+1]
   #print("E is now $e\n")
  e = e%4
  if (e < 0) 
    e+=4
  end
   #print("modded E is now $e\n")
  
  return e
end 


# This one is called by measure
function rowmult!(alteredState,i,k) # does the multiplication in place in a state array.
  # left multiply row i by row k
  # println("Rowmult! with $i and $k")
  n=div(size(alteredState,1),2)
  

  # This is part of a TERRIBLE hack.
  # For the purpose of the *tableau* we only need to record an r which is 0 or 1 (in place 2*n+1)
  # And so rather than allow it four values, the software was coded using a 0 or 1
  # This makes matching the stuff on page 4 column 2 to the code (imo) easier, it directly matches the maths.
  # HOWEVER, when pritinging out the basis (getStatesFor) we need all four possible states +1,i,-1,-i
  # So we can't just divide the Clifford phase by 2.
  # All the syplectic stuff built on top of this assumes r = 0 or 1, so I don't want to rever to the Aaronson code
  # which has r ranging from 0->3 (and divides by 2 normally)
  # But afaik this only impact when we getStatesFor (to print out the basis states).

  # So.. (the HACK) ... when you call cliffordPhase it assumes the first row has an r that has been multiplied by 2
  # i.e. can be between 0..3, but the second will be either 0 or 1. That is how rowmult works (basically it has a different
  # scratch state and we 'reinterpret' the r so it can 0..3, and printBasisState deals with this.)
  # BUT here we are only interested in 0 or 1, but we have to multiply the r for the first state passed into cliffordPhase by 2
  # (i.e. convert to the 0..3 scale) and then divide the result by 2, to get from the 0..3 scale to the 0,1 scale.
  # as I said (terrible hack).

  
  # Multiply this one by two (0..3) scale
  alteredState[i,2*n+1] = alteredState[i,2*n+1]*2
 #   println("alteredState : $(alteredState[i,:])")
 #   println("secondState : $(alteredState[k,:])")

  phaseBack = cliffordPhase(alteredState[i,:],alteredState[k,:])
 #    println("Clifford phase returned $(phaseBack)")
  if (phaseBack  == 1 || phaseBack  ==3)
    println("phaseBack  WAS 1 OR 3")
  end
  # Divide answer by 2 (from 0..3->0..1 scale.)
  # alteredState[i,2*n+1] = alteredState[i,2*n+1]/2
#  println("Changed to $(phaseBack/2)")
  
 #println("Setting altered state to $(alteredState[i,2*n+1])")
  #println("back")
  alteredState[i,:] = xor.(alteredState[i,:],alteredState[k,:])
  # we dont xor the r bits.
  alteredState[i,2*n+1] = phaseBack/2
# println("Returning an alterestate of $(alteredState[i,:])")
end


# This one is called by getStatesFor (which prints out the relevant basis states).
function rowmult(i,k) # supply two vectors and get the multiplied vector out.
  #println("Rowmult with $i and $k")
  n=div(length(i),2)
  #println("i is $i\nk is $k\n n is $n\n")
  temp = cliffordPhase(i,k)
  #if (temp == 1 || temp ==3)
  #  println("TEMP WAS 1 OR 3")
  #end
  #println("$i and $k and $(xor.(i,k))")
  tempState = [xor(i[c],k[c]) for c =1:length(i)] # $ is Julia's bitwise j_xor.
  tempState[length(i)]=temp
  #println("Returning a temp of $temp\n")
  # println("tempstate is $tempState")
  #  tempState[2*n+1,:]=temp/2
  return tempState
 end

#returns a tableau and a size where gaussian elimination has
#been carried out on the supplied state, in order to
#provide a minimal set of generators containing X's and Y's in 
#upper triangular form 

function gaussianElimination(state)
  alteredState = copy(state)
  #println(alteredState)
  n=div(size(alteredState,1),2)
  i = n+1 # first row of commuting generators
  k=0
  k2=0
  j=0
  for j=1:n
    found = 0;
  
    for k=i:2*n
      # Find a generator containing X in the jth column
      if alteredState[k,j] == 1
        found = 1
        break
      end
    end
    if found  == 1
      #swap row with the row pointed to by i
      if (i!=k) 
        rowswap!(alteredState,i,k)
        rowswap!(alteredState,i-n,k-n) # swap the non-commutators as well
      end
      # then use the row to eliminate X from that bit in the rest of the tableau
      for k2=i+1:2*n
        if alteredState[k2,j] == 1
          #println("Before, $(k2) and j is $j and i is $i")
          #println(alteredState)
          #print("ROWMULTA $k,$i\n")
          rowmult!(alteredState,k2,i)
          #println(alteredState)
          #println("Doing $(i-n) and $(k2-n)")
          #print("ROWMULTB $(i-n) $(k2-n)\n")
          rowmult!(alteredState,i-n,k2-n)
          #println("Row mult $k2, $i $alteredState")
        end
        #println(alteredState)
      end
      #println(alteredState)
      i+=1
    end
  end
  gen = i - n
  # The first gen generators are X/Ys and in quasi upper triangular form.
  for j=1:n
    found = 0
    for k=i:2*n
      if alteredState[k,n+j] == 1 # we have found a Z in the jth bit
        found = 1
        break
      end
    end
    if found == 1
      if (i!=k)
        rowswap!(alteredState,i,k)
        rowswap!(alteredState,i-n,k-n)
      end
      for k2=i+1:2*n
        if (alteredState[k2,j+n] ==1 ) # z in this "bit"
          #print("ROWMULTC $k2,$i\n")
         
          rowmult!(alteredState,k2,i)
          #print("ROWMULTD $(i-n),$(k2-n)\n")
         
          rowmult!(alteredState,i-n,k2-n)
        end
      end
      i+=1
    end
    
  end
  return (alteredState,gen)
end


function printBasisState(sS)
    nqubits=div(size(sS,1),2)
    #println("\npBS: $sS and $nqubits and $(size(sS)) e=$(sS[2*nqubits+1])\n")
    e=sS[2*nqubits+1]
    for i=1:nqubits
      if sS[i]==1 && sS[nqubits+i]==1 # its a Y
        e = (e+1) % 4
      end
    end
    if (e==0) print("+|")
    elseif e== 1 print("+i|")
    elseif e==2 print("-|")
    else print("-i|")
    end
    for i=1:nqubits
      if sS[i]==1
        print("1")
      else 
        print("0")
      end
    end
    print(">\n")
end  

function seed(state,n)
  nqubits=div(size(state,1),2)
  #println("\n\nQubits $nqubits going $(2*nqubits) to $(nqubits+n)\n")
  scratchState = zeros(Int32,1,2*nqubits+1)
   for i=2*nqubits:-1:(nqubits+n)
    min = nqubits
    f = 2*state[i,2*nqubits+1]   
    #println("\nF is $f\n\n");
    # count down through the qubits.
    for j=nqubits:-1:1
      # if the z part of this is 1
      #print("Checking bit $j of string \n$(state[i,:])\n")
      if state[i,nqubits+j] == 1
        min = j
        #print("Setting minimum to $min\n");
        if scratchState[1,j]==1 
          f = (f+2)%4
          #print("Increasing f to $f");
        end
      end
    end
    if f==2
      #print("F was 2 so we need to flip the x");
      #print("Scratch was $(scratchState)\n")
      scratchState[1,min] = xor(scratchState[1,min],1)
      #print("Now $(scratchState)\n")
      
    end
  end
  #print("Freturning $scratchState\n")
  return scratchState
end

function getStates()
  global state
  getStatesFor(state)
end 


function getStatesFor(originalState)
    ss1 = gaussianElimination(originalState)
    # print(ss1);
   # print("Working on $(ss1[2])\n")
    #output(ss1[1])
    #print("=====$(ss1[2])\n")
    state=ss1[1]
    # output(state)
    #println(ss1[1])
    n=ss1[2]
    #print("n is $n")
    start=seed(state,n)
    scratchState=copy(start) # the first nqubits are the destabilisers
    printBasisState(vec(scratchState))   
    # So the idea is we have the stabilisers with an X in them in upper triangular form
    # Each X represents a "1" in the Ket
    # But we need to generate each of the states stabilised by these Paulis.
    # For the first row, the ones generated are where there is an X in the first row.
    # On the second row, the ones generated are the second row, and the second row multiplied by the first.
    # On the third row, its the third row, third times first, third times second and third times first and second
    # etc.
    # Aaronsen has a really neat algorighm to do this, but its not easily understood.
    # I have used this, slighly less efficient, but (I hope) more comprehenisible algorithm
    # For a particular row (say we are looking at our third row with an X in it (n 3) 
    # loop through the numbers 4:7 (being 2^(3-1) to (2^3)-1)
    # In binary, this loops as 100, 101, 110, 111, So we can see This is applying
    # Row 3 with each possible combination of Rows 1 and 2.
    # 
    # Then we need to multiply the rows together if 2^(row-1) (i.e. row 2 - equates to bit 2) is not zero when
    # logically anded with the count.
    number=2^(n-1)
    nqubits=div(size(state,1),2)
    #println("nqubits is $nqubits")
    for count = 0:(2^(n-1)-2)
        c2=xor(count,(count+1))
        #println("count is now $count c2 is $c2\n")
        for rows=0:(n-2)
            #println("count: $count, c2 is $c2, rows is now $rows\n")
            if (2^(rows) & c2) > 0 
               #println("******Row mult \nscratchState is $scratchState\n State: $(state[nqubits+rows+1,:])\n")
               scratchState=rowmult(scratchState,state[nqubits+rows+1,:])
               #println("\n**scratchState is now $scratchState\n\n")
            end
        end
       # println("Done")
        # Note the use of vec here, just makes sure that the scratch state is a vectdor
        # This was only an issue for the |00> one.
        printBasisState(vec(scratchState))
    end
end

