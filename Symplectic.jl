# This julia file builds on the Initial.julia
# It is based on the "How to efficiently select an arbitrary Clifford group element"
# Robert Koenig and John Smolin
# arXivv:1406.2170v1

# The premise here is that in order to generate a random Clifford, we can be sure of 
# proper (Haar) randomness if either we generate all the cliffords for a certain quibit size
# and select one randomly OR if we can have a 1-1 mapping between the integers and the cliffords
# and then we can just randomly select an integer. The latter is the one adopted in Koenig's paper.
# getNumberOfCliffords(n) returns the number of cliffords with n qubits. 



# This ties into the Aaronsen/Gottesman paper becauase the generated cliffords are generated 
# in such a way it is possible to create the tableau that would result IF the start myVector
# |000000> had been acted on by the chosen clifford unitary.
# We can then decompose the tableau using the algorithms specified in the Aaronson/Gottesman paper
# (see Initial.jl) to work out what combination of one-qubit (phase and hadmard) and two quibit 
# (cnot) gates would create the unitary in the first place.

# The initial part of this file is a port of the python code in the Koenig/Smolin paper
# followed by an implementation of the decomposition algorithm in the Aaronson/Gottesman paper.
# There is also code to:
# -  Draw a decomposed circuit 
# -  Display the kets of a |00..0> vector transformed by the circuit 
# -  (Todo) form a 2^n 2^n complex matrix representing the "raw circuit" for use outside this formalism
# -  There is brute force method for finding the minimum set of gates to construct a given unitary 
#    (even with only 3qubits this takes a tediously long time to run)
# -  There is the start of rationalise method, just now it eliminates chains of 4 phase gates or 2 hadamard gates.
#     Maybe the way to go here is to eliminate recognised patterns. [REF] uses a reverse tree building approach
#     which is interesting.
#
# Use is made of the ability of Julia to manipulate its own commands (a sort of hybrid functional/imperative)
# approach. With decomposition I am using two global arrays (so remember to copy them if you want to save them)
# The first "commands" is a text represenetation of the commands needed to rebuild a decomposed end state.
# the second "executeCommands" contains the actual Julia instructions to rebuild the state.
# This is better documented in the appropriate methods.

# The generation of the Clifford requires an element from the symplectic group


using ImageView,Images

# This just shows how big the groups get, Julia overflow is going to be a problem
function getNumberOfCliffords(n)
	return 2^(n^2+2*n)*prod([4^x-1 for x =1:n])
end

# This returns the number without taking into account +/- stuff
function getNumberOfSymplecticCliffords(n)
	return 2^(n^2)*prod([4^x-1 for x =1:n])
end

function getNumberOfBitStringsCliffords(n)
	return 2^(2*n)
end


#takes arrays and places them in a larger array, as follows
#   m1 0
#   0 m2
function directsum(m1,m2)
  n1=size(m1)[1]
  n2=size(m2)[1]
  out = zeros(n1+n2,n1+n2)
  for i = 1:n1
    for j=1:n1
      out[i,j]=m1[i,j]
    end
  end
  for i=1:n2
    for j=1:n2
      out[i+n1,j+n1]=m2[i,j]
     end
   end
  return out
end

# Helper function used in decomposing the clifford
# builds up an array of the commands used to call the hadamard/phase/cnot.
# We are using two globals here for the sake of simplicity.
# execute commands is the Julia code required to call the relevant functions.

function addCommand(comToAdd,comm)
	global commands
	global executeCommands
	commands = vcat(commands,comToAdd)
	executeCommands=append!(executeCommands,[comm]);
end


# Returns the symplectic inner product of two vectors.

function inner(v,w)
	t=0
	for i in 1:(size(v)[1]>>1)
		 t+=v[2*i-1]*w[2*i]
		 t+=w[2*i-1]*v[2*i]
    end
    return t%2
 end


 function transvection(k,v)
 	if (k==0) return mod(v,2)
 	end
 	return mod((v+inner(k,v)*k),2)
 end

 function int2bits(i,n)
 	out = zeros(n)
 	for j in 1:n
 		out[j]=i&1
 		i>>=1
 	end
 	return out
 end

 function findtransvection(x,y)
 	out = zeros(2,size(x)[1])
 	if (x==y)
 		return out
 	end
 	if inner(x,y) ==1
 		out[1,:]=mod((x+y),2)
 		return out
 	end
 	z=zeros(size(x)[1])
 	for i in 1:(size(x)[1]>>1)
 		ii=2*i
 		if ((x[ii-1] + x[ii]) != 0) && ((y[ii-1]+y[ii])!=0)
 			# found the pair
 			z[ii-1]=(x[ii-1]+y[ii-1])%2
 			z[ii]=(x[ii]+y[ii])%2
 			if (z[ii-1]+z[ii])==0 #the same
 				z[ii]=1
 				if (x[ii-1] != x[ii])
 					z[ii-1]=1
 				end
 			end
 			out[1,:]=mod((x+z),2)
 			out[2,:]=mod((y+z),2)
 			return out
 		end
 	end
 	#ok we didn't find a pair
 	# need places where x is 00 and y does not
 	# and vice versa.
 	for i in 1:(size(x)[1]>>1)
 		ii=2*i
 		if ((x[ii-1]+x[ii]) !=0) && ((y[ii-1]+y[ii]) == 0)
 			if x[ii-1]==x[ii]
 				z[ii]=1
 			else
 				z[ii]=x[ii-1]
 				z[ii-1]=x[ii]
 			end
 			break
 		end
 	end
 	for i in 1:(size(x)[1]>>1)
 		ii=2*i
 		if ((x[ii-1]+x[ii]) ==0) && ((y[ii-1]+y[ii]) !=0 )
 			if y[ii-1]==y[ii]
 				z[ii]=1
 			else
 				z[ii]=y[ii-1]
 				z[ii-1]=y[ii]
 			end
 			break
 		end
 	end
 	# (x+z)%2 doesn't work on arrays, so mod function instead.
 	out[1,:]=mod((x+z),2)
 	out[2,:]=mod((y+z),2)
 	return out
 end

# The steps are as set out on page 5 of the paper.
function symplectic(i,n)
	nn=2*n
	s=((1<<nn)-1)
	k=(i%s)+1
	i=int(i/s)
	# step 2
	f1=int2bits(k,nn)
	#step 3
	e1=zeros(nn)
	#define the first basis 
	e1[1]=1
	T=findtransvection(e1,f1)
	# step 4
	if (n > 63) 
		println("Currently we only support up to 63 qubits bits (n) as I have to create a 2*n int and Julia only has 128 bits\n")
		println("The only relevant code is in the function Symplectic, next to this comment. It is easy to change ")
		println("If the need arises");
		return;
	end
	bits=int2bits(i%(convert(Uint128,1)<<(nn-1)),nn-1)

	#step 5
	# Note that in Julia we need to expressly copy
	eprime=copy(e1)
	for j in  3:nn
		eprime[j]=bits[j-1]
	end
	h0=transvection(T[1,:]',eprime)
	h0=transvection(T[2,:]',h0)

	#step 6
	if bits[1]==1
		f1*=0
	end

	#step 7
	id2=[1 0;0 1]
	
	if n!=1
		g=directsum(id2,symplectic(i>>(nn-1),n-1))
	else
		g=id2
	end
	for j in 1:nn
		g[j,:]=transvection(T[1,:]',g[j,:]')
		g[j,:]=transvection(T[2,:]',g[j,:]')
		g[j,:]=transvection(h0,g[j,:]')
		g[j,:]=transvection(f1,g[j,:]')
	end
	return g
end


#Splits up the symplectic into the alpha,beta,gamma and delta arrays that 
#specify the action of the clifford unitary on the X and Z Paulis respectively

function parseSymplectic(symp)
	#note that here we have (so far) ignored that for a,b,c,d we need to multiply an r and s
	#then I am guessing we just use rj as bit string on rhs. ie r1 to rn
	# so with one bit, we have and additional 2*2 * 6
	# then with two bits its 2^2 * 2^2 = 16 
	s2 = int(size(symp)[1]/2)
	a = zeros(s2,s2)
	b = zeros(s2,s2)
	c = zeros(s2,s2)
	d = zeros(s2,s2)
	for i=1:s2
		for j = 1:s2
			a[i,j]=symp[(i-1)*2+1,(j-1)*2+1]
			b[i,j]=symp[(i-1)*2+1,j*2]
			c[i,j]=symp[i*2,(j-1)*2+1]
			d[i,j]=symp[i*2,j*2]
		end
	end
	return (a,b,c,d)
end

# Takes the symplectic, parses it and uses it to create an Aaronson/Gottesman tableau
#The bits specify the 'bit' pattern that controls the sign
#For example in the 1 qubit there are only 6 unique 'symplectic' patterns
#but each of the two rows in the "tableau" can have bit signs of 0 or 1, leading 
# to 6*4 different cliffords.
#this is what the bits is.


function stabiliseSymp(symp,bits)
	(a,b,c,d)=parseSymplectic(symp)

	n=size(a)[1]
	top = hcat(a,b,[ (bits >> (x-1)) %2 for x=1:n])
	bits = div(bits,2)
  	bottom = hcat(c,d,[ (bits >> (x-1)) %2 for x=1:n])
 	state = convert(Array{Int32,2},vcat(top,bottom))
 end


# Takes tableau state and works out what state we are in the the 11 step deomposition
# where we are taking an arbitrary tableau state and "decomposing" it with basic gates
# back to the initial |00...00> state.
function getState(state)
	as = state[:,1:end] 
	# calculate the states, backwards, so we get the "most refined" first
	#define some useful comparison matrices
	# The tableau is split up as follows:
	# A | B
	# ______
	# C | D
	#
	n=div(size(state,1),2) # half the dimension of this 2n x (2n+1) matrix
	state_i = eye(Int32,n)
	state_z = zeros(Int32,n,n)
	state_r = zeros(Int32,n,1)
	if (as==vcat(hcat(state_i,state_z,state_r),hcat(state_z,state_i,state_r)))
		return 11
	end
	#split up into four matrices
	A=as[1:n,1:n]
	B=as[1:n,n+1:2*n]
	C=as[n+1:2*n,1:n]
	D=as[n+1:2*n,n+1:2*n]
	R=as[1:end,end]
	R2=as[n+1:end,end:end]
	if (B==state_z && C==state_z && R2==state_r)
		return 10
	end
	if (C==state_z && A==B  && R2==state_r)
		return 9
	end
	if (A==state_i && C==state_z && D==state_i && R2==state_r)
		return 7
	end
	if (C==state_i && D==state_z && as[n+1:2*n,end:end] == state_r)
		return 6
	end
	if (D==state_z )
		return 5
	end
	if (C==D)
		return 4
	end
	if (C==state_i)
		return 2
	end
	if (rank(C) == n)
		return 1
	end
	return 0
end

# The next series of vectors each implement one of the manipulations specified in Aaronson/Gottesman
# Section VI Canonical form.

function getFullRank(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		r1 = rank(svec[n+1:2*n,1:n]) 
		if r1 == n 
			return
		end
		hadamard(svec,i,false)
		r2 = rank(svec[n+1:2*n,1:n])
		if (r2 == n )
			#println("hadamard(",i,")")
			addCommand("hadamard($i)",Expr(:call,:hadamard,:svec,i))
			return
		end
		if (r2 > r1) 
			#println("hadamard(",i,")",Expr(:call,:hadamard,:svec,i))
			addCommand("hadamard($i)",Expr(:call,:hadamard,:svec,i))
			continue
		else
			hadamard(svec,i,false)
		end
	end
end

function getfirstOne(myVector)
	#println(myVector)
	for i=1:size(myVector,2)
		if myVector[i]==1
			#println("return ",i) 
			return i
		end
	end
	return 0
end


function makeCtheI(svec,offset=-1)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	if offset==-1 offset=n # no offset means we assume its the "C" we are making Int32
	end
	allDone = false;
	while !allDone
		for i=1:n
			allDone = identifytheRow(svec,offset,i)
			if allDone == false
				break
			end
		end
	end
	# now we need double phases to make r_n+1..r_2n to be zero
end

function zapPhase(svec,offset=-1)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	if offset==-1 offset=n # no offset means we assume its the "C" we are making Int32
	end
	for i=1:n
		if (svec[offset+i,end] ==1) 
			phase(svec,i,false)
			addCommand("phase($i)",Expr(:call,:phase,:svec,i))
			phase(svec,i,false)
			addCommand("phase($i)",Expr(:call,:phase,:svec,i))
		end 
	end 
end

# What I mean is take the passed in row, and use cnots to make it an identity element.
function identifytheRow(svec,offset,i)
	#println("Identify row",offset,i)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	#println("N is ",n)
	for j=1:i-1
		#println("Check ",i,",",j)
		if svec[offset+i,j]!= 0
			#println("Wasnt zero")
			if (svec[offset+i,i]==0)
				#println("The diagonal is zero")
					# make the diagonal one first
					# using a hadamard can zap the rank
					# so we should first check if there is a later bit we can cnot into place.
					#t = getfirstOne(svec[offset+i,i+1:n])
					#if (t==0) 
					#	println("OOOPS we couldn't get a one on the diagonl")
						#println("Trying hardy hardimard on bit ",i)
				        #hadamard(svec,i,false)
						#if svec[offset+i,i]==0 #didnt work reverse and not
						#println("that didn't work")
						#if (i==n) #no point in reversing it out if it was the last row
					    #basically we will have flipped a n bit somewhere, so repeat the whole process. 
					#	return false
					#end

					#println(svec)
					#hadamard(svec,i,false)
					#println("Reverse it out and try a cnot from a later one");
					t=getfirstOne(svec[offset+i,i+1:end])
					if (t==0) 
						println("OOOOOOPPPS DEBUG TIME")
						return;
					end
					#println("a: cnot(",t+i,",",i,")")
					cnot(svec,t+i,i,false)
					addCommand("cnot($(t+i),$i)",Expr(:call,:cnot,:svec,j,i))
					#println(svec)
				#else 
					#println("hadmard(svec,",i,")")
				#	addCommand("hadamard($i)",Expr(:call,:hadamard,:svec,i))
				#end
			else
				#println("Diagonal was already one")
			end
			#println("Using the diagonal to zap");
			#cnot from this found bit, to bit i, to turn it to zero
			#println("b: cnot(",i,",",j,")")
			addCommand("cnot($i,$j)",Expr(:call,:cnot,:svec,i,j))
			cnot(svec,i,j,false)
			#println(svec)
		end
	end
	#so up to the "diagonal element" we have them all 0
	#println("Checking Diagonal ",i,",",i)
	if svec[offset+i,i] ==0 # need to get this diagonal equal to 1
		#println("It wasn't one we have offset ",offset, " i ", i)
		t = getfirstOne(svec[offset+i,i+1:end]) + i
		#println("Got a t back of ",t)
		if t== 0
			#println("Couldnt make this the identity\n")
			return false
		end
		#println("cnot(",t,",",i,")")
		cnot(svec,t,i,false)
		addCommand("cnot($t,$i)",Expr(:call,:cnot,:svec,t,i))
		#println(svec)
	end
	# so the diagonal element is now 1, make the rest zero
	for j=i+1:n
		#println("CHecking ",i,",",j)
		if svec[offset+i,j]==1
			#println("It was one")
			#println("cnot(",i,",",j,")")
			cnot(svec,i,j,false)
			addCommand("cnot($i,$j)",Expr(:call,:cnot,:svec,i,j))
			#println(svec)
		end
	end
	return true
end

function diagonaliseD(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		if svec[n+i,n+i]== 0
			phase(svec,i,false)
			#println("phase(",i,")");
			addCommand("phase($i)",Expr(:call,:phase,:svec,i))
		end
	end
end

function diagonaliseB(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		if svec[i,n+i]== 0
			phase(svec,i,false)
			#println("phase(",i,")");
			addCommand("phase($i)",Expr(:call,:phase,:svec,i))
		end
	end
end

function cmdm(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		for j=1:n
			if svec[n+i,j] == svec[n+i,n+j] 
				continue
			end
			if svec[n+i,j] == 0
				cnot(svec,j,i,false)
				#println("cnot(",j,",",i,")")
				addCommand("cnot($j,$i)",Expr(:call,:cnot,:svec,j,i))
			else 
				phase(svec,j,false)
				#println("phase(",j,")")
				addCommand("phase($j)",Expr(:call,:phase,:svec,i))
			end
		end
	end
end
function ambm(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		for j=1:n
			if svec[i,j] == svec[i,n+j] 
				continue
			end
			if svec[i,j] == 0
				cnot(svec,j,i,false)
				#println("cnot(",j,",",i,")")
				addCommand("cnot($j,$i)",Expr(:call,:cnot,:svec,j,i))
			else 
				phase(svec,j,false)
				#println("phase(",j,")")
				addCommand("phase($j)",Expr(:call,:phase,:svec,j))
			end
		end
	end
end

function zapD(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		for j=1:n
			if svec[n+i,n+j]==1
				phase(svec,j,false)
				#println("phase(",j,")")
				addCommand("phase($j)",Expr(:call,:phase,:svec,j))
			end
		end
	end
end

function zapB(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		for j=1:n
			if svec[i,n+j]==1
				phase(svec,j,false)
				#println("phase(",j,")")
				addCommand("phase($j)",Expr(:call,:phase,:svec,j))
			end
		end
	end
end

function hadamardHard(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		hadamard(svec,i,false)
		#println("hadamard(",i,")")
		addCommand("hadamard($i)",Expr(:call,:hadamard,:svec,i))
	end
end

function makeItR0(svec,offset=-1)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	if offset<0
		offset=n
	end
	for i=1:n
		if svec[offset+i,2*n+1]== 1
			phase(svec,i,false)
			phase(svec,i,false)
			#println("phase(",i,")")
			#println("phase(",i,")")
			addCommand("phase($i)",Expr(:call,:phase,:svec,i))
			addCommand("phase($i)",Expr(:call,:phase,:svec,i))
		end
	end
end

# decomposes the "arbitrary" state
# into a series of cnot, hadamard and phase gates.
# the gates are contained in text in the vector of strings, "commands"
# and as Expressions (and thus executable in Julia) in executeCommands

function decomposeState(state,supressOutput = false,rationalise=true)
	global commands
	global executeCommands
	commands=["output(svec)"]
	executeCommands = append!(Expr[],[Expr(:call,:output,:svec)])
    ss1=state
	if (!supressOutput) 
		println("Tableau for unitary: ");
		output(ss1);
	end
	while getState(ss1) < 11
		nextStep(ss1)
	end
	j=div(size(state)[1],2)
  	addCommand("setup($j)",Expr(:call,:setup,j))
	reverse!(commands)
	if rationalise==true
		removeRedundancy(j)
	end

	if (!supressOutput) 
		for i = 1:size(commands,1) 
			println(commands[i]) 
		end
	end

end


# decomposes the "arbitrary" clifford i, with qubits j
# into a series of cnot, hadamard and phase gates.
# the gates are contained in text in the vector of strings, "commands"
# and as Expressions (and thus executable in Julia) in executeCommands

function decompose(i,bits,j=4,supressOutput = false,rationalise=true)
	global commands
	global executeCommands
	commands=["output(svec)"]
	executeCommands = append!(Expr[],[Expr(:call,:output,:svec)])
	ss1=stabiliseSymp(symplectic(i,j),bits)
	if (!supressOutput) 
		println("Tableau for unitary: ");
		output(ss1);
	end
	while getState(ss1) < 11
		nextStep(ss1)
	end
	addCommand("setup($j)",Expr(:call,:setup,i))
	reverse!(commands)
	if rationalise==true
		removeRedundancy(j)
	end

	if (!supressOutput) 
		for i = 1:size(commands,1) 
			println(commands[i]) 
		end
	end

end


function drawCircuit()
	global commands
	currentDir = pwd()
	try
	cOut = open("qasm/temp.qasm","w")
	for i = 1:size(commands,1)
		m = match(r"setup\((.*)\)",commands[i])
		if (m!=nothing)
			for i=1:int(m.captures[1])
    			println(cOut,"    qubit q$i")
			end
    	else
    		m=match(r"hadamard\((.*)\)",commands[i])
    		if (m!=nothing)
    			println(cOut,"     h q",m.captures[1])
    		else 
    			m=match(r"phase\((.*)\)",commands[i])
    			if (m!=nothing) # QASM doesn't appear to have Phase gates, so output an S
    				# I altered QASM to print a P instead of an S.
    				println(cOut,"     S q",m.captures[1])
    			else
    				m=match(r"cnot\((.*),(.*)\)",commands[i])
    				if (m!=nothing)
    					println(cOut,"     cnot q",m.captures[1],",q",m.captures[2])
    				end
    			end
    		end
    	end
    end
    close(cOut)
    cd("qasm")
	test = `csh qasm2png temp.qasm`
	temp = readall(test) 
	finally
	  cd(currentDir)
	end
	img = imread("qasm/temp.png")
end

function nextStep(ss1)
	currentState = getState(ss1)
	if (currentState == 0)
		getFullRank(ss1)
	elseif currentState == 1
		makeCtheI(ss1)
	elseif currentState < 4
		diagonaliseD(ss1)
		cmdm(ss1)
	elseif currentState == 4 
		zapD(ss1)
		makeItR0(ss1)
	elseif currentState == 5
		makeCtheI(ss1)
		zapPhase(ss1)
	elseif currentState == 6
		hadamardHard(ss1)
	elseif currentState < 9
		diagonaliseB(ss1)
		ambm(ss1)
	elseif currentState == 9
		zapB(ss1)
		makeItR0(ss1,0)
	elseif currentState == 10  
		makeCtheI(ss1,0) # 0 = actually its A
		zapPhase(ss1,0)
	elseif currentState == 11
		println("Done")
	end
	return 
end

# maxGates is the number of gates we are allowed (an integer offset to an array)
# gates is the "stack" <- an array of gates we are applying.
# increment index 1, if its greater then maxGate, increment the one above it.
function incrementGate!(gates,maxGate,offset=1)
	if (size(gates,1) < offset ) 
		push!(gates,1)
		return
	end
	currentValue = gates[offset]
	currentValue = currentValue +1
	if (currentValue > maxGate)
		gates[offset]=1
		incrementGate!(gates,maxGate,offset+1)
	else 
		gates[offset]=currentValue
	end
end


function bruteForceBreadthFirst(clifford)
	global svec
	count = 0
	n=div(size(clifford,1),2) # half the dimension of this 2n x (2n+1) matrix
	if (n > 4) 
		println("This might take some time!, there are ",n*(n+1) + n + n, " different gates , we might need ", n*n/log(n), " of them so I *might* have to check ", n*n/log(n)^(n*(n+1) + n + n), "combos !")
		println("Ctrl-C could be your friend")
	end
	gates = (Expr)[]
	svec = []
	# If we have a n quibit system then there are a total of 
	# n*(n+1) possible cnots, n phase gates and n hadamards we can apply
	for i = 1:n
	  for j = 1:n
	  	if i!=j
	  		push!(gates,Expr(:call,:cnot,:svec,i,j,false))
	  	end
	  end
	  push!(gates,Expr(:call,:hadamard,:svec,i,false))
	  push!(gates,Expr(:call,:phase,:svec,i,false))
	end
	# so we now have an array of all possible gates we can apply to the system
	# we want to search through the gates, we "should" need at most n^2 of them.
	# The process will be to have a stack of gates to apply, check them out
	# if we match getState(11) i.e. they reverted the clifford to the initial state
	# we have found the gates needed
	# otherwise we increment our "bottom" gate, and if we have gone though them, recursively increment the 
	# gate above it (incrementGate! does this for us)
	state = 0
	maxGates = size(gates,1);
	gatesToApply=(Int32)[]
	if getState(clifford) == 11
		println("Very funny you don't need to do anything")
		return
	end
	currentAt = 1
	while (state!=11 && size(gatesToApply,1) < n*n*n*n*n*n*n*n*n*n)
		svec = copy(clifford)
		incrementGate!(gatesToApply,maxGates)
		count += 1
		if (size(gatesToApply,1) > currentAt)
			print("Trying: ")
				for i=1:size(gatesToApply,1)
				print(gatesToApply[i]," ")
			end
			println("")
			currentAt+=1
		end
		for index=1:size(gatesToApply,1)
			#println("Going to apply ",gates[gatesToApply[index]])
			eval(gates[gatesToApply[index]])
		end
		state = getState(svec)
	end

	if (state == 11)
		println("Found it after (only) $count permuations")
		println("Gates applied")
		reverse!(gatesToApply)
		for i = 1:size(gatesToApply,1)
			println(gates[gatesToApply[i]])
		end
	else
		println("Didn't find it!")
	end
end

#coding for current bits
# 1 = phase, 2 = hadamard
# 2 + x  = cnot control (with target at x)
# 2 + n + x = cnot target (with control at x)

function checkGate(index,currentBits,n)
	global commands
	global executeCommands
	global toDelete
	checking = commands[index]
	m = match(r"setup\((.*)\)",checking)
	if (m!=nothing)
    	return
    end
    m=match(r"hadamard\((.*)\)",checking)
    if (m!=nothing)
    	bit = int(m.captures[1])
    	if size(currentBits[bit],1) > 0
    		if currentBits[bit][end][1] == 2
    			# two hadamards make a nothing
    			(been,gone) = pop!(currentBits[bit])
    			push!(toDelete,gone)
    			push!(toDelete,index)
			else
				push!(currentBits[bit],(2,index))
			end
		else
				push!(currentBits[bit],(2,index))
		end
		return
	end
	m=match(r"phase\((.*)\)",checking)
    if (m!=nothing) # Its a phase, we need 4 of these
    	bit = int(m.captures[1])
    	if size(currentBits[bit],1) > 2
    		if currentBits[bit][end][1] == 1 && currentBits[bit][end-1][1] == 1 && currentBits[bit][end-2][1] == 1
    			# four phases make a nothing.
    			(b,g) = pop!(currentBits[bit])
    			push!(toDelete,g)
    			(b,g) = pop!(currentBits[bit])
    			push!(toDelete,g)
    			(b,g) = pop!(currentBits[bit])
    			push!(toDelete,g)
    			push!(toDelete,index)
			else
				push!(currentBits[bit],(1,index))
			end
		else
				push!(currentBits[bit],(1,index))
		end
		return
	end
	m=match(r"cnot\((.*),(.*)\)",checking)
    if (m!=nothing)
    		cbit = int(m.captures[1])
    		tbit = int(m.captures[2])
    		# so just now we are going to catch two cnots in a row
    		# needs to match both bits.
    		cbitNo = 2+tbit
    		tbitNo = 2+n+cbit
    		if (size(currentBits[cbit],1) > 0 && size(currentBits[tbit],1) > 0 && currentBits[cbit][end][1] == cbitNo &&
    			currentBits[tbit][end][1] == tbitNo)
    			(b,g)=pop!(currentBits[cbit])
    			push!(toDelete,g)
    			(b,g)=pop!(currentBits[tbit])
    			# dont need to push the deletion of this gate twice.
    			push!(toDelete,index)
    		else
    			push!(currentBits[cbit],(cbitNo,index))
    			push!(currentBits[tbit],(tbitNo,index))
    		end
    end
end


function removeRedundancy(n)
	global commands
	global executeCommands
	global toDelete
	toDelete = (Int32)[]
	bitsOf = (Array{(Int32,Int32)})[]
	for i=1:n
    	push!(bitsOf,[])
	end
	for i=1:size(commands,1)
		checkGate(i,bitsOf,n)
	end
	deleteat!(commands,sort!(toDelete))
	deleteat!(executeCommands,toDelete)
end







