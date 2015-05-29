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

function addCommand(comToAdd,comm)
	global commands
	global executeCommands
	commands = vcat(commands,comToAdd)
	executeCommands=append!(executeCommands,[comm]);
end



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
	bits=int2bits(i%(convert(Uint128,1)<<(nn-1)),nn-1)

	#step 5
	# deep copy? Check if needed
	eprime=deepcopy(e1)
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

function parseSymplectic(symp)
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

function stabiliseSymp(symp)
	(a,b,c,d)=parseSymplectic(symp)
	n=size(a)[1]
	top = hcat(a,b,zeros(Int32,n,1))
  	bottom = hcat(c,d,zeros(Int32,n,1))
 	state = convert(Array{Int32,2},vcat(top,bottom))
 end

function getState(state)
	as = state[:,1:(size(state,2)-1)] # ignore the sign bits for just now
	# calculate the states, backwords, so we get the "most refined" first
	#define some useful comparison matrices
	n=div(size(state,1),2) # half the dimension of this 2n x (2n+1) matrix
	state_i = eye(Int32,n)
	state_z = zeros(Int32,n,n)
	if (as==vcat(hcat(state_i,state_z),hcat(state_z,state_i)))
		return 11
	end
	#split up into four matrices
	A=as[1:n,1:n]
	B=as[1:n,n+1:2*n]
	C=as[n+1:2*n,1:n]
	D=as[n+1:2*n,n+1:2*n]
	if (B==state_z && C==state_z)
		return 10
	end
	if (C==state_z && A==B)
		return 9
	end
	if (A==state_i && C==state_z && D==state_i)
		return 7
	end
	if (C==state_i && D==state_z)
		return 6
	end
	if (D==state_z)
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
			println("hadamard(svec,",i,")")
			addCommand("hadamard(svec,$i)",Expr(:call,:hadamard,:svec,i))
			return
		end
		if (r2 > r1) 
			println("hadamard(svec,",i,")",Expr(:call,:hadamard,:svec,i))
			addCommand("hadamard(svec,$i)",Expr(:call,:hadamard,:svec,i))
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
end

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
					# some difficulty with losing rank, so the first 
					# thing to try is a hadmard on this quibit
				#println("Trying hardy hardimard on bit ",i)
				hadamard(svec,i,false)
				if svec[offset+i,i]==0 #didnt work reverse and not
					#println("that didn't work")
					if (i==n) #no point in reversing it out if it was the last row
					#basically we will have flipped a n bit somewhere, so repeat the whole process. 
					#	println("Going to return false")
						return false
					end
					#println(svec)
					hadamard(svec,i,false)
					#println("Reverse it out and try a cnot from a later one");
					t=getfirstOne(svec[offset+i,i+1:end])
					if (t==0) 
						println("OOOOOOPPPS DEBUG TIME")
						return;
					end
					cnot(svec,t+i,i,false)
					println("cnot(svec,",j,",",i,")")
					addCommand("cnot(svec,$j,$i,)",Expr(:call,:cnot,:svec,j,i))
				else 
					println("hadmard(svec,",i,")")
					addCommand("hadamard(svec,$i)",Expr(:call,:hadamard,:svec,i))
				end
			else
				#println("Diagonal was alread one")
			end
			#println("Using the diagonal to zap");
				#cnot from this found bit, to bit i, to turn it to zero
			println("cnot(svec,",i,",",j,")")
			addCommand("cnot(svec,$i,$j)",Expr(:call,:cnot,:svec,i,j))
			cnot(svec,i,j,false)
			#println(svec)
		end
	end
	#so up to the "diagonal element" we have them all 0
	#println("Checking Diagonal ",i,",",i)
	if svec[offset+i,i] ==0 # need to get this diagonal equal to 1
		#println("It wasn't one")
		t = getfirstOne(svec[offset+i,i+1:end]) + i
		if t== 0
			#println("Couldnt make this the identity\n")
			return false
		end
		cnot(svec,t,i,false)
		println("cnot(svec,",t,",",i,")")
		addCommand("cnot(svec,$t,$i)",Expr(:call,:cnot,:svec,t,i))
		#println(svec)
	end
	# so the diagonal element is now 1, make the rest zero
	for j=i+1:n
		#println("CHecking ",i,",",j)
		if svec[offset+i,j]==1
			#println("It was one")
			cnot(svec,i,j,false)
			println("cnot(svec,",i,",",j,")")
			addCommand("cnot(svec,$i,$j)",Expr(:call,:cnot,:svec,i,j))
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
			println("phase(svec,",i,")");
			addCommand("phase(svec,$i)",Expr(:call,:phase,:svec,i))
		end
	end
end

function diagonaliseB(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		if svec[i,n+i]== 0
			phase(svec,i,false)
			println("phase(svec,",i,")");
			addCommand("phase(svec,$i)",Expr(:call,:phase,:svec,i))
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
				println("cnot(svec,",j,",",i,")")
				addCommand("cnot(svec,$j,$i)",Expr(:call,:cnot,:svec,j,i))
			else 
				phase(svec,j,false)
				println("phase(svec,",j,")")
				addCommand("phase(svec,$j)",Expr(:call,:phase,:svec,i))
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
				println("cnot(svec,",j,",",i,")")
				addCommand("cnot(svec,$j,$i)",Expr(:call,:cnot,:svec,j,i))
			else 
				phase(svec,j,false)
				println("phase(svec,",j,")")
				addCommand("phase(svec,$j)",Expr(:call,:phase,:svec,j))
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
				println("phase(svec,",j,")")
				addCommand("phase(svec,$j)",Expr(:call,:phase,:svec,j))
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
				println("phase(svec,",j,")")
				addCommand("phase(svec,$j)",Expr(:call,:phase,:svec,j))
			end
		end
	end
end

function hadamardHard(svec)
	n=div(size(svec,1),2) # half the dimension of this 2n x (2n+1) matrix
	for i=1:n
		hadamard(svec,i,false)
		println("hadamard(svec,",i,")")
		addCommand("hadamard(svec,$i)",Expr(:call,:hadamard,:svec,i))
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
			println("phase(svec,",i,")")
			println("phase(svec,",i,")")
			addCommand("phase(svec,$i)",Expr(:call,:phase,:svec,i))
			addCommand("phase(svec,$i)",Expr(:call,:phase,:svec,i))
		end
	end
end


function temp()

	#println(ss1)
	println("Initial state",getState(ss1))
	getFullRank(ss1)
	#println(ss1)
	println("State now",getState(ss1))
	makeCtheI(ss1)
	#println(ss1)
	println("State now",getState(ss1))
	diagonaliseD(ss1)
	#println(ss1)
	println("State now",getState(ss1))
	cmdm(ss1)
	println("State now",getState(ss1))
	zapD(ss1)
	println("State now",getState(ss1))
	makeItR0(ss1)
	println("State now",getState(ss1))
	makeCtheI(ss1)
	println("State now",getState(ss1))
	hadamardHard(ss1)
	println("State now",getState(ss1))
	diagonaliseB(ss1)
	println("State now",getState(ss1))
	ambm(ss1)
	println("State now",getState(ss1))
	zapB(ss1)
	println("State now",getState(ss1))
	makeItR0(ss1,0)
	makeCtheI(ss1,0) # 0 = actually its A
	println("State now",getState(ss1))
	end

function sofar(i,j=4)
	global commands
	global executeCommands
	commands=["output(svec)"]
	executeCommands = append!(Expr[],[Expr(:call,:output,:svec)])
	ss1=stabiliseSymp(symplectic(i,j))
	while getState(ss1) < 11
		nextStep(ss1)
	end
	addCommand("setup($j)",Expr(:call,:setup,i))
	reverse!(commands)
	for i = 1:size(commands,1) println(commands[i]) end
	return ss1
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
	elseif currentState == 11
		println("Done")
	end
	println("State now",getState(ss1))
	return 
end


	





