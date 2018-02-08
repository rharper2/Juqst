using Cliffords
#import Base.writemime
import Base.show
# Changed slightly from the version on BitBucket
# I added a translate in the paulliou conversion functions, so that the basis was IXYZ (rather than IXZY)
# I removed a division by d when converting from liou 2 choi (and vice versa) (slightly different conventions)
# The choi2liou uses a different involution - keeping as column stacked.
# (esp for 1 quibit - the X0,0 should be 1, not 0.5)

export mat,
       liou,
       choi_liou_involution,
       swap_involution,
       choi2liou,
       choiX2liou,
       choi2kraus,
       kraus2choi,
       kraus2liou,
       liou2choi,
       liou2choiX,
       liou2kraus,
       liou2pauliliou,
       pauliliou2liou,
       depol,
       istp,
       iscp,
       ischannel,
       isunital,
       nearestu,
       nicePrint,
       writemime

mat( v::Vector, r=round(Int,sqrt(length(v))), c=round(Int,sqrt(length(v))) ) = reshape( v, r, c )
liou{T<:AbstractMatrix}( left::T, right::T ) = kron( transpose(right), left )
liou{T<:AbstractMatrix}( m::T ) = kron( conj(m), m )

function nicePrint(m)
    start = string("\\begin{equation*} \\left(\\begin{array}{*{11}c}")
    for i = 1:size(m,1)
        for j = 1:size(m,2)
            start = "$start $(m[i,j]) &"
        end
        start= "$start \\\\"
    end
    start = "$start \\end{array}\\right)\\end{equation*}"
    return start
end

show(stream, ::MIME"text/latex", x::AbstractMatrix) =write(stream, nicePrint(x)) 

#writemime(io, ::MIME"text/latex", x::Matrix) = write(io, nicePrint(x)) 

# NOTE USING the Rc involution. 1,3,2,4 or 4,2,1,3?
function choi_liou_involution( r::Matrix )
  d = round(Int, sqrt(size(r,1)) )
  rl = reshape( r, (d, d, d, d) )
  rl = permutedims( rl, [1,3,2,4] )
  reshape( rl, size(r) )
end


# Changed again this was 3, 4, 1 2 - changed to 2, 1, 4, 3 i.e. L(x \otimes y)->L(y \otimes x)
function swap_involution( r::Matrix )
  d = round(Int, sqrt(size(r,1)) )
  rl = reshape( r, (d, d, d, d) )
  rl = permutedims( rl, [2,1,4,3] )
  reshape( rl, size(r) )
end

function choi2liou( r::Matrix  )
  choi_liou_involution( r )
end

function liou2choi( r::Matrix )
  choi_liou_involution( r )
end

function liou2choiX(r::Matrix)
  d = round(Int, sqrt(size(r,1)) )
  rl = reshape( r, (d, d, d, d) )
  rl = permutedims( rl, [3,1,4,2] )
  reshape( rl, size(r) )'
end

function choiX2liou(r::Matrix)
  d = round(Int, sqrt(size(r,1)) )
  rl = reshape( r, (d, d, d, d) )
  rl = permutedims( rl, [3,1,4,2] )
  reshape( rl, size(r) )
end

# This stays the same, save that we need to take into account dimensional factors. 
function choi2kraus{T}( r::Matrix{T}  )
  d = sqrt(size(r,1))
  (vals,vecs) = eig( d*r )
  #vals = eigvals( sqrt(size(r,1))*r )
  kraus_ops = Matrix{T}[]
  for i in 1:length(vals)
    push!(kraus_ops, sqrt(vals[i])*mat(vecs[:,i]))
  end
  factor = trace(sum([k'*k for k in kraus_ops]))
  kraus_ops = kraus_ops/sqrt(factor/d)
  kraus_ops
end



function kraus2liou{T}( k::Vector{Matrix{T}} )
  l = zeros(T,map(x->x^2,size(k[1])))
  for i in 1:length(k)
    l = l + liou(k[i],k[i]')
  end
  l
end

function liou2kraus( l::Matrix )
  choi2kraus( liou2choi( l ) )
end

function kraus2choi{T}( k::Vector{Matrix{T}} )
  liou2choi(kraus2liou(k))
  # Seems to fail for some kraus vectors.
  #c = zeros(T,map(x->x^2,size(k[1])))
  #for i in 1:length(k)
  #  c = vec(k[i])*vec(k[i])'
  #end
  #c/sqrt(size(c,1))
end


_translate(v) = [x==2?3:x==3?2:x for x in v] 

_num2quat(n,l) = _translate(map(s->parse(Int,s),collect(base(4,n,l))))

function pauliliou2liou( m::Matrix )
  if size(m,1) != size(m,2)
    error("Only square matrices supported")
  elseif size(m,1) != 4^(floor(log2(size(m,1))/2))
    error("Only matrices with dimension 4^n supported.")
  end
  dsq = size(m,1)
  res = zeros(Complex128,size(m))
  l = round(Int,log2(dsq)/2)
  for i=1:dsq
    for j=1:dsq
      res += m[i,j] * vec(complex(Pauli(_num2quat(i-1,l)))) * vec(complex(Pauli(_num2quat(j-1,l))))' / sqrt(dsq)
    end
  end
  res
end

function liou2pauliliou{T}( m::Matrix{T} )
  if size(m,1) != size(m,2)
    error("Only square matrices supported")
  elseif size(m,1) != 4^(floor(log2(size(m,1))/2))
    error("Only matrices with dimension 4^n supported.")
  end
  dsq = size(m,1)
  res = zeros(Complex128,size(m))
  l = round(Int,log2(dsq)/2)
  for i=1:dsq
    for j=1:dsq
      res[i,j] += trace( m * vec(complex(Pauli(_num2quat(j-1,l)))) * vec(complex(Pauli(_num2quat(i-1,l))))' / sqrt(dsq) )
    end
  end
  res
end

function choi2chi{T}( m::Matrix{T} )
  if size(m,1) != size(m,2)
    error("Only square matrices supported")
  elseif size(m,1) != 4^(floor(log2(size(m,1))/2))
    error("Only matrices with dimension 4^n supported.")
  end
  dsq = size(m,1)
  dim = round(Int(log2(dsq)/2))
  thepaulis= allpaulis(dim)
  #using the transpose .' to conform with python row major.
  stat=vec(complex(thepaulis[1]).')
  for i= 2:4^dim
    stat = hcat(stat,vec(complex(thepaulis[i]).'))
  end
  stat'*m*stat

end

function pauliliou2chi{T}( m::Matrix{T})
  choi2chi(liou2choi(pauliliou2liou(m)))
end

"""
Returns a superoperator that replaces the input with a maximally
mixed state with probability p, and leaves it unchanged with probability (1-p).
"""
function depol( d::Int, p=1.0 )
  choi2liou( p * eye(d^2)/d^2 + (1-p) * projector(_max_entangled_state(d)) )
end

"""
Given a superoperator, it extracts the closest superoperator (in Frobenius norm)
that is unital. The result may not be completely positive.
"""
function unitalproj{T}( m::Matrix{T} )
  d2 = size(m,1)
  d  = round(Int,sqrt(d2))
  id = projector(normalize(vec(eye(d))))
  id*m*id + (I-id)*m*(I-id)
end

function iscp(m; tol=0.0)
    evs = eigvals(liou2choi(m))
    tol = tol==0.0 ? eps(abs.(one(eltype(m)))) : tol
    all(real(evs) .> -tol) && all(abs.(imag(evs)) .< tol)
end

function istp(m; tol=0.0)
    tol = tol==0.0 ? eps(abs.(one(eltype(m)))) : tol
    dsq = size(m,1)
    d = round(Int,sqrt(dsq))
    norm(m'*vec(eye(d))-vec(eye(d)),Inf) < tol
end

function ischannel(m; tol=0.0)
    #println(iscp(m,tol=tol))
    #println(istp(m,tol=tol))
    iscp(m,tol=tol) && istp(m,tol=tol)
end

function isunital(m; tol=0.0)
    tol = tol==0.0 ? eps(abs(one(eltype(m)))) : tol
    dsq = size(m,1)
    d = round(Int,sqrt(dsq))
    norm(m*vec(eye(d))-vec(eye(d)),Inf) < tol
end

"""
Computes the unitary CP map closest (interferometrically) to a given CP map. 
See D. Oi, [Phys. Rev. Lett. 91, 067902 (2003)](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.91.067902)
"""
function nearestu(l)
    c = liou2choi(l)
    vals,vecs = eig(Hermitian(c))
    imax = indmax(vals)
    Λ = mat(vecs[:,imax])
    U,Σ,V = svd(Λ)
    W = U*V'
    return kron(conj(W),W)
end
