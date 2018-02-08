# Generate random channels - 

using Distributions
# Note this is for reproducibility --- for multiple runs or spreadsheets you will need some 
# other way to set.
srand(4321)

d=Normal()

"""
    Simple function that takes a normal distribution and returns a random variable
    that is bounded (non-inclusively) by min and max.
"""
function genTruncated(min,max,average,sigma,normalDistribution)
    x=0.0
    while (x < min) || (x > max)
        x=rand(normalDistribution,1)[1]*sigma+average
    end
    return x
end
    
"""
    Given the eigenvalues λ return a
    vector of random τ, constrained by the λ
""" 
function getτ(λ1,λ2,λ3)
    Maxτ1 = 1-abs(λ1)
    Maxτ2 = 1-abs(λ2)
    Maxτ3 = 1-abs(λ3)
    τ1 = (Maxτ1)*(2*rand()-1)
    τ2 = (Maxτ2)*(2*rand()-1)
    τ3 = (Maxτ3)*(2*rand()-1)
    return [τ1 τ2 τ3]
end

"""
    q(l) and Zη(t,l) See equations 20 and 21 of 1707.06926 (Rudnicki et al)
"""
function q(l)
 return (1+l[1]+l[2]+l[3])*(1 + l[1]-l[2]-l[3])*(1-l[1]+l[2]-l[3])*(1-l[1]-l[2]+l[3])
end

"""
    q(l) and Zη(t,l) See equations 20 and 21 of 1707.06926 (Rudnicki et al)
"""
function Zη(t,l)
    norm(t)^4-2*norm(t)^2-2*sum([l[i]^2*(2*(t[i]^2)-norm(t)^2) for i=1:3]) + q(l)
end


"""
    Main function to generate the 'internal' matrix of the channel i.e. the one 
    in the form where the unital part consists of diagonal numbers, and the we have
    a non-unital part. This can then be acted on by unitaries either side (rotations) 
    to generate an arbitarty CPTP channel
"""
function genSigma(min,max,average,sigma)
    λ1=genTruncated(min,max,average,sigma,d)
    λ2=genTruncated(min,max,average,sigma,d)
    λ3=genTruncated(min,max,average,sigma,d)
    while 1+λ3 < abs(λ1+λ2) || 1-λ3 < abs(λ2-λ1) 
        λ3=genTruncated(min,max,average,sigma,d)
    end
    τ = getτ(λ1,λ2,λ3)
    while norm(τ)^2 > 1-sum([x^2 for x in [λ1 λ2 λ3]])+2*λ1*λ2*λ3 || Zη(τ,[λ1 λ2 λ3]) < 0
        τ = getτ(λ1,λ2,λ3)
    end
    return [1 0 0 0;τ[1] λ1 0 0;τ[2] 0 λ2 0;τ[3] 0 0 λ3]
end
    
""" 
    Rotate the sigma channel randomly to generate an arbitrary channel
"""
function genChannel(r1,sm,sv,r2)
    createRotation(rand(d,3)*r1) * genSigma(0.9,1,sm,sv) * createRotation(rand(d,3)*r2)'
end

"""
    given a vector of three angles, generate a random rotation channel
    by using the vector to control a Z*X*Z rotation
"""
function createRotation(angles)
    rotation = eye(4)
    rotation[2:4,2:4]=
    [cos(angles[1]) sin(angles[1]) 0;-sin(angles[1]) cos(angles[1]) 0;0 0 1] *
    [1 0 0;0 cos(angles[2]) sin(angles[2]);0 -sin(angles[2]) cos(angles[2])] *
    [cos(angles[3]) sin(angles[3]) 0;-sin(angles[3]) cos(angles[3]) 0;0 0 1] 
    return rotation
end

"""
    returns the fidelity of the 1 qubit channel (in PauliLiouville basis)
"""
function fidelity(channel)
    (1+1/3*(channel[2,2]+channel[3,3]+channel[4,4]))/2
end
    
"""
    Some helper funcation to generate typical, high fidelity channels and lower fidelity
    state preperation and measurement noise channels.
"""

function randomFidelityNoise()
    return genChannel(0.06,0.998,0.04,0.06)
end

function randomPrepNoise()
    return genChannel(0.05,0.985,0.15,0.05)
end

function randomMeasureNoise()
    return genChannel(0.05,0.98,0.15,0.05)
end


""" 
    Simple helper function - designed to help generate a map of a specific fidelity
    In general though you would want to rewrite this to give genChannel
    parameters that will result in random channels that often bracket the desired fidelity
"""
function genChannelMap(f)
    while true
        candidate = randomFidelityNoise()
        if abs(f-fidelity(candidate)) < 0.0001
            return candidate
        end
    end
end