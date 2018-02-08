
function makeSuper(u)
    return round(real(liou2pauliliou(liou(u))),10)
end

function findInThis(needle,haystack)
    tofind = round(needle,10)
    return findfirst(x->round(x,10)==tofind,haystack)
end


#Super operator version
function sfindInThis(needle,haystack)
    tofind = makeSuper(needle)
    return findfirst(x->round(x,10)==tofind,haystack)

end

function checkFrame(x)
    sum = 0
    for i=1:length(x)
        for j=1:length(x)
            sum += abs(trace(x[i]'*x[j]))^4
        end
    end
    sum/(length(x)^2)
end

# Just some convenience functions
operatorPaulis= Array{Complex{Float64},2}[[1 0;0 1], [0 1;1 0],[0 -im;im 0], [1 0;0 -1]];
superPaulis = Array{Float64,2}[makeSuper(i) for i in operatorPaulis];
pI=operatorPaulis[1]
pX=operatorPaulis[2]
pY=operatorPaulis[3]
pZ=operatorPaulis[4];

piBy8 = exp(im*π/8)*[exp(-im*π/8) 0;0 exp(im*π/8)]
superPiBy8 = makeSuper(piBy8)


superCliffords = Array{Float64,2}[round(generateClifford(i,0,0,0),5) for i=1:24];


function showUseful() 

print("Defined: makeSuper\n")
print("Defined: superCliffords\n")
print("Defined: findInThis(needle,haystack)\n")
print("Defined: sfindInThis(needle,haystack) (throws into makeSuper first, haystack in super basis)\n")
print("Defined: checkFrame(Array) -- assumes non-super form\n")
print("Defined: operatorPaulis\n")
print("Defined: superPaulis\n")
print("Defined: pI,pX,pY,pZ\n")
print("Defined: piBy8 and superPiBy8\n")

end