#SLFunctions.py

from finitefield import *
import numpy as np

# just enumerate all possible field elements (size 2^n)
def getExpansion(nos):
    if nos==1:
        return [[0],[1]];
    else:
        list = []
        sofar=getExpansion(nos-1)
        for i in sofar:
            add1 = i + [0]
            add2 = i + [1]
            list = list + [add1] + [add2]
        return list
 
# This will equal an array of 2^(3n)-2^n
def setUpSLGroups(x):
  # Just now assuming GF(2^n)
  count=0;
  elements=getExpansion(x.degree)
  fieldElements=[]
  for i in elements:
    fieldElements.append(x(i))
  SLMembers=[]
  for i,ie in enumerate(fieldElements):
    for j,je in enumerate(fieldElements):
        for k,ke in enumerate(fieldElements):
            for l,le in enumerate(fieldElements):
                if ie*le+je*ke == x([1]):
                    SLMembers.append([i,j,k,l])
  return fieldElements,SLMembers

def getPolyBasis(x):
    a=[1]
    basis=[]
    for i in range(0,x.degree):
        basis = basis + [x(a)]
        a = [0]+a
    return basis

def getWForPoly(x):
    basis = getPolyBasis(x)
    xlen = len(basis)
    w = np.zeros((xlen,xlen))
    for i,ip in enumerate(basis):
        for j,jp in enumerate(basis):
            if ((ip*jp).trace() == x([1])):
              w[i,j]=1
            else:
              w[i,j]=0
    return w

def getDualBasisForP(x):
    dualBasis=[]
    primal = getPolyBasis(x)
    w=getWForPoly(x,primal)
    for pr in primal:
        dualBasis = dualBasis + [x(translateWithW(w,pr))]
    return  dualBasis

def translateWithW(W,anX):
    mod = anX.modulus
    modP = IntegersModP(mod)
    listx=anX.poly.coefficients
    if (len(W) != anX.degree):
        raise Exception('The supplied W and the degreee of the field must agree')
    xlen=len(W)
    retC=[]
    for i in range(0,xlen):
        sum = modP(0)
        for j,lx in enumerate(listx):
            sum+=W[i,j]*lx
        retC =  retC + [sum]
    return (retC)
  
#Not quite sure why I needed this
#Takes the field (.e.g F23) and the modded list from translate and
# returns the polynomial        
def toPoly(x,tl):
    a=[]
    my1=x([1])
    for t in tl:
        if t == my1:
            a=a+[1]
        else:
            a=a+[0]
    return x(a)  


            