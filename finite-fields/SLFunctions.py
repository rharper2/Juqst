#SLFunctions.py

from finitefield import *
import numpy as np
import numpy.linalg as lg 

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
    return w,abs(lg.inv(w))

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

def getAlpha(field,no):
    """ Returns the actualy block [alpha,beta//gamma,delta]
        nicely wrapped up in a numpy array. Which allows us
        to nicely np.dot() [array multiply] it later
        Because of the way types work it all just happens in the correct
        field.
    """
    return np.array([ [field[no[0]],field[no[1]] ], [ field[no[2]],field[no[3]]]])


def createStabiliserBasis(x):
    """ Returns the matrices for the field x
        That are required to form the stabiliser states
        You would normally multiply this by the relevant SL matrix
     to get the stabiliser state that represents the clifford you want
    
    """    
    matrices=[]
    w,wb=getWForPoly(x)
    for i in range(0,x.degree):
        makingIt=[]
        for pre in range(0,i):
            makingIt = makingIt+[0]
        makingIt=makingIt+[1]
        if (i+1<= x.degree):
            for aft in range(i+1,x.degree):
                makingIt=makingIt+[0]
        matrices = matrices + [makingIt]
    samples=[]
    for i in matrices:
        samples = samples + [np.array([[x(i)],[x([0])]])]
    for i in matrices:
        samples = samples + [np.array([[x(0)],[toPoly(x,translateWithW(wb,x(i)))]])]
    return samples

def polyToList(x):
    """
        Takes a field and returns a list of coefficients, including redundant zeros
    """
    deg = x.degree
    listx = x.poly.coefficients
    toRet=[]
    for i in listx:
        if i==1:
            toRet.append(1)
        else:
            toRet.append(0)
    for i in range(len(listx),deg):
        toRet.append(0)
    return toRet

def stabiliserFori(fieldType,field,slm,no,stabsA):
    fullStabiliser = []
    w,wb=getWForPoly(fieldType)
    mat=getAlpha(field,slm[no])
    for i in stabsA:
          temp =  np.dot(mat,i)
          fullStabiliser.append(polyToList(temp[0,0])+polyToList(toPoly(fieldType,translateWithW(w,temp[1,0]))))
     #     fullStabiliserZ.append(polyToList(toPoly(fieldType,translateWithW(w,temp[1,0]))))
    return fullStabiliser

def getAllStabilisers(fieldType):
    field,slm = setUpSLGroups(fieldType)
    stabsA = createStabiliserBasis(fieldType)
    stabiliserList = []
    for i,idx in enumerate(slm):
        stabiliserList.append(stabiliserFori(fieldType,field,slm,i,stabsA))
    return stabiliserList

            
def getF2SLGroups():
    field,slm = setUpSLGroups(FiniteField(2,2))
    return (field,slm)

def getF3SLGroups():
    field,slm = setUpSLGroups(FiniteField(2,3))
    return (field,slm)


def getF2Stabilisers():
    return getAllStabilisers(FiniteField(2,2))

def getF3Stabilisers():
    return getAllStabilisers(FiniteField(2,3))