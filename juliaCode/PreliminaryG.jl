
"""
Just implements a Z rotation in CHP gates
"""
function Z(s,b)
    phase(s,b,false)
    phase(s,b,false)
    return s
end

"""
Implement an X rotation in CHP gates
"""
function X(s,b)
    hadamard(s,b,false)
    phase(s,b,false)
    phase(s,b,false)
    hadamard(s,b,false)
    return s
end

"""
Implements a two logical qubit block
"""
type G4 <: Integer
    state:: Array{Int32,2}
    qubits:: Array{Int32,1}
    function G4()
        new(gSetup(),[1,2,3,4]);
    end
    function G4(x::G4)
        # Note this copies the state, the state might be bigger than 4 qubits
        state = x.state
        qubits = x.qubits
        new(state,qubits)
    end
end

function allStates(q::G4)
    print(getStatesFor(q.state))
end

function gSetup()
    s = setup(4);
    hadamard(s,2,false)
    cnot(s,2,3,false)
    cnot(s,2,1,false)
    cnot(s,3,4,false)
    return s
end


function init(q::G4)
    q.state = gSetup()
end


function logicalZ(q::G4,b)
    if b==1
        q.state = Z(q.state,1)
        q.state = Z(q.state,2)
    else
        q.state=Z(q.state,1)
        q.state=Z(q.state,3)
    end
end

function logicalX(q::G4,a)
    if a == 1
        q.state = X(q.state,1)
        q.state = X(q.state,3)
    else
        q.state = X(q.state,1)
        q.state = X(q.state,2)
    end
end

function logicalHHSwap(q::G4)
    hadamard(q.state,1,false)
    hadamard(q.state,2,false)
    hadamard(q.state,3,false)
    hadamard(q.state,4,false)
end

function logicalCPHASEZZ(q::G4)
    phase(q.state,1,false)
    phase(q.state,2,false)
    phase(q.state,3,false)
    phase(q.state,4,false)
end



function logicalZ(s,a)
    if a == 1
        s = Z(s,1)
        s = Z(s,2)
    else
        s = Z(s,1)
        s = Z(s,3)
    end
    return s
end


function logicalX(s,a)
    if a == 1
        s = X(s,1)
        s = X(s,3)
    else
        s = X(s,1)
        s = X(s,2)
    end
    return s
end

function logicalHH(s)
    hadamard(s,1,false)
    hadamard(s,2,false)
    hadamard(s,3,false)
    hadamard(s,4,false)
    return s
end
