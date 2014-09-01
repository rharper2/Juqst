# Juqst JUlia Quantum Simulator Toolbox
##Beginning of a quantum simulator toolbox, primarily written in the Julia Language

###Done:

* Implement the simulation of Stabiliser Circuits (Aaronson/Gottesman arXiv:quant-ph/0406196)
* Implement ability to select an arbitrary Clifford group element (Koenig/Smolin arXiv:quant-ph/1406.2170) (not quite complete - phases not implemented properly yet)
* Implement the ability to decompose an arbitrary clifford unitary into a quantum circuit consistiting of hadamard, phase and two-qubit cnot gates.(Aaronson/Gottesman arXiv:quant-ph/0406196)

###Working on:

- Draw the quantum circuit resulting from the decomposition of the clifford/unitary
- Shadow the stabiliser state with the exact density matrix representing the state (the **"base representation"**)
- Output the "ket" state represented by the  (Aaronson/Gottesman arXiv:quant-ph/0406196) tabelau (the **"tableau"**)
- Integrate earlier work pre-defining the steane code generators and logical operators, allowing an arbitrary qubit to be projected into a steane code stabilised state.

###To do:
- randomised benchmarking from arbitrary cliffords
- introduce the ability to add noise at the level of the "base representation" as well as in the stabilised state
- explore noisy quantum channels
- integrate tomography work

# To install

This has been developed on Julia 3.0. Currently there are two files that need to be loaded.

Move to the directory containing these files 
    cd("juqst")
Then 
    require("Initial.jl")
and 
    require("Symplectic.jl")


# Sample use

## Stabiliser Circuits

    state = setup(number_ofQubits)

prepares the stabiliser state for the correct number of qubits in the |000..000> basis state

The state is represented internally as a matrix of the form:

<TO BE INSERTED>

Currently I am just using Int32 Arrays, although binary arrays would save space (if it ever becomes necessary).

    output(state)

Prints the state in a human readable form. The states above the line are the 'destabiliser' state, below the line are the 'stabiliser' states. 

So in a 3 qubit system the initial state of |000> is coded as 

```
XII
IXI
IIX
---
ZII
IZI
IIZ
```

The following commands are defined

    hadamard(state,qubit)  # apply a hadamard to the relevant qubit
    phase(state,qubit)     # apply a phase gate to the relevant qubit
    cnot(state,control,target) # apply a controlled not from control qubit to target qubit

Output of the resultant state can be supressed by adding an extra false parameter
    hadamard(state,qubit,false) # hadamard as before, but supress output

**NOTE! that these commands alter the state passed into them. I have broken Julia convention which requires functions 
with side effects to be written thus - hadamard!(state,qubit).**

## Arbitrary cliffords

(Koenig/Smolin arXiv:quant-ph/1406.2170)

The basis behind this papers is that we can implement a one-to-one mapping between the cliffords and an integer (plus a random phase string).

The mapping is as follows:

<img src="Clifford Mapping.png">Mapping png</img>

We can generate the alpha,beta,gamma and delta via

`symplectic(i,n) # i = integer represting the clifford, n is the number of qubits

Which returns the nxn arrays (alpha->delta) coded as follows:

<img src="coding.png">coding png</img>

More usefully these can be placed into a stabiliser tableau (that is the equivlent of passing the state |0000> through a gate that implements the unitary in question as follows:

    stabiliseSymp(symp) # where symp is the symplectic of the clifford generated.

e.g.

    state = stabiliseSymp(symplectic(23,4)) # for the tableau of clifford '23' in a 4 qubit system

Of course there are actually 4^4 version of of symplectic 23, because of the different phases that we can have - this will be implemented shortly.

# Decomposing a tableau (such as clifford)

This will be made more general, but just now it decomposes an arbitrary clifford

    decompose(clifford_number, qubits)

This prints out the elementary gates that would reconstruct the relevant clifford unitary.

The commands are stored as string in the vector commands
The commands are also stored as Julia code in the vector executeCommands (so you can for instance execute them to re-create the tableau)

Shortly I will add instructions as to how to draw the circuit represented by these gates.

    getState(state) 

Is a simple funciton that returns the state the tableau is in vis-a-vis the Aarosnon/Gottesman decomposition algorithm.






