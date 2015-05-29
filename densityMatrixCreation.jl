
#creates a normalised vector
#representing n qubits
#TODO some type of low rank vector?

function randomNormalisedVector(n)
  size = 2^n # so an n bit qubit needs 2^n vectors to represent it.
  preNorm = randn(size)+rand(size)*im
  #I think this all we need to normalise it, i.e all the probablilities sum to 1.
  return preNorm/norm(preNorm)
end

