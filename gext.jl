module GExt
""" Perform extrapolation of density matrices using their geometric structure with
  least-squares method.

  nbas: number of basis functions
  nocc: number of occupied orbitals (half the number of electrons)
  nmat: number of elements used to perform the extrapolation plus one (the first element is
        supposed to contain the point we use as a reference, and we do not need its
        associated positions)
  nqm:  number of QM atoms
"""

using LinearAlgebra, DelimitedFiles, Printf
# To make sure the guesses we write are precise enough
Base.show(io::IO, f::Float64) = @printf(io, "%24.16e", f)

export gext

function gext(; ε=1e-5)
  """ Main function which does the extrapolation """
  cmat = get_cmat("cmat.dat")
  smat = get_smat("smat.dat")
  hpos = get_hpos("hpos.dat")

  # Make sure that we have the correct number of elements
  @assert size(smat, 2) == size(cmat, 3) == size(hpos, 3)
  nbas, nocc, nmat = size(cmat)

  # Perform orthonormalisation
  C = Array{Float64}(undef, nbas, nocc, nmat)
  for i in 1:nmat
    C[:, :, i] = sqrt(tril_to_sym(smat[:, i]))*cmat[:, :, i]
  end

  # The reference point
  γ_ref = C[:, :, 1]

  # Project the points on the tangent space at γ_ref
  Γs = logGrs(γ_ref, C)

  # Create a matrix of previous descriptors
  poly = [get_coulomb(hpos[:, :, imat]) for imat in 1:nmat-1]
  P = [poly[i][j] for i in eachindex(poly), j in eachindex(poly[1])]

  # Descriptor of the current point
  poly_test = get_coulomb(hpos[:, :, nmat])

  # Optional stablisation
  ds = hcat(poly_test, zeros(nmat-1))
  Ps = hcat(P, ε*Matrix(I, nmat-1, nmat-1))

  # Least-squares on the Coulomb descriptors
  α = ds'*pinv(Ps)

  # Extrapolation on the tangent space
  Γguess = sum(α[i]*Γs[:, :, i+1] for i in 1:nmat-1)
  # Retract the point back to the Grassmannian
  Cguess = expGr(γ_ref, Γguess)

  # Compute the idempotent guess
  pguess = 2*Cguess*Cguess'
  lguess = vec_tril(pguess)

  # Write the guess on file
  writedlm("guess", lguess)
end

function get_coulomb(hpos)
  """ Return the Coulomb descriptor, which is the lower triangular part of the symmetric
  Coulomb matrix """
  charges = hpos[1, :]
  positions = hpos[2:4, :]

  nqm = size(positions, 2)
  G = Matrix{Float64}(undef, nqm, nqm)
  for i in 1:nqm
    for j in 1:nqm
      if i == j
        G[i, i] = 0.5*charges[i]^(2.4)
      else
        G[i, j] = charges[i]*charges[j]/norm(positions[:, i] - positions[:, j])
      end
    end
  end

  return vec_tril(G)
end

#############################################
#### Exponential and logarithm functions ####
#############################################

function expGr(γ_ref, Γ)
  """ Perform the exponential of the tangent vector Γ at base point γ_ref """
  d = svd(Γ)
  Y = γ_ref * d.V * cos(Diagonal(d.S)) + (d.U) * sin(Diagonal(d.S))
  # Optional, to make sure that rounding errors do not prevent Y from being orthogonal
  return Matrix(qr(Y).Q)
end

function logGr(γ_ref, γ)
  """ Perform the logarithm of the point γ at base point γ_ref """
  Z = γ'*γ_ref
  W = Z \ (γ'-Z*γ_ref')
  d = svd(W, full = false)
  A = (d.V)*atan(Diagonal(svd(W, full = false).S))*(d.U)'
  return A
end

function logGrs(γ_ref, γs)
  """ Perfor the logarithms on a list of points """
  nbas, nocc, nmat = size(cmat)

  Γs = Array{Float64}(undef, nbas, nocc, nmat)
  for iγ in 1:nmat
    Γs[:, :, iγ] = logGr(γ_ref, γs[:, :, iγ])
  end

  return Γs
end

#############################################
####           Helper functions          ####
#############################################

function tril_to_sym(v)
  """ Convert a vector representing the lower triangular part of a symmetric matrix to
  a symmetric matrix """
  n = Int(floor(sqrt(2*length(v))))
  M = zeros(n, n)

  k = 0
  for i in 1:n
    for j in 1:i
      k += 1
      M[i, j] = v[k]
      M[j, i] = v[k]
    end
  end

  return M
end

function vec_tril(M)
  """ Extract the the lower triangular part of a square matrix """
  m, n = size(M)
  @assert m == n
  l = n*(n+1) ÷ 2
  v = zeros(l)

  k = 0
  for i in 1:n
    for j in 1:i
      v[k + j] = M[i, j]
    end
    k += i
  end

  return v
end

#############################################
####      Data processing functions      ####
#############################################

function get_cmat(file)
  """ Placeholder for a function that returns an array cmat of size nbas×nocc×nmat
  containing molecular orbitals matrices """
end

function get_smat(file)
  """ Placeholder for a function that returns an array smat of size nbas×(nbas+1)/2×nmat
  containing lower diagonal of overlap matrices """
end

function get_hpos(file)
  """ Placeholder for a function that returns an array smat of size 4×nqm×nmat
  containing charges and positions of the QM atoms for the previous positions as well as
  the current for which we want the extrapolation """
end

end
