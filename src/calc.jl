using Magnus
using SimpleFields

# Calculation in the dipole approximation, linear
# polarization. Integration is done in the eigenbasis of the atomic
# Hamiltonian.
function calc{M<:AbstractMatrix{Float64}}(observe::Function,
                                          E::Vector{Vector{Float64}},
                                          V::Vector{M},
                                          D::Operator,
                                          field::Field, ndt::Integer,
                                          cutoff::Real = 0;
                                          magnus_kwargs...)
    npartial = length(E)
    H₀,D,gst = hamiltonian(E, V, D, cutoff, verbose)

    ψ₀ = zeros(Complex128, size(H₀,1))
    ψ₀[gst] = 1.0

    N = ceil(Int, ndt*field.tmax)

    f = t -> field(t/field.T)

    Ψ = integrate(observe,
                  ψ₀, field.tmax*field.T, N,
                  H₀, f, D, -im;
                  magnus_kwargs...)

    diag(H₀), Ψ
end
