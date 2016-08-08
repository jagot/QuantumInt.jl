using Magnus
using SimpleFields

# Calculation in the dipole approximation, linear
# polarization. Integration is done in the eigenbasis of the atomic
# Hamiltonian.
function calc{M<:AbstractMatrix{Float64}}(observe::Function,
                                          E::Vector{Vector{Float64}},
                                          V::Vector{M},
                                          D::Operator,
                                          field::Field, ndt::Integer;
                                          cutoff::Real = 0,
                                          magnus_kwargs...)
    npartial = length(E)
    H₀,D,gst = hamiltonian(E, V, D, cutoff)

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

calc{M<:AbstractMatrix{Float64}}(E::Vector{Vector{Float64}},
                                 V::Vector{M},
                                 D::Operator,
                                 field::Field, ndt::Integer;
                                 kwargs...) =
                                     calc((Ψ,i,τ) -> (),
                                          E, V, D,
                                          field, ndt;
                                          kwargs...)

function calc{M<:AbstractMatrix{Float64}}(E::Vector{Vector{Float64}},
                                          V::Vector{M},
                                          D::Operator,
                                          field::Field, ndt::Integer,
                                          observables::Vector{Symbol},
                                          extra_observables::Function = (Ψ,i,τ,field) -> ();
                                          kwargs...)

    N = ceil(Int, ndt*field.tmax)

    obss = []
    for o in observables
        o in keys(observable_types) || error("Unknown observable, $(t)")
        push!(obss, (o,observable_types[o](N)))
    end

    E, Ψ = calc(E, V, D,
                field, ndt;
                kwargs...) do Ψ,i,τ
                    for o in obss
                        o[2](Ψ,i,τ,field)
                    end
                    extra_observables(Ψ,i,τ,field)
                end
    E, Ψ, obss
end

export calc
