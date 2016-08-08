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
                                          mode = :cpu,
                                          magnus_kwargs...)
    npartial = length(E)
    H₀,D,gst = hamiltonian(E, V, D, cutoff)

    if mode == :gpu
        H₀ *= (1.0 + 0im)
        D *= (1.0 + 0im)
    end

    ψ₀ = zeros(Complex128, size(H₀,1))
    ψ₀[gst] = 1.0

    N = ceil(Int, ndt*field.tmax)

    f = t -> field(t/field.T)

    Ψ = integrate(ψ₀, field.tmax*field.T, N,
                  H₀, f, D, -im;
                  mode = mode,
                  magnus_kwargs...) do Ψ,i,τ
                      observe(Ψ,i,τ,field)
                  end

    diag(H₀), Ψ
end

calc{M<:AbstractMatrix{Float64}}(E::Vector{Vector{Float64}},
                                 V::Vector{M},
                                 D::Operator,
                                 field::Field, ndt::Integer;
                                 kwargs...) =
                                     calc((Ψ,i,τ,field) -> (),
                                          E, V, D,
                                          field, ndt;
                                          kwargs...)

function calc{M<:AbstractMatrix{Float64}}(observe::Function,
                                          E::Vector{Vector{Float64}},
                                          V::Vector{M},
                                          D::Operator,
                                          field::Field, ndt::Integer,
                                          observables::Vector{Symbol};
                                          kwargs...)

    N = ceil(Int, ndt*field.tmax)

    obss = []
    for o in observables
        o in keys(observable_types) || error("Unknown observable, $(t)")
        push!(obss, (o,observable_types[o](N)))
    end

    E, Ψ = calc(E, V, D,
                field, ndt;
                kwargs...) do Ψ,i,τ,field
                    observe(Ψ,i,τ,field)
                    for o in obss
                        o[2](Ψ,i,τ,field)
                    end
                end
    E, Ψ, obss
end

calc{M<:AbstractMatrix{Float64}}(E::Vector{Vector{Float64}},
                                 V::Vector{M},
                                 D::Operator,
                                 field::Field, ndt::Integer,
                                 observables::Vector{Symbol};
                                 kwargs...) =
                                     calc((Ψ,i,τ,field) -> (),
                                          E, V, D,
                                          field, ndt,
                                          observables;
                                          kwargs...)

export calc
