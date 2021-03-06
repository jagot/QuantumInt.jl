using Magnus
using SimpleFields

# Calculation in the dipole approximation, linear
# polarization.
function calc{T<:AbstractFloat}(H₀::AbstractMatrix{T},
                                D::AbstractMatrix,
                                ψ₀::Vector{Complex{T}},
                                field::Field, ndt::Integer;
                                gobble!::Function = Ψ -> (),
                                observe::Function = (Ψ,i,τ,field) -> (),
                                observables::Vector{Symbol} = Symbol[],
                                mode = :cpu,
                                verbose = true,
                                magnus_kwargs...)
    if verbose
        H = H₀ + D
        @printf("Hamiltonian dimensions: %dx%d, sparsity: %03.3f%%\n",
                size(H,1), size(H,2),
                100(1.0-length(nonzeros(H))/length(H)))
    end

    H₀ *= one(Complex{T})
    D *= one(Complex{T})

    N = ceil(Int, ndt*field.tmax)

    f = t -> field(t/field.T)

    obss = Dict{AbstractString,Any}()
    for o in observables
        os = string(o)
        os in keys(observable_types) || error("Unknown observable, $(os)")
        obss[os] = observable_types[os]{typeof(field.tmax)}(N)
    end

    results = integrate(ψ₀, T(field.tmax*field.T), N,
                        H₀, f, D, -im;
                        mode = mode,
                        verbose = verbose,
                        magnus_kwargs...) do Ψ,i,τ
                            gobble!(Ψ)
                            observe(Ψ,i,τ,field)
                            for o in values(obss)
                                o(Ψ,i,τ,field)
                            end
                        end

    Dict("E" => real(diag(H₀)), "psi" => results[:V],
         "milliseconds" => results[:milliseconds],
         "performance" => results[:performance],
         [k => v.v for (k,v) in obss]...)
end


# Calculation in the dipole approximation, linear
# polarization. Integration is done in the eigenbasis of the atomic
# Hamiltonian.
function calc{T<:AbstractFloat,
              M<:AbstractMatrix}(E::Vector{Vector{T}},
                                 V::Vector{M},
                                 D::Operator,
                                 field::Field, ndt::Integer;
                                 observe::Function = (Ψ,i,τ,field) -> (),
                                 observables::Vector{Symbol} = Symbol[],
                                 cutoff::Real = 0,
                                 mask_ratio::Real = 0,
                                 mode = :cpu,
                                 verbose = true,
                                 magnus_kwargs...)
    npartial = length(E)
    H₀,D,gst = hamiltonian(E, V, D, cutoff)
    ψ₀ = zeros(Complex{T}, size(H₀,1))
    ψ₀[gst] = 1.0
    gobble! = Egobbler(H₀, npartial, mask_ratio*cutoff, mode)

    calc(H₀, D, ψ₀, field, ndt;
         gobble! = gobble!,
         observe = observe,
         observables = observables,
         mode = mode, verbose = verbose,
         magnus_kwargs...)
end

export calc
