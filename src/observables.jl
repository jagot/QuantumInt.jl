abstract AbstractObservable
using SimpleFields

import Base: call

call(obs::AbstractObservable, Ψ, i, τ, field::Field) = (obs.v[i] = value(obs, Ψ, i, τ, field))

type TimeObservable <: AbstractObservable
    v::Vector{Float64}
    TimeObservable(N::Int) = new(Vector{Float64}(N))
end
value(obs::TimeObservable, Ψ, i, τ, field::Field) = (i-1)*τ

type NormObservable <: AbstractObservable
    v::Vector{Float64}
    NormObservable(N::Int) = new(Vector{Float64}(N))
end
value(obs::NormObservable, Ψ, i, τ, field::Field) = norm(Ψ)

type FieldObservable <: AbstractObservable
    v::Vector{Float64}
    FieldObservable(N::Int) = new(Vector{Float64}(N))
end
value(obs::FieldObservable, Ψ, i, τ, field::Field) = field((i-1)*τ/field.T)

const observable_types = Dict("time" => TimeObservable,
                              "norm" => NormObservable,
                              "field" => FieldObservable)
