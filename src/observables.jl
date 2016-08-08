abstract AbstractObservable
using SimpleFields

import Base: call

type TimeObservable <: AbstractObservable
    v::Vector{Float64}
    TimeObservable(N::Int) = new(Vector{Float64}(N))
end
call(obs::TimeObservable, Ψ, i, τ, field::Field) = (obs.v[i] = (i-1)*τ)

type NormObservable <: AbstractObservable
    v::Vector{Float64}
    NormObservable(N::Int) = new(Vector{Float64}(N))
end
call(obs::NormObservable, Ψ, i, τ, field::Field) = vecnorm(Ψ)

type FieldObservable <: AbstractObservable
    v::Vector{Float64}
    FieldObservable(N::Int) = new(Vector{Float64}(N))
end
call(obs::FieldObservable, Ψ, i, τ, field::Field) = field((i-1)*τ)

const observable_types = Dict(:time => TimeObservable,
                              :norm => NormObservable,
                              :field => FieldObservable)
