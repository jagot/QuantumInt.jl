abstract AbstractObservable
using SimpleFields

import Base: call

call{T<:AbstractFloat}(obs::AbstractObservable,
                       Ψ, i::Integer, τ::T, field::Field) = (obs.v[i] = value(obs, Ψ, i, τ, field))

type TimeObservable{T<:AbstractFloat} <: AbstractObservable
    v::Vector{T}
    TimeObservable(N::Int) = new(Vector{T}(N))
end
value{T<:AbstractFloat}(obs::TimeObservable{T},
                        Ψ, i::Integer, τ::T, field::Field) = (i-1)*τ

type NormObservable{T<:AbstractFloat} <: AbstractObservable
    v::Vector{T}
    NormObservable(N::Int) = new(Vector{T}(N))
end
value{T<:AbstractFloat}(obs::NormObservable{T},
                        Ψ, i::Integer, τ::T, field::Field) = norm(Ψ)

type FieldObservable{T<:AbstractFloat} <: AbstractObservable
    v::Vector{T}
    FieldObservable(N::Int) = new(Vector{T}(N))
end
value{T<:AbstractFloat}(obs::FieldObservable{T},
                        Ψ, i::Integer, τ::T, field::Field) = field((i-1)*τ/field.T)

const observable_types = Dict("time" => TimeObservable,
                              "norm" => NormObservable,
                              "field" => FieldObservable)
