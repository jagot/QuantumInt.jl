# Subspace matrix describing the coupling between partial waves k1 and
# k2.
type SubOperator
    m::AbstractMatrix
    k1::Integer
    k2::Integer
end
type Operator
    o::Vector{SubOperator}
end

export SubOperator, Operator
