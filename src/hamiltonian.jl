"""Given a matrix of energies `E` (in Hartrees), where each column
corresponds to a different angular symmetry, a vector of matrices `V`
containing the basis vectors for each angular symmetry, generate a
field-free Hamiltonian on diagonal form and transform all operators
`O` to the field-free basis. Optionally, `cutoff` (in eV) limits the
spectrum. Also returned is the index of the lowest, i.e., groundstate
energy."""
function hamiltonian{M<:AbstractMatrix{Float64}}(E::Vector{Vector{Float64}},
                                                 V::Vector{M},
                                                 O::Vector{Operator},
                                                 cutoff::Real = 0,
                                                 verbose::Bool = false)
    n = length(E)
    sel = [(1:length(E[j])) for j = 1:n]

    # Find range of elements with energies below the cutoff energy.
    if cutoff != 0
        cutoff /= 27.211
        for j = 1:n
            i = findlast(e -> e <= cutoff, E[j])
            verbose && i == sel[j][end] && warn("Largest energy $(E[j][end]) smaller than cutoff $(cutoff).")
            i == 0 && error("No energy smaller than cutoff in partial wave $j.")
            sel[j] = 1:i
        end
    end

    # Find block coordinates
    coords = Vector{UnitRange{Int64}}(n)
    coords[1] = sel[1]
    for j = 2:n
        coords[j] = sel[j] + coords[j-1][end]
    end

    # Reduce the space
    E = vcat([E[j][sel[j]] for j = 1:n]...)
    V = Matrix{Float64}[V[j][:,sel[j]] for j in 1:n]

    # Field-free Hamiltonian in its eigenbasis
    H = sparse(Diagonal(E))

    # Transform all operators O into the eigenbasis
    O = map(O) do o
        ot = spzeros(eltype(first(o.o).m), size(H)...)
        for oo in o.o
            o12 = sparse(V[oo.k1]'*oo.m*V[oo.k2])
            ot[coords[oo.k1],coords[oo.k2]] += o12
            ot[coords[oo.k2],coords[oo.k1]] += o12'
        end
        ot
    end

    # Find groundstate
    gst = indmin(E)

    H,O,gst
end

function hamiltonian{M<:AbstractMatrix{Float64}}(E::Vector{Vector{Float64}},
                                                 V::Vector{M},
                                                 O::Operator,
                                                 cutoff::Real = 0,
                                                 verbose::Bool = false)
    H,O,gst = hamiltonian(E, V, [O], cutoff, verbose)
    H,O[1],gst
end

export hamiltonian
