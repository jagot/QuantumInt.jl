mask(dE,dEmax) = cospi(dE/(2dEmax))
function mask_v(E, cutoff, i, j)
    dE = E[i:(j-1)] - cutoff
    mask(dE, maximum(E)-cutoff)
end

function mask(H, npartial, cutoff, power = 1.0)
    cutoff /= 27.211
    E = diag(H)
    i = 1
    T = eltype(E)
    cutoff = T(cutoff)
    m = ones(T, length(E))
    for k = 1:npartial-1
        i = findnext(e -> e > cutoff, E, i)
        j = findnext(e -> e <= cutoff, E, i)
        j == 0 && (j=length(E)+1)
        m[i:(j-1)] = mask_v(E, cutoff, i, j)
        i = j
    end
    m[i:end] = 0
    m.^T(power)
end

# This function generates a diagonal matrix which, when applied to a
# wavefunction, will eat up all energies larger than cutoff, and kill
# the highest partial wave completely.
function Egobbler{T<:AbstractFloat}(H::AbstractMatrix{T}, npartial, cutoff, mode)
    if cutoff == 0
        Ψ -> ()
    else
        d = speye(T,size(H,1)) - sparse(Diagonal(mask(H, npartial, cutoff, T(1/8))))
        if mode == :gpu && Magnus.cuda
            d = Magnus.upload(d*one(Complex{T}))
        end
        Ψ -> A_mul_B!(-one(Complex{T}), d, Ψ, one(Complex{T}), Ψ)
    end
end
