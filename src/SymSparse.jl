module SymSparse

# Parametric so we can use with Int32 or Int64 types.
struct SymbolicRow{T}
    nnz::Int
    cols::Vector{T}
end

end # module SymSparse