module Mesh2D

"""
The information needed to represent a structured mesh of rectangular elements. An instance
represents a grid of nx x ny rectangular elements with dimensions dx x dy.

Convention: element #1 is the bottom left and then numbering goes from left to right, bottom
to top. The base coordinates x0 and y0 are the bottom left corner of this element.
Reference node numbers start at the bottom left and run from left to right, bottom to top.

`node_map[i]` is the true node number (i.e. used for finite element DOFs) for reference node
number `i`. These are computed as a reverse Cuthill-McKee ordering.
"""
struct Mesh{IntType, FloatType}
    nx::IntType
    ny::IntType
    node_map::Vector{IntType}
    dx::FloatType
    dy::FloatType
    x0::FloatType
    y0::FloatType

    Mesh{IT, FT}(nx::IT, ny::IT, dx::FT, dy::FT, x0::FT, y0::FT) where {IT <: Integer, FT <: AbstractFloat} =
        new(nx, ny, renumber_nodes(nx, ny), dx, dy, x0, y0)
end

numelements(mesh::Mesh) = mesh.nx * mesh.ny
numnodes(mesh::Mesh) = (mesh.nx + 1) * (mesh.ny + 1)

"""
Return the node numbers in the mesh of each of the 4 vertices of given element.
Ordering is starting with bottom left and going counterclockwise.
"""
function element_reference_nodes(nx, el_index)
    # row and col are 0-based
    row = (el_index - 1) / nx
    col = (el_index - 1) % nx
    bot_left = row * (nx + 1) + col
    top_left = bot_left + nx + 1
    (bot_left, bot_left + 1, top_left - 1, top_left)
end

"""
Get all nodes that share an element with node 'ni' given mesh dimensions nx and ny (in elements)
"""
function get_adjacent_nodes(ni, nx, ny, scratch_space)
    # row and col are 0-based
    row = (ni - 1) ÷ (nx + 1)
    col = (ni - 1) % (nx + 1)

    if isempty(scratch_space)
        adjacent = Vector{typeof(ni)}(undef, 8)
    else
        adjacent = pop!(scratch_space)
        resize!(adjacent, 8)
    end
    # ccw from bottom left
    @inbounds begin
        adjacent[1] = (ni - (nx + 1) - 1) * (col > 0)
        adjacent[2] = ni - (nx + 1)
        adjacent[3] = (ni - (nx + 1) + 1) * (col < nx)
        adjacent[4] = (ni + 1) * (col < nx)
        adjacent[5] = (ni + (nx + 1) + 1) * (col < nx)
        adjacent[6] = ni + (nx + 1)
        adjacent[7] = (ni + (nx + 1) - 1) * (col > 0)
        adjacent[8] = (ni - 1) * (col > 0)
    end

    filter!(i -> i > 0 && i <= (nx + 1) * (ny + 1), adjacent)
end

"""
Assign nodes in reverse Cuthill-McKee ordering.
"""
function renumber_nodes(nx, ny)
    nodes = [one(nx)]
    scratch_space = Vector{typeof(nx)}[]
    Ai_with_adjacencies = Tuple{typeof(nx), Vector{typeof(nx)}, Int}[]
    node_map = zeros(typeof(nx), (nx + 1) * (ny + 1))
    node_map[1] = 1

    compare_nodes = (tup1, tup2) -> begin
        if tup1[2] > tup2[2]
            return true
        elseif tup1[2] < tup2[2]
            return false
        elseif tup1[3] < tup2[3]
            return true
        elseif tup1[3] > tup2[3]
            return false
        else
            return tup1[1] < tup2[1]
        end
    end

    for i in 1:length(node_map)
        empty!(Ai_with_adjacencies)
        Ai = filter!(n -> node_map[n] == 0, get_adjacent_nodes(nodes[i], nx, ny, scratch_space))
        for n in Ai
            Aj = get_adjacent_nodes(n, nx, ny, scratch_space)
            degree = length(Aj)
            @inbounds for j in eachindex(Aj)
                Aj[j] = node_map[Aj[j]]
            end
            filter!(!=(0), Aj)
            sort!(Aj, rev=true)
            push!(Ai_with_adjacencies, (n, Aj, degree))
        end
        push!(scratch_space, Ai)
        sort!(Ai_with_adjacencies, lt=compare_nodes)
        for (n, A, _) in Ai_with_adjacencies
            push!(scratch_space, A)
            push!(nodes, n)
            node_map[n] = length(nodes)
        end
        length(nodes) == length(node_map) && break
    end
    node_map
end

"""
sparsity(mesh)

Returns just the sparsity pattern, i.e. the I, J arrays for a CSR format matrix,
given the mesh. The returned sparsity has only the lower triangular part, since
the stiffness matrix is symmetric.

This accounts for the 8 degrees of freedom collocated at each node.
"""
function sparsity(mesh::Mesh{IntType}) where IntType
    IJ = Tuple{IntType, IntType}[]
    scratch_space = Vector{IntType}[]
    for ni in 1:(mesh.nx + 1) * (mesh.ny + 1)
        i = (mesh.node_map[ni] - 1) * 8 + 1
        for j in i:i+7
            for k in j:i+7
                push!(IJ, (k, j))
            end
        end
        Ai = get_adjacent_nodes(ni, mesh.nx, mesh.ny, scratch_space)
        for nj in Ai
            j = (mesh.node_map[nj] - 1) * 8 + 1
            if i > j
                for k in i:i+7
                    for l in j:j+7
                        push!(IJ, (k, l))
                    end
                end
            end
        end
    end
    
    sort!(IJ)
    I = IntType[]
    J = IntType[]
    first = 1
    second = 2
    n = length(IJ)
    while first <= n
        # Advance the second pointer
        while second <= n && IJ[second] == IJ[first]
            second += 1
        end
        push!(J, IJ[first][2])
        if isempty(I) || IJ[first][1] != IJ[first-1][1]
            push!(I, length(J))
        end
        first = second
    end
    push!(I, length(J) + 1)
    I, J
end

"""
stencil(mesh)

Return a 528 * num_elements array. For each element the column here has the index
in the assembled sparse matrix for each lower triangular entry of the stiffness
matrix. So if `V` is the array of non-zero values for the stiffness matrix, we can
add element `i`'s contribution by: `V[stencil[:, i]] .+= K_map * λ`.
"""

end # module Mesh2D
