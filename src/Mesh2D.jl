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
struct Mesh{IntType,FloatType}
    nx::IntType
    ny::IntType
    node_map::Vector{IntType}
    dx::FloatType
    dy::FloatType
    x0::FloatType
    y0::FloatType

    Mesh{IT,FT}(nx::IT, ny::IT, dx::FT, dy::FT, x0::FT, y0::FT) where {IT<:Integer,FT<:AbstractFloat} =
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
    Ai_with_adjacencies = Tuple{typeof(nx),Vector{typeof(nx)},Int}[]
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

mutable struct NodeMap{IntType}
    mapped::Dict{IntType,IntType}
    reference_to_real::Vector{IntType}
    max_node_examined::IntType

    NodeMap{IT}(node_map) where {IT<:Integer} = new(Dict{IT,IT}(), node_map, zero(IT))
end

"""
Given ni, a real node index (i.e. in the reordered order), look up the index of its reference
node.
"""
function lookup!(node_map::NodeMap, ni)
    if ni in keys(node_map.mapped)
        pop!(node_map.mapped, ni)
    else
        next_to_examine = node_map.max_node_examined + 1
        while node_map.reference_to_real[next_to_examine] != ni
            node_map.mapped[node_map.reference_to_real[next_to_examine]] = next_to_examine
            next_to_examine += 1
        end
        node_map.max_node_examined = next_to_examine
    end
end

"""
sparsity(mesh)

Returns just the sparsity pattern, i.e. the I, J arrays for a CSR format matrix,
given the mesh. The returned sparsity has only the lower triangular part, since
the stiffness matrix is symmetric.

This accounts for the 8 degrees of freedom collocated at each node.
"""
function sparsity(mesh::Mesh{IntType}) where {IntType}
    I = [one(IntType)]
    J = IntType[]
    scratch_space = Vector{IntType}[]
    real_to_reference = NodeMap{IntType}(mesh.node_map)

    for real_index in 1:(mesh.nx+1)*(mesh.ny+1)
        ni = lookup!(real_to_reference, real_index)
        # Get the adjacencies for this node.
        Ai = get_adjacent_nodes(ni, mesh.nx, mesh.ny, scratch_space)
        # Map adjacent nodes to their reordered indices
        for i in eachindex(Ai)
            Ai[i] = mesh.node_map[Ai[i]]
        end
        # Filter previously considered nodes (we create only the lower triangular part,
        # or upper triangular if in CSR format). 
        filter!(>(real_index), Ai)
        # Sort in ascending order so that indices are in the correct order.
        sort!(Ai)

        outer_first_dof = (real_index - 1) * 8 + 1
        outer_last_dof = outer_first_dof + 7
        for i in outer_first_dof:outer_last_dof
            # Indices for degrees of freedom associated to this node.
            for j in i:outer_last_dof
                push!(J, j)
            end

            # Indices for each adjacent node
            for nj in Ai
                inner_first_dof = (nj - 1) * 8 + 1
                inner_last_dof = inner_first_dof + 7
                for j in inner_first_dof:inner_last_dof
                    @assert j > i
                    push!(J, j)
                end
            end

            # Update the colptr array
            push!(I, length(J) + 1)
        end
        push!(scratch_space, Ai)
    end
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
