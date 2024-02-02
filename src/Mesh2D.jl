module Mesh2D

using DataStructures: Queue, enqueue!, dequeue!

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
function get_adjacent_nodes(ni, nx, ny)
    # row and col are 0-based
    row = (ni - 1) ÷ (nx + 1)
    col = (ni - 1) % (nx + 1)

    adjacent = Vector{typeof(ni)}(undef, 8)
    # ccw from bottom left
    adjacent[1] = (ni - (nx + 1) - 1) * (col > 0)
    adjacent[2] = ni - (nx + 1)
    adjacent[3] = (ni - (nx + 1) + 1) * (col < nx)
    adjacent[4] = (ni + 1) * (col < nx)
    adjacent[5] = (ni + (nx + 1) + 1) * (col < nx)
    adjacent[6] = ni + (nx + 1)
    adjacent[7] = (ni + (nx + 1) - 1) * (col > 0)
    adjacent[8] = (ni - 1) * (col > 0)

    filter!(i -> i > 0 && i <= (nx + 1) * (ny + 1), adjacent)
end

"""
Assign nodes in reverse Cuthill-McKee ordering.
"""
function renumber_nodes(nx, ny)
    node_map = zeros(typeof(nx), (nx + 1) * (ny + 1))
    node_map[1] = 1
    i = 2
    node_queue = Queue{Int}()
    enqueue!(node_queue, 1)

    max_neighbor = n -> begin
        adj = get_adjacent_nodes(n, nx, ny)
        (maximum(@view(node_map[adj])), length(adj))
    end

    while !isempty(node_queue)
        ni = dequeue!(node_queue)
        adjacent = get_adjacent_nodes(ni, nx, ny)
        filter!(n -> node_map[n] == 0, adjacent)
        sort!(adjacent, by=max_neighbor, rev=true)
        for n in adjacent
            node_map[n] = i
            i += 1
            enqueue!(node_queue, n)
        end
    end

    node_map
end

end # module Mesh2D
