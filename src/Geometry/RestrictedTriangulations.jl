"""
"""
struct RestrictedTriangulation{Dc,Dp,G,A} <: Triangulation{Dc,Dp}
  parent_trian::G
  cell_to_parent_cell::A
  @doc """
  """
  function RestrictedTriangulation(
    parent_trian::Triangulation,
    cell_to_parent_cell::AbstractVector{<:Integer})

    Dc = num_cell_dims(parent_trian)
    Dp = num_point_dims(parent_trian)
    G = typeof(parent_trian)
    A = typeof(cell_to_parent_cell)
    new{Dc,Dp,G,A}(parent_trian,cell_to_parent_cell)
  end
end

# Constructors

function RestrictedTriangulation(
  parent_trian::Triangulation,
  parent_cell_to_mask::AbstractArray{Bool})

  cell_to_parent_cell = findall(collect1d(parent_cell_to_mask))
  RestrictedTriangulation(parent_trian,cell_to_parent_cell)
end

function RestrictedTriangulation(
  parent_trian::Triangulation,
  parent_cell_to_mask::AbstractVector{Bool})

  cell_to_parent_cell = findall(parent_cell_to_mask)
  RestrictedTriangulation(parent_trian,cell_to_parent_cell)
end

function Triangulation(
  parent_trian::Triangulation,
  cell_to_parent_cell::AbstractVector{<:Integer})

  RestrictedTriangulation(parent_trian,cell_to_parent_cell)
end

function Triangulation(
  parent_trian::Triangulation,
  parent_cell_to_mask::AbstractArray{Bool})

  RestrictedTriangulation(parent_trian,parent_cell_to_mask)
end

function Triangulation(
  parent_model::DiscreteModel,
  cell_to_parent_cell::AbstractVector{<:Integer})

  parent_trian = Triangulation(parent_model)
  RestrictedTriangulation(parent_trian,cell_to_parent_cell)
end

function Triangulation(
  parent_model::DiscreteModel,
  parent_cell_to_mask::AbstractArray{Bool})

  parent_trian = Triangulation(parent_model)
  RestrictedTriangulation(parent_trian,parent_cell_to_mask)
end

function Triangulation(parent_model::DiscreteModel, labels::FaceLabeling; tags)
  parent_trian = Triangulation(parent_model)
  parent_cell_to_mask = get_face_mask(labels,tags,num_cell_dims(parent_model))
  RestrictedTriangulation(parent_trian,parent_cell_to_mask)
end

# Triangulation API


TriangulationStyle(::Type{<:RestrictedTriangulation}) = SubTriangulation()

get_background_triangulation(trian::RestrictedTriangulation) = get_background_triangulation(trian.parent_trian)

get_node_coordinates(trian::RestrictedTriangulation) = get_node_coordinates(trian.parent_trian)

get_reffes(trian::RestrictedTriangulation) = get_reffes(trian.parent_trian)

#GridTopology(trian::RestrictedTriangulation) = GridTopology(trian.parent_trian)
function GridTopology(trian::RestrictedTriangulation)
  parent_trian = trian.parent_trian
  parent_cell_to_vertices, parent_vertex_to_node, = _generate_cell_to_vertices_from_grid(parent_trian)
  cell_to_vertices = Table(lazy_map(Reindex(parent_cell_to_vertices),trian.cell_to_parent_cell))
  vertex_to_node = parent_vertex_to_node

  @notimplementedif (! is_regular(parent_trian)) "Extrtacting the GridTopology form a Grid only implemented for the regular case"
  node_to_coords = get_node_coordinates(trian)
  if vertex_to_node == 1:num_nodes(trian)
    vertex_to_coords = node_to_coords
  else
    vertex_to_coords = node_to_coords[vertex_to_node]
  end
  cell_to_type = collect(get_cell_type(trian))
  polytopes = map(get_polytope, get_reffes(trian))

  UnstructuredGridTopology(vertex_to_coords, cell_to_vertices, cell_to_type, polytopes, OrientationStyle(parent_trian))
end

function get_cell_coordinates(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_coordinates(trian.parent_trian)
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_type(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_type(trian.parent_trian)
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_node_ids(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_node_ids(trian.parent_trian)
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_map(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_map(trian.parent_trian)
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_to_bgcell(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_to_bgcell(trian.parent_trian)
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_facet_normal(trian::RestrictedTriangulation)
  parent_cell_data = get_facet_normal(trian.parent_trian)
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end

function get_cell_ref_map(trian::RestrictedTriangulation)
  parent_cell_data = get_cell_ref_map(trian.parent_trian)
  lazy_map(Reindex(parent_cell_data),trian.cell_to_parent_cell)
end
