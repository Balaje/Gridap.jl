module DivConformingFESpacesTests

using Test
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Io

function test_div_v_q_equiv(U,V,P,Q,Ω)
  v=get_fe_basis(V)
  u=get_trial_fe_basis(U)

  q=get_fe_basis(Q)
  p=get_trial_fe_basis(P)

  dΩ=Measure(Ω,1)
  dΩᵣ=Measure(Ω,1,ReferenceDomain())

  a1(p,v)=∫(divergence(v)*p)dΩ
  a2(p,v)=∫(DIV(v)*p)dΩᵣ

  tol=1.0e-12
  assem = SparseMatrixAssembler(P,V)
  data = collect_cell_matrix(P,V,a1(p,v))
  A1 = assemble_matrix(assem,data)
  data = collect_cell_matrix(P,V,a2(p,v))
  A2 = assemble_matrix(assem,data)
  @test norm(A1-A2) < tol

  a3(u,q)=∫(q*divergence(u))dΩ
  a4(u,q)=∫(q*DIV(u))dΩᵣ
  assem = SparseMatrixAssembler(U,Q)
  data = collect_cell_matrix(U,Q,a3(u,q))
  A3 = assemble_matrix(assem,data)
  data = collect_cell_matrix(U,Q,a4(u,q))
  A4 = assemble_matrix(assem,data)
  @test norm(A3-A4) < tol
end

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 1

u(x) = x

reffe = ReferenceFE(raviart_thomas,order)

V = TestFESpace(model,reffe,dirichlet_tags = [1,6])
test_single_field_fe_space(V)
U = TrialFESpace(V,u)

reffe = ReferenceFE(lagrangian,Float64,order)
Q = TestFESpace(model,reffe,conformity=:L2)
P = TrialFESpace(Q)

uh = interpolate(u,U)

e = u - uh

Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)

el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

test_div_v_q_equiv(U,V,P,Q,Ω)

#using Gridap.Visualization
#
#writevtk(Ω,"trian",nsubcells=10,cellfields=["uh"=>uh])

order = 1

reffe = ReferenceFE(TET,raviart_thomas,order)

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = simplexify(CartesianDiscreteModel(domain,partition))

labels = get_face_labeling(model)
dir_tags = Array{Integer}(undef,0)

V = FESpace(model,reffe,conformity=DivConformity())
U = TrialFESpace(V,u)

reffe = ReferenceFE(lagrangian,Float64,order)
Q = TestFESpace(model,reffe,conformity=:L2)
P = TrialFESpace(Q)

v(x) = VectorValue(-0.5*x[1]+1.0,-0.5*x[2],-0.5*x[3])
vh = interpolate(v,V)

e = v - vh

Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)

el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

test_div_v_q_equiv(U,V,P,Q,Ω)

#using Gridap.Visualization
#
#writevtk(trian,"test",order=3,cellfields=["vh"=>vh])

# Tests on manifold

# Create domain
domain = (0,1,0,1,0,1)
cells  = (2,2,2)
model  = CartesianDiscreteModel(domain,cells)

# Restrict model to cube surface (using new BoundaryDiscreteModel)
labels = get_face_labeling(model)
bgface_to_mask = get_face_mask(labels,"boundary",2)
Γface_to_bgface = findall(bgface_to_mask)
Dc2Dp3model = BoundaryDiscreteModel(Polytope{2},model,Γface_to_bgface)

order  = 0
degree = 1

reffe_rt = ReferenceFE(raviart_thomas,Float64,order)
V  = FESpace(Dc2Dp3model, reffe_rt ; conformity=:HDiv)
U = TrialFESpace(V,u)
reffe = ReferenceFE(lagrangian,Float64,order)
Q = TestFESpace(Dc2Dp3model,reffe,conformity=:L2)
P = TrialFESpace(Q)
uh = FEFunction(V,rand(num_free_dofs(V)))
vh = interpolate_everywhere(uh,V)

Ω = Triangulation(Dc2Dp3model)
dΩ = Measure(Ω,2*order)
e=sqrt(sum(∫((uh-vh)⋅(uh-vh))dΩ))
@test e < 1.0e-12

test_div_v_q_equiv(U,V,P,Q,Ω)


end # module
