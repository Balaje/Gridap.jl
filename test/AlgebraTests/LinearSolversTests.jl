module LinearSolversTests

using Test
using Gridap.Algebra

using LinearAlgebra
using SparseArrays


diff1(M) = [ [1.0 zeros(1,M-1)]; diagm(1=>ones(M-1)) - Matrix(1.0I,M,M) ]

sdiff1(M) = sparse(diff1(M))

function Laplacian(Nx, Ny, Lx, Ly)
   dx = Lx / (Nx+1)
   dy = Ly / (Ny+1)
   Dx = sdiff1(Nx) / dx
   Dy = sdiff1(Ny) / dy
   Ax = Dx' * Dx
   Ay = Dy' * Dy
   return kron(sparse(1.0I,Ny,Ny), Ax) + kron(Ay, sparse(1.0I,Nx,Nx))
end

n = 10
A = Laplacian(n,n,1,1)
x = rand(n^2)
b = A*x

op = AffineOperator(A,b)
test_nonlinear_operator(op,x,zeros(size(x)),jac=A)

using SymRCM: symrcm
ls = LUSolver(reorder=symrcm)
test_linear_solver(ls,A,b,x)

#using UnicodePlots
#display(spy(A))
#p = symrcm(A)
#display(spy(A[p,p]))

ls = LUSolver()
test_linear_solver(ls,A,b,x)

n = 10
A = Laplacian(n,n,1,1)
x = rand(n^2)
b = A*x

ls = BackslashSolver()
test_linear_solver(ls,A,b,x)

nls = NLSolver(show_trace=false,method=:newton)
x0 = zeros(length(x))
op = AffineOperator(A,b)
@test jacobian(op,x) === op.matrix

cache = solve!(x0,nls,op)
test_nonlinear_solver(nls,op,x0,x)

x1 = copy(x0)
cache = solve!(x1,nls,op,cache)
@test all(x1 .≈ x0)

op = AffineOperator(A,2*b)
solve!(x1,nls,op,cache)
@test all(x1 .≈ 2*x0)

op = AffineOperator(A,b)
x0 = copy(x)
cache = solve!(x0,ls,op,nothing,true)
@test all(x .≈ x0)
cache = solve!(x0,ls,op,nothing)
@test all(x .≈ x0)
cache = solve!(x0,ls,op)
@test all(x .≈ x0)

x0 = copy(x)
cache = solve!(x0,ls,op,cache)
@test all(x .≈ x0)
cache = solve!(x0,ls,op,cache,false)
@test all(x .≈ x0)

op2 = AffineOperator(A,2*b)
x0 = copy(x)
cache = solve!(x0,ls,op2,cache)
@test all(2*x .≈ x0)
cache = solve!(x0,ls,op2,cache,false)
@test all(2*x .≈ x0)
cache = solve!(x0,ls,op2,nothing,false)
@test all(2*x .≈ x0)

end # module
