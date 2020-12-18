using ModelingToolkit, StaticArrays, SparseArrays
#
# # Define the constants for the PDE
# const α₂ = 1.0
# const α₃ = 1.0
# const β₁ = 1.0
# const β₂ = 1.0
# const β₃ = 1.0
# const r₁ = 1.0
# const r₂ = 1.0
# const _DD = 100.0
# const γ₁ = 0.1
# const γ₂ = 0.1
# const γ₃ = 0.1
# const N = 32
# const X = reshape([i for i in 1:N for j in 1:N],N,N)
# const Y = reshape([j for i in 1:N for j in 1:N],N,N)
# const α₁ = 1.0.*(X.>=4*N/5)
#
# const Mx = Array(Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1]))
# const My = copy(Mx)
# Mx[2,1] = 2.0
# Mx[end-1,end] = 2.0
# My[1,2] = 2.0
# My[end,end-1] = 2.0
#
# # Define the discretized PDE as an ODE function
# function f(u,p,t)
#     A = u[:,:,1]
#     B = u[:,:,2]
#     C = u[:,:,3]
#     MyA = My*A
#     AMx = A*Mx
#     DA = @. _DD*(MyA + AMx)
#     dA = @. DA + α₁ - β₁*A - r₁*A*B + r₂*C
#     dB = @. α₂ - β₂*B - r₁*A*B + r₂*C
#     dC = @. α₃ - β₃*C + r₁*A*B - r₂*C
#     cat(dA,dB,dC,dims=3)
# end
#
# @variables u[1:N,1:N,1:3]
# du = simplify.(f(u,nothing,0.0))
#
# fastf = eval(ModelingToolkit.build_function(du,u,
#             parallel=ModelingToolkit.MultithreadedForm())[2])
#
# # uc = randn(N,N,3)
# jac = ModelingToolkit.sparsejacobian(vec(du),vec(u))
# fjac = eval(ModelingToolkit.build_function(jac,u,
#             parallel=ModelingToolkit.MultithreadedForm())[2])
#
#
# using OrdinaryDiffEq
# u0 = zeros(N,N,3)
# MyA = zeros(N,N);
# AMx = zeros(N,N);
# DA = zeros(N,N);
# # prob = ODEProblem(f!,u0,(0.0,10.0))
# fastprob = ODEProblem(ODEFunction((du,u,p,t)->fastf(du,u),
#                                    jac = (du,u,p,t) -> fjac(du,u),
#                                    jac_prototype = similar(jac,Float64)),
#                                    u0,(0.0,10.0))
#
# using BenchmarkTools
# @btime solve(fastprob, TRBDF2())



# general test
function dyn(x)
    # J = Diagonal([1;2;3])
    ω = x[1:3]
    return inv(J)*(-ω×(J*ω))
end

J = Diagonal([1;2;3])


@variables x_sym[1:3]

sym_out = dyn(x_sym)

f_expr = build_function(sym_out,x_sym)

myf = eval(f_expr[1])


@btime dyn(randn(3))

@btime myf(randn(3))

f_expr2 = build_function(sym_out,x_sym,parallel=ModelingToolkit.MultithreadedForm())

myf2 = eval(f_expr2[1])

@btime myf2(randn(3))
