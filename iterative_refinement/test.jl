using LinearAlgebra

function ss()
# n= 10
# A = [10000000*randn(n,n) randn(n,n);
#      randn(n,n) 0.0000001*randn(n,n)]
# # A = [3.9 1.6;6.8 2.9]
# # A[2,:] = (ones(2) + 1e-10*randn(2)) .* A[1,:]
# # A = float([64       -2016       20160      -92400      221760     -288288      192192      -51480;
# # -2016       84672     -952560     4656960   -11642400    15567552   -10594584     2882880 ;
# # 20160     -952560    11430720   -58212000   149688000  -204324120   141261120   -38918880 ;
# # -92400     4656960   -58212000   304920000  -800415000  1109908800  -776936160   216216000 ;
# # 221760   -11642400   149688000  -800415000  2134440000 -2996753760  2118916800  -594594000 ;
# # -288288    15567552  -204324120  1109908800 -2996753760  4249941696 -3030051024   856215360 ;
# # 192192   -10594584   141261120  -776936160  2118916800 -3030051024  2175421248  -618377760 ;
# # -51480     2882880   -38918880   216216000  -594594000   856215360  -618377760   176679360 ])
#
# # b = zeros(8)
# # b[3] = 1
# b = randn(n*2)
# x = A\b
#
# @show A*x - b
#
# @show cond(A)
#
# println("Naive Guess")
#
# # Ap = A + 1e-5*I
# @show norm(A*x - b)

n = 6
m = 4
P = randn(n,n);P = P'*P;
Al = randn(m,n)
Au = randn(m,n)

K = [P Al' Au';Al zeros(m,m) zeros(m,m); Au zeros(m,m) zeros(m,m)]

@show cond(K)

δ = 1e-6
K2 = K + Diagonal([δ*ones(n);δ*ones(2*m)])

g = randn(n +2*m)
print("naive")
x_n = K\g
@show norm(K*x_n - g)

t = K2\g

for i = 1:4
    Δt = K2\(g-K*t)
    # @show norm(Δt)
    t = t + Δt
    @show norm(K*t - g)
end

# @show size(K)



# Ap = A + 1e-6*I
# F = factorize(Ap)
#
# x = F\b
# @show norm(A*x - b)
# for i = 1:5
#     Δx = F\(b - A*x)
#     x -= Δx
#     @show norm(A*x - b)
# end

#
# Af = factorize(A)
# # r = zeros(2*n)
# for i = 1:5
#     r = b-A*x
#     # r .= A*x - b
#     # x -= A\r
#     d = Af\r
#     x = x +d
#
#     @show norm(A*x -b)
# end
#
# @show norm(A*x -b)

end

ss()


# regularize it, solve
# then take the residual and solve with non-regularized matrix
#
# use backsolves with regularized matrix
#
# residuals with regular matrix
#
# regularize kkt system
# use backslack
# iterative refinement
#
# try to figure out what the ricatti is
#
# OSQP defaults to 1e-5 for altro
