using ModelingToolkit, SparseArrays, ForwardDiff, StaticArrays
const FD = ForwardDiff
const MTK = ModelingToolkit
function L_fx(u)

    return @SMatrix [u[1] -u[2] -u[3]  u[4];
                     u[2]  u[1] -u[4] -u[3];
                     u[3]  u[4]  u[1]  u[2];
                     u[4] -u[3]  u[2] -u[1]]

end
function Lt_fx(u)
    return [u[1]  u[2]  u[3]  u[4];
         -u[2]  u[1]  u[4] -u[3];
         -u[3] -u[4]  u[1]  u[2];
          u[4] -u[3]  u[2] -u[1]]
end
function dumb_transpose(X)
    a,b = size(X)
    out = zeros(eltype(X),b,a)
    for i = 1:a
        for j = 1:b
            out[j,i] = X[i,j]
        end
    end
    return out
end


# function L_fx_2(u)
#
#     return  [u[1] -u[2] -u[3]  u[4];
#                      u[2]  u[1] -u[4] -u[3];
#                      u[3]  u[4]  u[1]  u[2];
#                      u[4] -u[3]  u[2] -u[1]]
#
# end
# function dynamics_vec(x_in)
#     """ Here is the state
#             p  ∈ R⁴ [1:4]
#             p′ ∈ R⁴ [5:8]
#             h  ∈ R  [9]
#            |r| ∈ R  [10]
#             r₃ ∈ R  [11]
#             u  ∈ R³ [12:14]
#             u̇  ∈ R³ [15:17]
#             t  ∈ R  [18]
#     """
#     x = x_in[1:14]
#     u_dot = x_in[15:17]
#     t = x_in[18]
#     # u_dot = @views x[15:17]
#     # t = x[18]
#     # unpack state
#     # p = SVector(x[1],x[2],x[3],x[4])
#     # p_prime = SVector(x[5],x[6],x[7],x[8])
#     # h = x[9]
#     p = x[1:4]
#     p_prime = x[5:8]
#     h = x[9]
#
#
#
#     P_j2 = [x[12:14];0]*(params.tscale^2/params.dscale)/params.uscale
#     # Levi-Civita transformation matrix
#     # L = L_fx(p)
#     L = [u[1] -u[2] -u[3]  u[4];
#         u[2]  u[1] -u[4] -u[3];
#         u[3]  u[4]  u[1]  u[2];
#         u[4] -u[3]  u[2] -u[1]]
#
#     # append 0 to end of P_j2 acceleration vector
#     # P_j2 = SVector{4}([P_j2;0]*params.tscale^2/params.dscale)
#
#     p_dp = -(h/2)*p + ((dot(p,p))/(2))*transpose(L)*P_j2
#
#     # return [p_prime; p_dp;-2*dot(p_prime,transpose(L)*P_j2);
#     #          2*dot(p,p_prime);
#     #          2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]);
#     #          u_dot]
#     return [2*(p[1]*p_prime[1];
#     2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4])]
#
#     # SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
#     #                p_dp[1],p_dp[2],p_dp[3],p_dp[4],
#     #                -2*dot(p_prime,transpose(L)*P_j2),
#     #                2*dot(p,p_prime),
#     #                2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]),
#     #                u_dot[1],
#     #                u_dot[2],
#     #                u_dot[3])
#
# end
function dynamics(x,u_dot,t)
    """ Here is the state
            p  ∈ R⁴ [1:4]
            p′ ∈ R⁴ [5:8]
            h  ∈ R  [9]
           |r| ∈ R  [10]
            r₃ ∈ R  [11]
            u  ∈ R³ [12:14]
            u̇  ∈ R³ [15:17]
            t  ∈ R  [18]
    """

    p = x[1:4]
    p_prime = x[5:8]
    h = x[9]



    P_j2 = [x[12:14];0]*(params.tscale^2/params.dscale)/params.uscale
    # Levi-Civita transformation matrix
    Lt = Lt_fx(p)

    p_dp = -(h/2)*p + ((p[1]^2 + p[2]^2 + p[3]^2 + p[4]^2)/(2))*Lt*P_j2

    return [p_prime; p_dp;(-2*[p_prime[1] p_prime[2] p_prime[3] p_prime[4]]*Lt*P_j2);
             2*(p[1]*p_prime[1] + p[2]*p_prime[2] + p[3]*p_prime[3] + p[4]*p_prime[4]);
             2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]);
             u_dot]
end

function rk4(x)

    # x_n = SVector{6}(x_n)
    x_n =x[1:14]
    u = x[15:17]
    t = x[18]
    dt = x[19]

    k1 = dt*dynamics(x_n,u,t)
    k2 = dt*dynamics(x_n+k1/2, u,t+dt/2)
    k3 = dt*dynamics(x_n+k2/2, u,t+dt/2)
    k4 = dt*dynamics(x_n+k3, u,t+dt)


    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

end

@variables x_sym[1:19]
sym_out = rk4(x_sym)

sym_out2 = simplify.(sym_out)
f_expr = build_function(sym_out2,x_sym)
fast_rk4 = eval(f_expr[1])

@btime fast_rk4(randn(19))

MTK.sparsejacobian(sym_out2,x_sym)
# @btime rk4(randn(19))
#
# @variables x_sym[1:19]
#
# sym_out = rk4(x_sym)
#
# f_expr = build_function(sym_out,x_sym)
#
# fast_rk4 = eval(f_expr[1])
#
#
#
# @variables x_sym[1:18]
#
# sym_out = dynamics_vec(x_sym)
# @btime fast_rk4(randn(19))

# function discrete_dynamics_fast(x,u,t,dt)
#     # return rk4(dynamics,x,u,t,dt)
#     return fast_rk4([x;u;t;dt])
# end


# x_test = [randn(14);randn(3);0;2.0]
# J_test = FD.jacobian(fast_rk4,x_test)
# @btime FD.jacobian(rk4,randn(19))

# @btime FD.jacobian(fast_rk4,randn(19))
# jac = MTK.sparsejacobian(sym_out,x_sym)
