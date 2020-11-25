
function L_fx(u)

    return @SMatrix [u[1] -u[2] -u[3]  u[4];
                     u[2]  u[1] -u[4] -u[3];
                     u[3]  u[4]  u[1]  u[2];
                     u[4] -u[3]  u[2] -u[1]]

end
function dyno(t, x, u)

    # unpack state
    p = SVector(x[1],x[2],x[3],x[4])
    p_prime = SVector(x[5],x[6],x[7],x[8])
    h = x[9]

    # scale u
    # u /= model.uscale
    # r = x_from_u(p)*dscale                       # ECI, m
    # mu = μ



    P_j2 = SVector(u[1],u[2],u[3])
    # Levi-Civita transformation matrix
    L = L_fx(p)

    # append 0 to end of P_j2 acceleration vector
    P_j2 = SVector{4}([P_j2;0]*tscale^2/dscale)

    p_dp = -(h/2)*p + ((dot(p,p))/(2))*L'*P_j2

    return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
                   p_dp[1],p_dp[2],p_dp[3],p_dp[4],
                   -2*dot(p_prime,L'*P_j2),
                   dot(p,p),
                   2*p'*p_prime,
                   2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]))

end
# function dyno(t, x, u)
#
#     # unpack state
#     p = SVector(x[1],x[2],x[3],x[4])
#     p_prime = SVector(x[5],x[6],x[7],x[8])
#     h = x[9]
#
#     # scale u
#     # u /= model.uscale
#     # r = x_from_u(p)*dscale                       # ECI, m
#     # mu = μ
#
#     # position vector in ECI from earth center to spacecraft
#
#     P_j2 = SVector(u[1],u[2],u[3])
#     # Levi-Civita transformation matrix
#     L = L_fx(p)
#
#     # append 0 to end of P_j2 acceleration vector
#     P_j2 = SVector{4}([P_j2;0]*tscale^2/dscale)
#     # P_j2 = SVector{4}([P_j2;0]*0.005)
#     # P_j2 = @SVector zeros(4)
#     L = L_fx(p)
#     # P_j2 = SVector(u[1],u[2],u[3],0)/uscale
#     # @infiltrate
#     p_dp = -(h/2)*p + ((dot(p,p))/(2))*L'*P_j2
#
#     return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
#                    p_dp[1],p_dp[2],p_dp[3],p_dp[4],
#                    -2*dot(p_prime,L'*P_j2),
#                    dot(p,p),
#                    2*p'*p_prime,
#                    2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]))
#
# end

dscale = 1e6
tscale = 200
x = float([i for i = 1:12])

u = [1;2;5]

xdot = dyno(0, x, u)
