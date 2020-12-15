using TrajectoryOptimization
using Altro
using RobotDynamics
using StaticArrays
using LinearAlgebra
const TO = TrajectoryOptimization
const RD = RobotDynamics

function L_fx(u)

    return @SMatrix [u[1] -u[2] -u[3]  u[4];
                     u[2]  u[1] -u[4] -u[3];
                     u[3]  u[4]  u[1]  u[2];
                     u[4] -u[3]  u[2] -u[1]]

end

struct KSopt <: TO.AbstractModel
    μ::Float64
    tscale::Float64
    dscale::Float64
    uscale::Float64
end
function RobotDynamics.dynamics(model::KSopt, x, u)

    # unpack state
    p = SVector(x[1],x[2],x[3],x[4])
    p_prime = SVector(x[5],x[6],x[7],x[8])
    h = x[9]


    P_j2 = SVector(u[1],u[2],u[3])/model.uscale
    # Levi-Civita transformation matrix
    L = L_fx(p)

    # append 0 to end of P_j2 acceleration vector
    P_j2 = SVector{4}([P_j2;0]*model.tscale^2/model.dscale)

    p_dp = -(h/2)*p + ((dot(p,p))/(2))*L'*P_j2

    return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
                   p_dp[1],p_dp[2],p_dp[3],p_dp[4],
                   -2*dot(p_prime,L'*P_j2),
                   dot(p,p),
                   2*p'*p_prime,
                   2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]))

end
# function RobotDynamics.dynamics(model::KSopt, x, u)
#
#     # unpack state
#     p = SVector(x[1],x[2],x[3],x[4])
#     p_prime = SVector(x[5],x[6],x[7],x[8])
#     h = x[9]
#
#     # scale u
#     # u /= model.uscale
#
#     L = L_fx(p)
#     P_j2 = SVector(u[1],u[2],u[3],0)/model.uscale
#     p_dp = -(h/2)*p + ((dot(p,p))/(2*dscale))*L'*P_j2
#
#     return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
#                    p_dp[1],p_dp[2],p_dp[3],p_dp[4],
#                    -2*dot(p_prime,L'*P_j2),
#                    dot(p,p),
#                    2*p'*p_prime,
#                    2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]))
#
# end

# scaling
dscale = 1e7 # m
tscale = 20000 # s
uscale = 10000.0
μ = 3.986004418e14 *(tscale^2)/(dscale^3)
model = KSopt(μ,tscale, dscale, uscale)
Base.size(::KSopt) = 12,3
n,m = size(model)

# initial conditions
r_eci = @SVector [6.578137000000001e6,0,0]
v_eci = @SVector [0,9111.210131355681,4642.393437617197]
R_0 = r_eci/dscale
V_0 = v_eci*tscale/dscale
p_0 = u_from_x(R_0)                          # u
p_prime_0 = uprime_from_xdot(V_0,p_0)     # du/ds
h_0 = μ/norm(R_0) - .5*norm(V_0)^2         # m^2/s^2
t_0 = 0                                    # seconds
x0 = SVector{12}([p_0; p_prime_0; h_0; t_0; norm(R_0); R_0[3]])

# time stuff
N = 2001
dt = 8e-2
tf = dt*(N-1)

# cost function
# Q = 10000*Diagonal(SVector{n}([zeros(10);1;0]))
# R = 1*Diagonal(@SVector fill(1.0, m))
# r_Desired = 42e6/dscale
# xf = SVector{n}([zeros(10);r_Desired;0])
# Qf = 10*Q
# obj = LQRObjective(Q, R, Qf, xf, N)
scaling_num = 10
r_Desired = 42e6/dscale#4.2#6.77814
R_scale = .1
Q = Diagonal(SVector{n}([1e-12*ones(10);scaling_num;1000]))
R = Diagonal(SVector(R_scale,R_scale,R_scale))
q = SVector{n}([zeros(10);-r_Desired*scaling_num;0])
costfun = QuadraticCost(Q,R,q=q)
costfun_term = QuadraticCost(Q,R,q=q)
obj = Objective(costfun, costfun_term, N)


# constraints
cons = ConstraintList(n,m,N-1)
u_bnd = 5
maxThrust = NormConstraint(n, m, u_bnd, TO.SecondOrderCone(), :control)
TO.add_constraint!(cons, maxThrust, 1:N-1)
# add_constraint!(cons, BoundConstraint(n,m, u_min=-1000, u_max=1000), 1:N-1)

xf = randn(n)
xf[11] = 4.2
xf = SVector{n}(xf)
prob = Problem(model, obj,xf, tf, x0=x0, constraints=cons,integration=RD.RK4)
# prob = Problem(model, obj,xf*1000, tf, x0=x0,integration=RD.RK4)
opts = SolverOptions(verbose=2,show_summary = true,iterations = 20000,
projected_newton = false, constraint_tolerance = 1e-4)
solver = ALTROSolver(prob,opts)

solve!(solver)


X = states(solver)
function mat_from_vec2(a)
    "Turn a vector of vectors into a matrix"


    rows = length(a[1])
    columns = length(a)
    A = zeros(rows,columns)

    for i = 1:columns
        A[:,i] = a[i]
    end

    return A
end
r_eci_hist = mat_from_vec2([x_from_u(X[i]) for i = 1:length(X)])


mat"
figure
hold on
plot3($r_eci_hist(1,:),$r_eci_hist(2,:),$r_eci_hist(3,:))
axis equal
hold off
"

U = mat_from_vec2(controls(solver))

mat"
figure
hold on
plot($U')
hold off
"






# using TrajectoryOptimization
# using Altro
# using RobotDynamics
# using StaticArrays
# using LinearAlgebra
# const TO = TrajectoryOptimization
# const RD = RobotDynamics
#
# function L_fx(u)
#
#     return @SMatrix [u[1] -u[2] -u[3]  u[4];
#                      u[2]  u[1] -u[4] -u[3];
#                      u[3]  u[4]  u[1]  u[2];
#                      u[4] -u[3]  u[2] -u[1]]
#
# end
#
# struct KSopt <: TO.AbstractModel
#     μ::Float64
#     tscale::Float64
#     dscale::Float64
#     uscale::Float64
# end
# function RobotDynamics.dynamics(model::KSopt, x, u)
#
#     # unpack state
#     p = SVector(x[1],x[2],x[3],x[4])
#     p_prime = SVector(x[5],x[6],x[7],x[8])
#     h = x[9]
#
#
#     P_j2 = SVector(u[1],u[2],u[3])/model.uscale
#     # Levi-Civita transformation matrix
#     L = L_fx(p)
#
#     # append 0 to end of P_j2 acceleration vector
#     P_j2 = SVector{4}([P_j2;0]*model.tscale^2/model.dscale)
#
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
# # function RobotDynamics.dynamics(model::KSopt, x, u)
# #
# #     # unpack state
# #     p = SVector(x[1],x[2],x[3],x[4])
# #     p_prime = SVector(x[5],x[6],x[7],x[8])
# #     h = x[9]
# #
# #     # scale u
# #     # u /= model.uscale
# #
# #     L = L_fx(p)
# #     P_j2 = SVector(u[1],u[2],u[3],0)/model.uscale
# #     p_dp = -(h/2)*p + ((dot(p,p))/(2*dscale))*L'*P_j2
# #
# #     return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
# #                    p_dp[1],p_dp[2],p_dp[3],p_dp[4],
# #                    -2*dot(p_prime,L'*P_j2),
# #                    dot(p,p),
# #                    2*p'*p_prime,
# #                    2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]))
# #
# # end
#
# # scaling
# dscale = 1e7 # m
# tscale = 20000 # s
# uscale = 1000.0;
# μ = 3.986004418e14 *(tscale^2)/(dscale^3)
# model = KSopt(μ,tscale, dscale, uscale)
# Base.size(::KSopt) = 12,3
# n,m = size(model)
#
# # initial conditions
# r_eci = @SVector [6.578137000000001e6,0,0]
# v_eci = @SVector [0,9111.210131355681,4642.393437617197]
# R_0 = r_eci/dscale
# V_0 = v_eci*tscale/dscale
# p_0 = u_from_x(R_0)                          # u
# p_prime_0 = uprime_from_xdot(V_0,p_0)     # du/ds
# h_0 = μ/norm(R_0) - .5*norm(V_0)^2         # m^2/s^2
# t_0 = 0                                    # seconds
# x0 = SVector{12}([p_0; p_prime_0; h_0; t_0; norm(R_0); R_0[3]])
#
# # time stuff
# N = 501
# dt = 1.0
# tf = dt*(N-1)
#
# # cost function
# Q = 10000*Diagonal(SVector{n}([zeros(10);1;0]))
# R = 1*Diagonal(@SVector fill(1.0, m))
# r_Desired = 42e6/dscale
# xf = SVector{n}([zeros(10);r_Desired;0])
# Qf = 10*Q
# obj = LQRObjective(Q, R, Qf, xf, N)
# # scaling_num = 10
# # r_Desired = 38e6/dscale#4.2#6.77814
# # R_scale = 1
# # Q = Diagonal(SVector{n}([1e-6*ones(10);scaling_num;100]))
# # R = Diagonal(SVector(R_scale,R_scale,R_scale))
# # q = SVector{n}([zeros(10);-r_Desired*scaling_num;0])
# # costfun = QuadraticCost(Q,R,q=q)
# # costfun_term = QuadraticCost(Q,R,q=q)
# # obj = Objective(costfun, costfun_term, N)
#
#
# # constraints
# cons = ConstraintList(n,m,N-1)
# u_bnd = 5.0
# maxThrust = NormConstraint(n, m, u_bnd, TO.SecondOrderCone(), :control)
# TO.add_constraint!(cons, maxThrust, 1:N-1)
#
#
# # xf = randn(n)
# # xf[11] = 4.2
# # xf = SVector{n}(xf)
# prob = Problem(model, obj,xf, tf, x0=x0, constraints=cons,integration=RD.RK4)
# # prob = Problem(model, obj,xf*1000, tf, x0=x0,integration=RD.RK4)
# opts = SolverOptions(verbose=2,show_summary = true,
# projected_newton = false)
# solver = ALTROSolver(prob,opts)
#
# solve!(solver)
#
# X = states(solver)
# function mat_from_vec2(a)
#     "Turn a vector of vectors into a matrix"
#
#
#     rows = length(a[1])
#     columns = length(a)
#     A = zeros(rows,columns)
#
#     for i = 1:columns
#         A[:,i] = a[i]
#     end
#
#     return A
# end
# r_eci_hist = mat_from_vec2([x_from_u(X[i]) for i = 1:length(X)])
#
#
# mat"
# figure
# hold on
# plot3($r_eci_hist(1,:),$r_eci_hist(2,:),$r_eci_hist(3,:))
# %axis equal
# hold off
# "
#
# U = mat_from_vec2(controls(solver))
#
# mat"
# figure
# hold on
# plot($U')
# hold off
# "
#
#
# # Q and R
# scaling_num = 10
# r_Desired = 42e6/dscale#4.2#6.77814
# # Q = Diagonal(1e-6*I,n)
# Q = Diagonal(SVector{n}([1e-6*ones(10);scaling_num;1]))
# # Q[11,11]=1*scaling_num
#
# # Q[12,12]= 1
# # Q = Diagonal(Q)
# #R = Diagonal(100000000000.0*I/uscale^2,m)
# R = Diagonal(SVector(1,1,1.0))
# # q = zeros(n)
# q = SVector{n}([zeros(10);-r_Desired*scaling_num;0])
# # q[11] = -r_Desired*scaling_num
