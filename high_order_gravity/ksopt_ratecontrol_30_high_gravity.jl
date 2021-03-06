using TrajectoryOptimization
using Altro
using RobotDynamics
using Attitude
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
function RobotDynamics.dynamics(model::KSopt, x, u_dot)

    # unpack state
    p = SVector(x[1],x[2],x[3],x[4])
    p_prime = SVector(x[5],x[6],x[7],x[8])
    h = x[9]
    #r = x[10]
    #z = x[11]
    u = SVector(x[12],x[13],x[14])
    t = x[15]

    r_eci = x_from_u(p)*model.dscale
    # this is iauEra00 but I wrote it in a diffable way
    # R_eci_ecef = myRz(1.0047517553493037 + 6.300387486754831*(mjd(epc0) + t*model.tscale/(24*3600)))
    R_eci_ecef = myRz(1.0047517553493037 + 6.300387486754831*(58484.5 + t*model.tscale/(24*3600)))
    a_eci = SVector{3}(spherical_harmonic_gravity_wrapper(r_eci,R_eci_ecef))

    u_eci = SVector(u[1],u[2],u[3])/model.uscale
    # Levi-Civita transformation matrix
    L = L_fx(p)

    # append 0 to end of P_j2 acceleration vector
    P_j2 = SVector{4}([u_eci + a_eci;0]*model.tscale^2/model.dscale)

    p_dp = -(h/2)*p + ((dot(p,p))/(2))*L'*P_j2

    return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],  # p'
                   p_dp[1],p_dp[2],p_dp[3],p_dp[4],              # p''
                   -2*dot(p_prime,L'*P_j2),                      # h'
                   2*p'*p_prime,                                 # r'
                   2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]), # z'
                   u_dot[1],u_dot[2],u_dot[3],
                   dot(p,p))                   # u'

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
μ = GM_EARTH *(tscale^2)/(dscale^3)
model = KSopt(μ,tscale, dscale, uscale)
Base.size(::KSopt) = 15,3
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
x0 = SVector{n}([p_0; p_prime_0; h_0; norm(R_0); R_0[3];zeros(3);t_0])

# time stuff
# epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0)
N = 1501
dt = 3e-2
tf = dt*(N-1)

scaling_num = 10
r_Desired = 42e6/dscale#4.2#6.77814
R_scale = .01
Q = Diagonal(SVector{n}([0*ones(9);scaling_num;100;R_scale*ones(3);0]))
R = 0.001*Diagonal(SVector(R_scale,R_scale,R_scale))
q = SVector{n}([zeros(9);-r_Desired*scaling_num;0;zeros(3);0])
costfun = QuadraticCost(Q,R,q=q)
costfun_term = QuadraticCost(Q,R,q=q)
obj = Objective(costfun, costfun_term, N)


# constraints
cons = ConstraintList(n,m,N)
u_bnd = 10.0
maxThrust = NormConstraint(n, m, u_bnd, TO.SecondOrderCone(), 12:14)
TO.add_constraint!(cons, maxThrust, 1:N)
# add_constraint!(cons, BoundConstraint(n,m, u_min=-1000, u_max=1000), 1:N-1)

xf = zeros(n)
xf[11] = 4.2
xf = SVector{n}(xf)
prob = Problem(model, obj,xf, tf, x0=x0, constraints=cons,integration=RD.RK4)
# prob = Problem(model, obj,xf*1000, tf, x0=x0,integration=RD.RK4)

using JLD2
@load "/Users/kevintracy/devel/Low-Thrust-TrajOpt/rate_control_transfers/30_v3_day_transfer.jld2" X U t_hist
U_old = copy(U)
initial_controls!(prob,U_old)
opts = SolverOptions(verbose=2,
                     show_summary = true,
                     iterations = 20000,
                     projected_newton = false,
                     constraint_tolerance = 1e-3,
                     penalty_scaling = 2.0,
                     penalty_initial = 5.0)
solver = ALTROSolver(prob,opts)

solve!(solver)

# get time stuff
function get_time_transfer(solver,dt,tscale)
    """get the time in days for a run. This allows us to remove t from state"""
    X = states(solver)
    t = 0.0
    t_hist = zeros(length(X))
    for i = 1:(length(X)-1)

        p0 = X[i][1:4]
        R0 = dot(p0,p0)

        p1 = X[i+1][1:4]
        R1 = dot(p1,p1)

        t += (R0 + R1)*dt/2
        t_hist[i+1] = copy(t)
    end
    t_days = t*tscale/(24*3600)
    t_hist = t_hist*tscale/(24*3600)
    return t_days, t_hist
end

t_days,t_hist = get_time_transfer(solver,dt,tscale)

X = states(solver)
U = controls(solver)

r_eci_hist = dscale*mat_from_vec([x_from_u(X[i][1:4]) for i = 1:length(X)])
v_eci_hist = (dscale/tscale)*mat_from_vec([xdot_from_u(X[i][1:4],X[i][5:8]) for i = 1:length(X)])

mat"
figure
hold on
plot3($r_eci_hist(1,:),$r_eci_hist(2,:),$r_eci_hist(3,:))
%axis equal
hold off
"

Xm = mat_from_vec(X)
U_real = Xm[12:14,:]

mat"
figure
hold on
title('Thrust Vector')
plot($t_hist,$U_real')
ylabel('Thrust Acceleration (0.0001m/s^2)')
xlabel('Time (days)')
hold off
xlim([0,$t_hist(end)])
legend('u_x','u_y','u_z')
%saveas(gcf,'long_traj.png')
"


U_real_vec = [X[i][12:14] for i = 1:length(X)]
Unorm = [norm(U_real_vec[i]) for i = 1:length(U)]

mat"
figure
hold on
plot($Unorm)
hold off
"


println("Time in days:  $t_days")# using TrajectoryOptimization


# now let's get orbital elements
const SD = SatelliteDynamics

oe_hist = [SD.sCARTtoOSC([r_eci_hist[:,i];v_eci_hist[:,i]]) for i = 1:size(r_eci_hist,2)]

a_hist = [oe_hist[i][1] for i = 1:length(oe_hist)]
e_hist = [oe_hist[i][2] for i = 1:length(oe_hist)]
i_hist = [rad2deg(oe_hist[i][3]) for i = 1:length(oe_hist)]

mat"
figure
hold on
subplot(3,1,1)
hold on
plot($t_hist,$a_hist)
subplot(3,1,2)
hold on
plot($t_hist,$e_hist)
subplot(3,1,3)
hold on
plot($t_hist,$i_hist)
hold off
"

# using JLD2
# X = states(solver)
# U = controls(solver)
# @save "rate_control_transfers/30_spherical_harmonic_day_transfer.jld2" X U t_hist
#
