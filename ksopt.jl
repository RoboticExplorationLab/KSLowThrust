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

    # scale u
    u /= model.uscale

    p_dp = -(h/2)*p + (R/(2*dscale))*L_fx(p)'*SVector(u[1],u[2],u[3],0)

    return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
                   p_dp[1],p_dp[2],p_dp[3],p_dp[4],
                   -2*dot(p_prime,L'*P_j2),
                   dot(p,p),
                   2*p'*p_prime,
                   2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]))

end

# scaling
dscale = 1e7 # m
tscale = 200 # s
uscale = 1000000.0;
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
N = 1001
dt = 2.0
tf = dt*(N-1)

# cost function
Q = Diagonal(SVector{n}([zeros(10);1;1.0]))
R = Diagonal(@SVector fill(1.0, m))
r_Desired = 42e6/dscale
xf = SVector{n}([zeros(10);r_Desired;0])
Qf = 10*Q
obj = LQRObjective(Q, R, Qf, xf, N)

# constraints
cons = ConstraintList(n,m,N-1)
u_bnd = 1000
maxThrust = NormConstraint(n, m, u_bnd, TO.SecondOrderCone(), :control)
TO.add_constraint!(cons, maxThrust, 1:N-1)

prob = Problem(model, obj, xf, tf, x0=x0, constraints=cons,integration=RD.RK4)
solver = ALTROSolver(prob)
