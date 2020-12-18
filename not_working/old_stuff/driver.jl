using LinearAlgebra, TrajectoryOptimization

# import KS transform functions
include("./KS_transform_functions.jl")

# initial condtions from orbit
Re =   6378137.0; #m

Rp = Re + 200*1000;
Ra = Re + 35000*1000;
i = deg2rad(27);
a = (Rp+Ra)/2;
e = (Ra-Rp)/(Ra+Rp);
OMEGA = 0;
omega = 0;
nu = 0;


r_eci, v_eci = OE2ECI(a,e,i,OMEGA,omega,nu)


# dynamics
function dynamics!(ẋ,x,u,p)

    # scaling info
    dscale = p.dscale
    tscale = p.tscale
    uscale = p.uscale

    # unpack state
    x = x[:]
    p = x[1:4]
    p_prime = x[5:8]
    h = x[9]
    #tee = state[10]
    #r = state[11]
    #k = state[12]

    # rescale u
    u = (u/uscale)

    #J2 stuff
    R_earth  = 6378137.0                   # m
    mu = 3.986004418e14                    # m^3 s^{}-2}
    X = x_from_u(p)*dscale                       # ECI, m
    J2 =  1.081874e-3

    # position vector in ECI from earth center to spacecraft
    r_i = X[1]                             # ECI, m
    r_j = X[2]                             # ECI, m
    r_k = X[3]                             # ECI, m
    R = norm(X)                          # ECI, m (dot(u,u)=norm(x))

    # pertubation acceleration in ECI
    # P_j2 = (-1.5*(mu*J2*R_earth^2)/(R^4)) * [(r_i/R)*(1-3*r_k^2/R^2);(r_j/R)*(1-3*r_k^2/R^2);(r_k/R)*(3-3*r_k^2/R^2)]
    # P_j2 = P_j2*tscale^2/dscale
    # # Levi-Civita transformation matrix
    L = L_fx(p)
    #
    # # append 0 to end of P_j2 acceleration vector
    # P_j2 = [P_j2;0]

    # no J2
    P_j2 = [u[1];u[2];u[3];0]

    # ODE's
    p_dp = -(h/2)*p + (R/(2*dscale))*L'*P_j2
    #p_p = p_prime
    tee_dot = dot(p,p)
    h_dot = -2*dot(p_prime,L'*P_j2)

    # index back into state_prime
    ẋ[1:4] = p_prime
    ẋ[5:8] = p_dp
    ẋ[9] = h_dot
    ẋ[10] = tee_dot
    ẋ[11] = 2*p'*p_prime
    ẋ[12] = 2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4])
    return ẋ

end

# states and control parameters
n, m = 12,3

# initial conditions

dscale = 1e7 # m
tscale = 200 # s
uscale = 1000000.0;
μ = 3.986004418e14 *(tscale^2)/(dscale^3)

# GTO
R_0 = r_eci/dscale
V_0 = v_eci*tscale/dscale


# R_0 = [0.0;6.77814e6/dscale;0]
# vcirc = .95*sqrt(μ/norm(R_0))
# a = .97
# b = sqrt(1-a^2)
# V_0 = [-a*vcirc;0;b*vcirc]


p_0 = u_from_x(R_0)                          # u
p_prime_0 = uprime_from_xdot(V_0,p_0)     # du/ds
h_0 = μ/norm(R_0) - .5*norm(V_0)^2         # m^2/s^2
t_0 = 0                                    # seconds
x0 = vec([p_0' p_prime_0' h_0 t_0 norm(R_0) R_0[3]])

# knot points
N = 1001
#dt_s = 60 # this is the real timestep  in seconds
dt = 10/2#dt_s /tscale # this is adjusted time step

# discretize model

p = (μ = μ, a= 1,tscale = tscale,dscale=dscale,uscale=uscale)
model = Model(dynamics!, n, m,p)
model_d = rk4(model)

# Q and R
scaling_num = 10
r_Desired = 42e6/dscale#4.2#6.77814
Q = Diagonal(1e-6*I,n)
Q[11,11]=1*scaling_num

Q[12,12]= 1
Q = Diagonal(Q)
#R = Diagonal(100000000000.0*I/uscale^2,m)
R = Diagonal(1.0*I,m)
q = zeros(n)
q[11] = -r_Desired*scaling_num

# Costfun and objective
#costfun = QuadraticCost(1e-6Q,1e-6R,q=1e-6q)
costfun = QuadraticCost(Q,R,q=q)
terminalcostfun = QuadraticCost(Q,R,q=q) # removed a q=q
#obj = Objective(costfun,terminalcostfun,N)
obj = Objective(costfun,N)

# constraints
max_accel_real = 1.8*2.6667e-04 #m/s^2
td_scaled_max_accel = max_accel_real*tscale^2/dscale
max_thrust = uscale*(td_scaled_max_accel)# N
c(v,x,u) = v[1] = (u[1]^2 + u[2]^2 + u[3]^2 - max_thrust^2)
p1 = 1
con = Constraint{Inequality}(c, n, m, p1, :mycon1)

# apply max_thrust constraint at every time step
constraints = Constraints(N) # define constraints at each time step
for k = 1:N
    constraints[k] += con
end

# construct problem
# prob = Problem(model_d, obj, constraints=constraints, x0=x0[:], N=N, dt=dt) # removed xf = xf
prob = Problem(model_d, obj, x0=x0[:], N=N, dt=dt)

# initialize controls
U0 = zeros(m,N-1)

# read in saved guess
# U_saved = CSV.read("Ueasydouble6.csv") |> Matrix
# u_1_vec = vec(U_saved[1,:])
# u_2_vec = vec(U_saved[2,:])
# u_3_vec = vec(U_saved[3,:])
# s = 0:10:2999*10
#
# spl_u1 = Spline1D(s,u_1_vec; k = 3)
# spl_u2 = Spline1D(s,u_2_vec; k = 3)
# spl_u3 = Spline1D(s,u_3_vec; k = 3)
#
# s_d = 0:5:5998*5
#
# u1_sd = spl_u1(s_d)
# u2_sd = spl_u2(s_d)
# u3_sd = spl_u3(s_d)
#
# U_saved = [u1_sd';u2_sd';u3_sd']


#guess_type = 1 # random
# guess_type = 2 # saved guess
# for k = 1:N-1
#
#     if guess_type == 1
#     # random guess
#         rand_thrust = randn(3)
#         U0[:,k] .= max_thrust*rand_thrust/norm(rand_thrust)
#     else
#     # previous guess
#     #U0[:,k] .= prob_opt.U[k]
#
#     # if k > size(U_saved)[2]
#     #     rand_thrust = randn(3)
#     #     U0[:,k] .= max_thrust*rand_thrust/norm(rand_thrust)
#     # else
#         # saved guess
#         if k > size(U_saved,2)
#             U0[:,k] .= [0.0;0.0;0.0]
#         else
#         U0[:,k] .= U_saved[:,k]
#         end
#     #end
#     end
# end
initial_controls!(prob,U0)
probsafe = prob
# solve
#opts = iLQRSolverOptions{Float64}(verbose = true)
#opts = ALTROSolverOptions{Float64}(verbose=true)
#opts.opts_al.opts_uncon.verbose = true
ilqropts = iLQRSolverOptions{Float64}(verbose = true, square_root=false)
opts = AugmentedLagrangianSolverOptions{Float64}(verbose=true,
    constraint_tolerance=1.0e-1,opts_uncon=ilqropts,penalty_max=1e12,dual_max=1e12)

prob_opt, solver = solve(prob, ilqropts)
println("solved")


plotting_x_opt = zeros(3,N)
state_mat = zeros(length(prob_opt.X[1]),N)
u_norm = zeros(N)
for i = 1:N
    plotting_x_opt[:,i]=x_from_u(prob_opt.X[i][1:4])
    state_mat[:,i] = prob_opt.X[i]
    if i != N
    u_norm[i] = norm(prob_opt.U[i])
end
end

# plot norm of thruster effort
#plot(1:length(u_norm),u_norm)

pyplot()
plot(plotting_x_opt[1,:],plotting_x_opt[2,:],plotting_x_opt[3,:],xlabel = "ECI I (10,000 km)",ylabel = "ECI J (10,000 km)",zlabel = "ECI K (10,000 km)",title = "Orbital Transfer Trajectory Optimization",label = "Spacecraft")

display(plot(plotting_x_opt[1,:],plotting_x_opt[2,:],plotting_x_opt[3,:],xlabel = "ECI I (10,000 km)",ylabel = "ECI J (10,000 km)",zlabel = "ECI K (10,000 km)",title = "Orbital Transfer Trajectory Optimization",label = "Spacecraft",size = (825, 600)))

plot(1:length(u_norm),u_norm)
save_or_nah = input()
if save_or_nah == nothing
    save_or_nah = input()
end

if save_or_nah == "yes"

    U_mat = zeros(3,N-1)
    for k = 1:N-1

        U_mat[:,k] = prob_opt.U[k]

    end
    println("Got a yes")
    u_write = DataFrame(U_mat)
    state_write = DataFrame(state_mat)
    println("Converted Data")
    CSV.write("Ueasydouble6.csv", u_write)
    CSV.write("state_easydouble6.csv", state_write)
    println("Controls Saved")
else
    println("Controls Aborted")
end
