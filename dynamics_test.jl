using SatelliteDynamics, MATLAB

function FODE(t,eci_state,u)

    # unpack state
    r_eci = eci_state[1:3]
    v_eci = eci_state[4:6]


    # J2 only
    # a_eci = FODE_J2(r_eci)
    a_eci = -GM_EARTH*r_eci/(norm(r_eci)^3)


    return [v_eci;a_eci+ u]
end

function FODE_J2(r_eci)

    # relavent parameters
    J2 = J2_EARTH
    # J2 = 0.0
    μ = GM_EARTH

    # eci position stuff
    r = norm(r_eci)
    x,y,z = r_eci

    # precompute repeated stuff
    Re_r_sqr = 1.5*J2*(R_EARTH/r)^2
    five_z_sqr = 5*z^2/r^2

    return  (-μ/r^3)*[x*(1 - Re_r_sqr*(five_z_sqr - 1));
                      y*(1 - Re_r_sqr*(five_z_sqr - 1));
                      z*(1 - Re_r_sqr*(five_z_sqr - 3))]

end


function rk4_orbital(f::Function, t_n, x_n, u, h::Real)
    """Runge-Kutta 4th order integration. Epoch for time.

    Args:
        ODE:           f(t,x,u)        | function
        initial time:  t_n             | epoch or float64
        initial cond:  x_n             | vec
        control input: u               | vec
        time step:     h               | scalar

    Returns:
        x_{n+1}
    """

    k1 = h*f(t_n,x_n,u)
    k2 = h*f(t_n+h/2,x_n+k1/2,u)
    k3 = h*f(t_n+h/2,x_n+k2/2,u)
    k4 = h*f(t_n+h,x_n+k3,u)

    x_np1 = (x_n + (1/6)*(k1+2*k2+2*k3 + k4))



    return x_np1
end
# Declare simulation initial Epoch
epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0)

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 500e3, 0.01, 75.0, 45.0, 30.0, 0.0]

# Convert osculating elements to Cartesean state
eci0 = sOSCtoCART(oe0, use_degrees=true)

N = 5*15000
dt = .1
X = fill(zeros(6),N)
X[1] = eci0

for i = 1:(N-1)
    u = [0.1;.2;.3]
    X[i+1] = rk4_orbital(FODE,0,X[i],u,dt)
end

Xm = mat_from_vec(X)


struct KSopt
    μ::Float64
    tscale::Float64
    dscale::Float64
    uscale::Float64
end
# scaling
# dscale = 1e7 # m
# tscale = 200 # s
# uscale = 1000000.0
dscale = 1
tscale = 1
# uscale = 1
# μ = 3.986004415e14 *(tscale^2)/(dscale^3)
μ = GM_EARTH
model = KSopt(μ,tscale, dscale, uscale)
R_earth = R_EARTH
J2 = J2_EARTH
function dyno(t, x, u)

    # unpack state
    p = SVector(x[1],x[2],x[3],x[4])
    p_prime = SVector(x[5],x[6],x[7],x[8])
    h = x[9]

    # scale u
    # u /= model.uscale
    r = x_from_u(p)*dscale                       # ECI, m
    mu = μ

    # position vector in ECI from earth center to spacecraft
    r_i = r[1]                             # ECI, m
    r_j = r[2]                             # ECI, m
    r_k = r[3]                             # ECI, m
    R = norm(r)                          # ECI, m (dot(u,u)=norm(x))

    # a_j2 = FODE_J2(r)
    # a_j2 = zeros(3)
    # pertubation acceleration in ECI
    # P_j2 = (-1.5*(mu*J2*R_earth^2)/(R^4)) * [(r_i/R)*(1-3*r_k^2/R^2);(r_j/R)*(1-3*r_k^2/R^2);(r_k/R)*(3-3*r_k^2/R^2)]
    # P_j2 = P_j2*tscale^2/dscale
    P_j2 = SVector(u[1],u[2],u[3])*tscale^2/dscale
    # Levi-Civita transformation matrix
    L = L_fx(p)

    # append 0 to end of P_j2 acceleration vector
    P_j2 = SVector{4}([P_j2;0]*tscale^2/dscale)
    # P_j2 = @SVector zeros(4)
    L = L_fx(p)
    # P_j2 = SVector(u[1],u[2],u[3],0)/uscale
    p_dp = -(h/2)*p + ((dot(p,p))/(2*dscale))*L'*P_j2

    return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
                   p_dp[1],p_dp[2],p_dp[3],p_dp[4],
                   -2*dot(p_prime,L'*P_j2),
                   dot(p,p),
                   2*p'*p_prime,
                   2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]))

end



function dynamics!(t,x,u)

    # scaling info
    # dscale = p.dscale
    # tscale = p.tscale
    # uscale = p.uscale

    # unpack state
    x = x[:]
    ẋ = zeros(length(x))
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
    P_j2 = (-1.5*(mu*J2*R_earth^2)/(R^4)) * [(r_i/R)*(1-3*r_k^2/R^2);(r_j/R)*(1-3*r_k^2/R^2);(r_k/R)*(3-3*r_k^2/R^2)]
    P_j2 = P_j2*tscale^2/dscale
    # # Levi-Civita transformation matrix
    L = L_fx(p)
    #
    # # append 0 to end of P_j2 acceleration vector
    P_j2 = [P_j2;0]

    # no J2
    # P_j2 = [u[1];u[2];u[3];0]

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

















r_eci = eci0[1:3]
v_eci = eci0[4:6]
R_0 = r_eci/dscale
V_0 = v_eci*tscale/dscale
p_0 = u_from_x(R_0)                          # u
p_prime_0 = uprime_from_xdot(V_0,p_0)     # du/ds
h_0 = μ/norm(R_0) - .5*norm(V_0)^2         # m^2/s^2
t_0 = 0                                    # seconds
x0 = SVector{12}([p_0; p_prime_0; h_0; t_0; norm(R_0); R_0[3]])

dt = 0.0000002
N = 5000
X2 = fill(zeros(12),N)

X2[1] = x0
for i = 1:(N-1)
    u = [0.1;.2;.3]
    X2[i+1] = rk4_orbital(dyno,0,X2[i],u,dt)
    # X2[i+1] = rk4_orbital(dynamics!,0,X2[i],0,dt)
end

Xm2 = mat_from_vec(X2)

eci2 = [x_from_u(X2[i][1:4]) for i = 1:length(X2)]

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

eci2m = dscale*mat_from_vec2(eci2)










mat"
figure
hold on
plot3($Xm(1,:),$Xm(2,:),$Xm(3,:))
plot3($eci2m(1,:),$eci2m(2,:),$eci2m(3,:),'r*')
hold off
"
