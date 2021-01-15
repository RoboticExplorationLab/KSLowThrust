using SatelliteDynamics, MATLAB, StaticArrays, Attitude

# const DE = OrdinaryDiffEq

# let's try the 2d case

#                       a e i Ω ω M
oe = [35000*1e3, 0.7, 0, 0,0,0]
rv_eci = sOSCtoCART(oe,use_degrees = true)

# Declare simulation initial Epoch
epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0)
T = orbit_period(oe[1])
# Create an EarthInertialState orbit propagagator
orb  = EarthInertialState(epc0, rv_eci, dt=10.0,
            mass=1.0, n_grav=0, m_grav=0,
            drag=false, srp=false,
            moon=false, sun=false,
            relativity=false
)

# Propagate the orbit
t, epc, eci_hist = sim!(orb, epc0 + T)

mat"
figure
hold on
plot($eci_hist(1,:),$eci_hist(2,:))
hold off
axis equal
"


# function star_conj(u)
#     return [u[1:3];-u[4]]
# end
# ks sim
include("/Users/kevintracy/devel/Low-Thrust-TrajOpt/ks_functions.jl")
u0 = u_from_x(rv_eci[1:3])
up0 = uprime_from_xdot(rv_eci[4:6],u0)

h0 = GM_EARTH/(norm(rv_eci[1:3])) - norm(rv_eci[4:6])^2/2

# function ksdyn(state,p,t)
#     u  = SA[state[1],state[2],state[3],state[4]]
#     up = SA[state[5],state[6],state[7],state[8]]
#     upp = -(h0/2)*u
#     return SA[state[5],state[6],state[7],state[8],upp[1],upp[2],upp[3],upp[4]]
# end

# state0 = SVector{8}([u0;up0])
# tspan = (0.0,18.65e-4)
# prob = ODEProblem(ksdyn,state0,tspan)
# sol = solve(prob,DE.RK4(),dt = 1e-6, adaptive = false)


# eci_hist2 = mat_from_vec([x_from_u(sol.u[i][1:4]) for i = 1:length(sol.u)])
#
# mat"
# figure
# hold on
# plot($eci_hist2(1,:),$eci_hist2(2,:))
# hold off
# axis equal
# "

# function true_anom_from_rv(r,v)
#     h = cross(r,v)
#     e = (v × h)/GM_EARTH - r/norm(r)
#
#     ν = acos(dot(e,r)/(norm(e)*norm(r)))
#     if dot(r,v)<0
#         ν = 2*pi - ν
#     end
#     return ν
# end
# tf = sol.t[end]
#
# oe_hist = cfill(6,20)
# anom_hist = cfill(3,20)
# s_hist = zeros(20)
# for i = 1:20
#     t = (i-1)*tf/19
#     s_hist[i] = t
#     state = sol(t)
#     u = state[1:4]
#     up = state[5:8]
#     r = x_from_u(u)
#     v = xdot_from_u(u,up)
#     oe = sCARTtoOSC([r;v])
#     M = oe[6]
#     e = oe[2]
#     E = anomaly_mean_to_eccentric(M,e,use_degrees = false)
#     ν = true_anom_from_rv(r,v)
#     anom_hist[i] = [M,E,ν]
#     oe_hist[i] = oe
# end
#
# a = mat_from_vec(anom_hist)
#
# mat"
# figure
# hold on
# plot($a')
# hold off
# "
#
#
# u69 = sol.u[end-34][1:4]
#
# x69 = x_from_u(u69)
#
# # Q = dcm_from_q(u69)
#
# nn = norm(u69)
#
# # Q = dcm_from_q(u69/nn)
# Q = dcm_from_q([u69[4];u69[1:3]]/nn)
#
#
# function dcm(q)
#     v = q[1:3]
#     s = q[4]
#     # s = q[1]
#     # v = q[2:4]
#     return I + 2*hat(v)*(s*I + hat(v))
# end


# u =
# SVector( u[1]^2 - u[2]^2 - u[3]^2 + u[4]^2,
#                  2*(u[1]*u[2] - u[3]*u[4]),
#                  2*(u[1]*u[3] + u[2]*u[4]))
# function dcm_from_q(q::Vec)::Mat
#     """DCM from quaternion, hamilton product, scalar last"""
#     q4,q1,q2,q3 = normalize(q)
#
#     # DCM
#     Q = @SArray [(2*q1^2+2*q4^2-1)   2*(q1*q2 - q3*q4)   2*(q1*q3 + q2*q4);
#           2*(q1*q2 + q3*q4)  (2*q2^2+2*q4^2-1)   2*(q2*q3 - q1*q4);
#           2*(q1*q3 - q2*q4)   2*(q2*q3 + q1*q4)  (2*q3^2+2*q4^2-1)]
# end
function imag_to_vec(x)
    return [real(x);imag(x)]
end
function vec_to_imag(r)
    return r[1] + r[2]im
end
function mysqr(x)
    return (x + norm(x))/(sqrt(2*(norm(x) + real(x))))
end
function u_from_x2(x)
    return sqrt(x)
end
function x_from_u2(u)
    return u^2
end
r = rv_eci[1:2]
v = rv_eci[4:5]
x = r[1] + r[2]im
xd = v[1] + v[2]im
r = norm(x)
u = sqrt(x)

up = r*xd/(2*u)

h = GM_EARTH/r - .5*norm(xd)^2


A = Array([zeros(2,2) I(2);(-h/2)*I(2) zeros(2,2)])

# statee = cfill(2,100)
# statee[1] = [imag_to_vec(u);imag_to_vec(up)]

ds = 1e-8

statee = [ exp(A*((i-1)*ds))*[imag_to_vec(u);imag_to_vec(up)] for i = 1:100]

xhist = [imag_to_vec( x_from_u2((vec_to_imag(statee[i][1:2])))) for i = 1:length(statee)]

xm = mat_from_vec(xhist)
mat"
figure
hold on
plot($xm(1,:),$xm(2,:),'o')
hold off
"

s = mat_from_vec(statee)
mat"
figure
hold on
plot($s')
hold off
"
