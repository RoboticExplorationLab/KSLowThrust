# There are two derrivatives of interest in the following functions, variables
# differentiated with respect to time (t) are expressed as dot, and variables
# differentiated with respect to ficticious time (s) are expressed as prime
#
# For a variable V:
#   dV/dt = V_dot
#   dV/ds = V_prime
using StaticArrays
using LinearAlgebra
function L_fx(u)

    return @SMatrix [u[1] -u[2] -u[3]  u[4];
                     u[2]  u[1] -u[4] -u[3];
                     u[3]  u[4]  u[1]  u[2];
                     u[4] -u[3]  u[2] -u[1]]

end

function uprime_from_xdot(xdot,u)
    # given a u and xdot, return u prime
    # xdot = float(xdot)
    # u = float(u)

    u1 = u[1]
    u2 = u[2]
    u3 = u[3]
    u4 = u[4]

    xd1 = xdot[1]
    xd2 = xdot[2]
    xd3 = xdot[3]

    up1 = (1/2)*( u1*xd1  +  u2*xd2  +  u3*xd3  )
    up2 = (1/2)*(-u2*xd1  +  u1*xd2  +  u4*xd3  )
    up3 = (1/2)*(-u3*xd1  -  u4*xd2  +  u1*xd3  )
    up4 = (1/2)*( u4*xd1  -  u3*xd2  +  u2*xd3  )

    uprime = [up1;up2;up3;up4]

    return uprime

    # return SVector((1/2)*( u[1]*xd[1]  +  u[2]*xd[2]  +  u[3]*xd3  ))

end


function u_from_x(x)
    # get KS u from x, u1 is always 1
    x = float(x)
    r = norm(x)

    u1 = .1
    u4 = sqrt(.5*(r+x[1]) - u1^2)
    u2 = (x[2]*u1 + x[3]*u4)/(r+x[1])
    u3 = (x[3]*u1 - x[2]*u4)/(r+x[1])

    u = [u1;u2;u3;u4]
    return u
end


function x_from_u(u)
    # gets x from u
    # u = float(u)

    return SVector( u[1]^2 - u[2]^2 - u[3]^2 + u[4]^2,
                     2*(u[1]*u[2] - u[3]*u[4]),
                     2*(u[1]*u[3] + u[2]*u[4]))

end

function xdot_from_u(u,u_prime)
    up = u_prime;
    u1 = u[1];u2 = u[2];u3 = u[3];u4 = u[4]
    up1 = up[1];up2 = up[2];up3 = up[3];up4 = up[4]
    r = dot(u,u)
    x1d = (2/r)*(u1*up1 - u2*up2 - u3*up3 + u4*up4)
    x2d = (2/r)*(u2*up1 + u1*up2 - u4*up3 - u3*up4)
    x3d = (2/r)*(u3*up1 + u4*up2 + u1*up3 + u2*up4)

    return [x1d;x2d;x3d]
end

# function KS_ODE_J2!(state_prime,state,p,t)
#     # ODE for orbital motion that is with respect to ficticious time (s)
#     # state_prime is d(state)/ds
#     # state is ordered [u1 u2 u3 u4 u1' u2' u3' u4' h t]
#     # time in seconds since epoch will be denoted as 'tee' to avoid confusion
#     # with the ODE time variable
#
#     # scaling info
#     dscale = 1e6
#     tscale = 200
#
#     # unpack state
#     state = state[:]
#     u = state[1:4]
#     u_prime = state[5:8]
#     h = state[9]
#     #tee = state[10]
#     #r = state[11]
#
#     #J2 stuff
#     R_earth  = 6378137.0                   # m
#     mu = 3.986004418e14                    # m^3 s^{}-2}
#     x = x_from_u(u)*dscale                       # ECI, m
#     J2 =  1.081874e-3
#
#     # position vector in ECI from earth center to spacecraft
#     r_i = x[1]                             # ECI, m
#     r_j = x[2]                             # ECI, m
#     r_k = x[3]                             # ECI, m
#     R = norm(x)                          # ECI, m (dot(u,u)=norm(x))
#
#     # pertubation acceleration in ECI
#     P_j2 = (-1.5*(mu*J2*R_earth^2)/(R^4)) * [(r_i/R)*(1-3*r_k^2/R^2);(r_j/R)*(1-3*r_k^2/R^2);(r_k/R)*(3-3*r_k^2/R^2)]
#     P_j2 = P_j2*tscale^2/dscale
#     # Levi-Civita transformation matrix
#     L = L_fx(u)
#
#     # append 0 to end of P_j2 acceleration vector
#     P_j2 = [P_j2;0]
#
#     # no J2
#     #P_j2 = [0;0;0;0]
#
#     # ODE's
#     u_dp = -(h/2)*u + (R/(2*dscale))*L'*P_j2
#     #u_p = u_prime
#     tee_dot = dot(u,u)
#     h_dot = -2*dot(u_prime,L'*P_j2)
#
#     # index back into state_prime
#     state_prime[1:4] = u_prime
#     state_prime[5:8] = u_dp
#     state_prime[9] = h_dot
#     state_prime[10] = tee_dot
#     state_prime[11] = 2*u'*u_prime
#
#     return state_prime
#
# end
#
# function KS_ODE_J2_test!(ẋ,x,p,t)
#     # scaling info
#     dscale = 1e6
#     tscale = 200
#
#     # unpack state
#     x = x[:]
#     p = x[1:4]
#     p_prime = x[5:8]
#     h = x[9]
#     #tee = state[10]
#     #r = state[11]
#     #k = state[12]
#
#     #J2 stuff
#     R_earth  = 6378137.0                   # m
#     mu = 3.986004418e14                    # m^3 s^{}-2}
#     X = x_from_u(p)*dscale                       # ECI, m
#     J2 =  1.081874e-3
#
#     # position vector in ECI from earth center to spacecraft
#     r_i = X[1]                             # ECI, m
#     r_j = X[2]                             # ECI, m
#     r_k = X[3]                             # ECI, m
#     R = norm(X)                          # ECI, m (dot(u,u)=norm(x))
#
#     # pertubation acceleration in ECI
#     # P_j2 = (-1.5*(mu*J2*R_earth^2)/(R^4)) * [(r_i/R)*(1-3*r_k^2/R^2);(r_j/R)*(1-3*r_k^2/R^2);(r_k/R)*(3-3*r_k^2/R^2)]
#     # P_j2 = P_j2*tscale^2/dscale
#     # # Levi-Civita transformation matrix
#     L = L_fx(p)
#     #
#     # # append 0 to end of P_j2 acceleration vector
#     # P_j2 = [P_j2;0]
#
#     # no J2
#     P_j2 = [0;0;0;0]
#
#     # ODE's
#     p_dp = -(h/2)*p + (R/(2*dscale))*L'*P_j2
#     #u_p = p_prime
#     tee_dot = dot(p,p)
#     h_dot = -2*dot(p_prime,L'*P_j2)
#
#     # index back into state_prime
#     ẋ[1:4] = p_prime
#     ẋ[5:8] = p_dp
#     ẋ[9] = h_dot
#     ẋ[10] = tee_dot
#     ẋ[11] = 2*p'*p_prime
#     ẋ[12] = 2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4])
#     return ẋ
# end
#
# function OE2ECI(a,e,i,OMEGA,omega,nu)
#     # inputs in meters and radians for all
#
#     mu_earth = 3.986004418e14#m
#
#     #semi latus rectum
#     p = a*(1-e^2)
#
#     #calculate radius
#     r = p/(1+e*cos(nu))
#
#     #perifocal coorinates, page 73 of the 1971 book
#     r_PQ = [r*cos(nu); r*sin(nu);0]
#
#     v_PQ = (sqrt(mu_earth/p))*[-sin(nu); (e+cos(nu));0]
#
#     #now lets transform using a rotaion matrix
#
#     rx_i = [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)]
#     rz_OMEGA = [cos(-OMEGA) sin(-OMEGA) 0 ; -sin(-OMEGA) cos(-OMEGA) 0; 0 0 1]
#     rz_omega = [cos(-omega) sin(-omega) 0 ; -sin(-omega) cos(-omega) 0; 0 0 1]
#
#     #now we transform
#     final_R = rz_OMEGA*rx_i*rz_omega
#
#     #now we get our final vectors
#     r_eci = final_R*r_PQ
#     v_eci = final_R*v_PQ
#
#     return r_eci, v_eci
#
# end
function get_time_transfer(X,dt,tscale)
    """get the time in days for a run. This allows us to remove t from state"""
    # X = states(solver)
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
function discretize_t(t_hist,n)
    """Take the time history vec, and discretize it even further"""
    # new time vec
    new_t = [0.0]
    for i = 1:(length(t_hist)-1)

        # time spacing betwen two knot points
        δt = t_hist[i+1] - t_hist[i]

        # add intermediate points between the knot points
        for j = 1:n
            push!(new_t,δt*j/n + t_hist[i])
        end
    end
    return new_t
end
function thrust_angle_from_rtn(u_rtn)
    u_rtn = normalize(u_rtn)
    α = atan(u_rtn[1],u_rtn[2])
    β = atan(u_rtn[3],norm(u_rtn[1:2]))
    return [α; β]
end
function get_rtn_in_out(U_real,r_eci_hist,v_eci_hist)

    return [thrust_angle_from_rtn(
                rECItoRTN([r_eci_hist[:,i];v_eci_hist[:,i]])*U_real[:,i] ) for i = 1:size(U_real,2)]
end
