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
    """Get rtn thrust angles from eci thrust vector"""
    return [thrust_angle_from_rtn(rECItoRTN(
    [r_eci_hist[:,i];v_eci_hist[:,i]])*U_real[:,i] ) for i = 1:size(U_real,2)]
end
