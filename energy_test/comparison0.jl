using LinearAlgebra, StaticArrays, Attitude, SatelliteDynamics
using Infiltrator, ProgressMeter
using MATLAB

function initial_conditions()

    # epc0 = Epoch(2019, 1, 1, 12, 0, 0, 0.0)

    # Declare initial state in terms of osculating orbital elements
    rp = R_EARTH + 200e3;ra = R_EARTH + 35000e3
    a = (ra+rp)/2
    e = (ra-rp)/(ra+rp)
    i = deg2rad(27)
    # i = 0
    ω = deg2rad(34)
    Ω = deg2rad(133)
    M = deg2rad(180)
    oe0  = [a,e,i,Ω,ω,M]

    # Convert osculating elements to Cartesean state
    eci0 = sOSCtoCART(oe0, use_degrees=false)

    # modified equinoctial elements
    p = a*(1-e^2)
    f = e*cos(ω + Ω)
    g = e*sin(ω + Ω)
    h = tan(i/2)*cos(Ω)
    k = tan(i/2)*sin(Ω)
    L = wrap_to_2pi(ω + Ω + M)

    eoe0 = [p;f;g;h;k;L]

    r_eci = eci0[1:3]
    v_eci = eci0[4:6]

    p_0 = u_from_x(r_eci)                          # u
    p_prime_0 = uprime_from_xdot(v_eci,p_0)     # du/ds
    h_0 = GM_EARTH/norm(r_eci) - .5*norm(v_eci)^2         # m^2/s^2
    t_0 = 0

    ks0 = [p_0;p_prime_0;h_0;t_0]

    return eci0, eoe0, ks0
end
function rk4(f::Function, x_n, dt)
    k1 = dt*f(x_n)
    k2 = dt*f(x_n+k1/2)
    k3 = dt*f(x_n+k2/2)
    k4 = dt*f(x_n+k3)
    return (x_n + (1/6)*(k1+ 2*k2 + 2*k3 + k4))
end
function fode_j2(r_eci)
    r = norm(r_eci)
    x,y,z = r_eci
    J2 = J2_EARTH
    μ = GM_EARTH

    # precompute repeated stuff
    Re_r_sqr = 1.5*J2*(R_EARTH/r)^2
    five_z_sqr = 5*z^2/r^2

    return  SA[(μ/r^3)*x*(Re_r_sqr*(five_z_sqr - 1)),
               (μ/r^3)*y*(Re_r_sqr*(five_z_sqr - 1)),
               (μ/r^3)*z*(Re_r_sqr*(five_z_sqr - 3))]
end
function cartesian_dynamics(x)
    r_eci = SA[x[1],x[2],x[3]]
    v_eci = SA[x[4],x[5],x[6]]

    # relavent parameters
    J2 = J2_EARTH
    μ = GM_EARTH

    # eci position stuff
    r = norm(r_eci)
    x,y,z = r_eci

    # precompute repeated stuff
    Re_r_sqr = 1.5*J2*(R_EARTH/r)^2
    five_z_sqr = 5*z^2/r^2

    return  SA[v_eci[1],v_eci[2],v_eci[3],
               (-μ/r^3)*x*(1 - Re_r_sqr*(five_z_sqr - 1)),
               (-μ/r^3)*y*(1 - Re_r_sqr*(five_z_sqr - 1)),
               (-μ/r^3)*z*(1 - Re_r_sqr*(five_z_sqr - 3))]
end
function J2_from_EOE(x)
    """J2 acceleration in RTN from equinoctial orbital elements.
    Args:
        p: a*(1-e^2) semi-latus rectum       (m)
        f: e*cos(ω + Ω)                      (non-dimensional)
        g: e*sin(ω + Ω)                      (non-dimensional)
        h: tan(i/2)*cos(Ω)                   (non-dimensional)
        k: tan(i/2)*sin(Ω)                   (non-dimensional)
        L: ω + Ω + θ argument of latitude    (rad)
    Returns:
        a_rtn: acceleration in the RTN frame  (m/s)
    """

    # unpack state
    p,f,g,h,k,L = x
    J2 =  J2_EARTH

    # precompute trig functions
    sinL = sin(L)
    cosL = cos(L)

    # constants for EOE
    mu = GM_EARTH
    w = 1+f*cosL + g*sinL
    r = p/w
    Re = R_EARTH

    # accelerations
    a_r = -((3*mu*J2*Re^2)/(2*r^4))*(1 -   (  (12*(h*sinL-k*cosL)^2 )
                                                    /   ((1+h^2 + k^2)^2))  )
    a_t = -((12*mu*J2*Re^2)/(r^4))*((  ((h*sinL - k*cosL)*(h*cosL+k*sinL) )
                                                    /  ((1+h^2 + k^2)^2))  )
    a_n = -((6*mu*J2*Re^2)/(r^4))*((  ((h*sinL - k*cosL)*(1-h^2-k^2) )
                                                    /  ((1+h^2 + k^2)^2))  )
    return a_r, a_t, a_n
end
function eoe_dynamics(x)
    """Equations of motion for equinoctial elements in 2-body problem.
    Args:
        epc: time(::epoch)
        x: state
            - p: a*(1-e^2) semi-latus rectum       (m)
            - f: e*cos(ω + Ω)                      (non-dimensional)
            - g: e*sin(ω + Ω)                      (non-dimensional)
            - h: tan(i/2)*cos(Ω)                   (non-dimensional)
            - k: tan(i/2)*sin(Ω)                   (non-dimensional)
            - L: ω + Ω + θ argument of latitude    (rad)
            - mass: spacecraft mass                (kg)
        u: thruster acceleration                   (m/s^2)
        thruster_on: bool for thruster             ()
    Returns:
        ̇x: state derivative
    """

    # unpack state
    p,f,g,h,k,L = x

    # equinoctial parameters for GVE's
    sinL = sin(L)
    cosL = cos(L)
    alpha_2 = h^2 - k^2
    s2 = 1+h^2 + k^2
    w = 1+f*cosL + g*sinL
    r = p/w

    # precompute to save calculations
    sqrtpu = sqrt(p/GM_EARTH)

    # J2
    u_r, u_t, u_n = J2_from_EOE(x)

    # GVE's for equinoctial orbital elements
    p_dot = (2*p/w)*sqrtpu*u_t
    f_dot = sqrtpu*( u_r*sinL    +    ((w+1)*cosL + f)*u_t/w    -
                                            (h*sinL - k*cosL)*g*u_n/w  )
    g_dot = sqrtpu*(-u_r*cosL   +    ((w+1)*sinL + g)*u_t/w    +
                                            (h*sinL - k*cosL)*f*u_n/w  )
    h_dot = sqrtpu*(s2*u_n/(2*w))*cosL
    k_dot = sqrtpu*(s2*u_n/(2*w))*sinL
    L_dot = sqrt(GM_EARTH*p)*(w/p)^2 + (1/w)*sqrtpu*(h*sinL - k*cosL)*u_n

    # return state derivative
    return SA[p_dot,f_dot,g_dot,h_dot,k_dot,L_dot]
end
function rv_from_eoe(eoe_vec)
    """Position and Velocity resolved in ECI
    Args:
        p: a*(1-e^2) semi-latus rectum     (m)
        f: e*cos(ω + Ω)                    (non-dimensional)
        g: e*sin(ω + Ω)                    (non-dimensional)
        h: tan(i/2)*cos(Ω)                 (non-dimensional)
        k: tan(i/2)*sin(Ω)                 (non-dimensional)
        L: ω + Ω + θ argument of latitude  (rad)
    Returns:
        r: position vector from eci origin to spacecraft, resolved in eci. (m)
        v: velocity vector of spacecraft in eci, resolved in eci.          (m/s)
    """

    # unpack eoe
    p,f,g,h,k,L = eoe_vec

    # precompute sin and cosine
    sinL = sin(L)
    cosL = cos(L)

    # various eoe specific parameters
    alf2 = h^2 - k^2
    s2 = 1 + h^2 + k^2
    w = 1 + f*cosL + g*sinL
    r = p/w

    # eci position vector
    r = SA[(r/s2) * (  cosL + alf2*cosL + 2*h*k*sinL),
         (r/s2) * (  sinL - alf2*sinL + 2*h*k*cosL),
         (2*r/s2)*(h*sinL -   k *cosL)]

    # ECI velocity vector
    v = sqrt(GM_EARTH/p)*SA[-(1/(s2)) * (sinL + alf2*sinL - 2*h*k*cosL + g - 2*f*h*k + alf2*g),
         -(1/(s2))* (-cosL + alf2*cosL + 2*h*k*sinL - f + 2*g*h*k + alf2*f),
          (2/s2) * (h*cosL +    k*sinL + f*h + g*k)]

    # @infiltrate
    return SA[r[1],r[2],r[3],v[1],v[2],v[3]]
end

function ks_dynamics(x)

    # unpack state
    p = SVector(x[1],x[2],x[3],x[4])
    p_prime = SVector(x[5],x[6],x[7],x[8])
    h = x[9]
    t = x[10]

    # position in eci
    r_eci = x_from_u(p)

    # j2 accel
    a_j2 = fode_j2(r_eci)

    L = L_fx(p)

    # append 0 to end of P_j2 acceleration vector
    P_j2 = SA[a_j2[1],a_j2[2],a_j2[3],0]

    p_dp = -(h/2)*p + ((dot(p,p))/(2))*L'*P_j2

    return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
                   p_dp[1],p_dp[2],p_dp[3],p_dp[4],
                   -2*dot(p_prime,L'*P_j2),
                   dot(p,p))
end
function rv_from_ks(p,p_prime)
    x = x_from_u(p)
    xdot = xdot_from_u(p,p_prime)
    return SA[x[1],x[2],x[3],xdot[1],xdot[2],xdot[3]]
end
function ϵ(rv)
    r = SA[rv[1],rv[2],rv[3]]
    v = SA[rv[4],rv[5],rv[6]]
    return dot(v,v)/2 - GM_EARTH/norm(r)
end
function P2(x)
    return .5*(3*x^2 -1)
end
function true_anom_from_rv(r,v)
    h = cross(r,v)
    e = (v × h)/GM_EARTH - r/norm(r)

    ν = acos(dot(e,r)/(norm(e)*norm(r)))
    if dot(r,v)<0
        ν = 2*pi - ν
    end
    return ν
end
function sma_from_eci(rv)
    rv = convert(Array{Float64},rv)

    try
        r = rv[1:3]
        v = rv[4:6]
        oe = sCARTtoOSC(rv,use_degrees = false)
        a,e,i,Ω,ω,M = oe
        θ = true_anom_from_rv(r,v)
        cosϕ = sin(i)*sin(θ + ω)
        U_j2 = - (GM_EARTH/norm(r))*J2_EARTH*(R_EARTH/norm(r))^2*P2(cosϕ)
        ϵ = dot(v,v)/2 - GM_EARTH/norm(r) - U_j2
        # return oe[1]
        return ϵ
    catch err
        return NaN
    end

    return oe[1]
end
function driver()
    eci0, eoe0, ks0 = initial_conditions()

    num_orbits = 6
    tf = num_orbits*36951.67397573928
    # N = 1000
    N = 600
    # N = Int(round((tf/(24*3600))*80))
    @show N
    dt = tf/N

    t_vec = [(i-1)*dt for i = 1:N]

    X_eci = [@SVector zeros(6) for i = 1:N]
    X_eci[1] = SVector{6}(eci0)
    X_eoe = [@SVector zeros(6) for i = 1:N]
    X_eoe[1] = eoe0
    X_ks = [@SVector zeros(10) for i = 1:N]
    X_ks[1] = ks0

    # one day in fictitious time (this was guess and check)
    sf = num_orbits*0.001563905
    ds = sf/N

    @inbounds for i = 1:(N-1)
        X_eci[i+1] = rk4(cartesian_dynamics,X_eci[i],dt)
        X_eoe[i+1] = rk4(eoe_dynamics,X_eoe[i],dt)
        X_ks[i+1] = rk4(ks_dynamics,X_ks[i],ds)
    end

    ksm = mat_from_vec(X_ks)
    goal_tf = 36951.67397573928

    @show goal_tf - X_ks[end][10]

    ks_t_vec = vec(ksm[10,:])

    X_eoe_eci = [rv_from_eoe(X_eoe[i]) for i = 1:length(X_eoe)]
    X_ks_eci = [rv_from_ks(X_ks[i][1:4],X_ks[i][5:8]) for i = 1:length(X_ks)]

    eoe_m = mat_from_vec(X_eoe_eci)
    eci_m = mat_from_vec(X_eci)
    ks_m = mat_from_vec(X_ks_eci)

    # SMA from all
    eoe_a= [sma_from_eci(X_eoe_eci[i]) for i = 1:length(X_eoe_eci)]
    eci_a = [sma_from_eci(X_eci[i]) for i = 1:length(X_eci)]
    ks_a = [sma_from_eci(X_ks_eci[i]) for i = 1:length(X_ks_eci)]

    T_period = 36951.67397573928
    t_vec /= T_period
    ks_t_vec /= T_period
    eps_0 = eci_a[1]
    pl1 = 1000*100*(1 .- eoe_a /eps_0)
    pl2 = 1000*100*(1 .- eci_a /eps_0)
    pl3 = 1000*100*(1 .- ks_a /eps_0)
    mat"
    figure
    hold on
    plot($t_vec , $pl1,'linewidth',3)
    plot($t_vec , $pl2)
    plot($ks_t_vec, $pl3)
    legend('EOE','Cartesian','KS')
    %xlim([0 5])
    %ylim([-1.5 2.8])
    %xticks([0 1 2 3 4 5])
    %xlabel('Number of Orbits')
    %ylabel('Specific Orbital Energy (kJ/kg)')
    %set(gca, 'XScale', 'log')
    addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
    matlab2tikz('energy_test3.tikz')
    hold off
    "
    # using JLD2
    # output_tup = (t_vec = t_vec, ks_t_vec = ks_t_vec, eoe_a, eci_a, ks_a)
    # @save "energy_data.jld2" output_tup
    # X_eoe_eci = [rv_from_eoe(X_eoe[i]) for i = 1:length(X_eoe)]
    # X_ks_eci = [rv_from_ks(X_ks[i][1:4],X_ks[i][5:8]) for i = 1:length(X_ks)]
    #
    # eoe_m = mat_from_vec(X_eoe_eci)
    # eci_m = mat_from_vec(X_eci)
    # ks_m = mat_from_vec(X_ks_eci)
    #
    # # SMA from all
    # eoe_a= [sma_from_eci(X_eoe_eci[i]) for i = 1:length(X_eoe_eci)]
    # eci_a = [sma_from_eci(X_eci[i]) for i = 1:length(X_eci)]
    # ks_a = [sma_from_eci(X_ks_eci[i]) for i = 1:length(X_ks_eci)]
    #
    # # energies
    # ks_ϵ = ϵ(X_ks_eci[end])
    # eci_ϵ = ϵ(X_eci[end])
    # eoe_ϵ = ϵ(X_eoe_eci[end])
    #
    # @show ks_ϵ
    # @show eci_ϵ
    # @show eoe_ϵ
    #
    #
    # mat"
    # figure
    # hold on
    # plot3($eci_m(1,:),$eci_m(2,:),$eci_m(3,:))
    # plot3($ks_m(1,:),$ks_m(2,:),$ks_m(3,:))
    # hold off
    # "
    #
    #
    # # # @show "yes"
    # # eci_m = mat_from_vec(X_eci)
    # # eoe_m = mat_from_vec(X_eoe_eci)
    # # mat"
    # # figure
    # # hold on
    # # plot3($eci_m(1,:),$eci_m(2,:),$eci_m(3,:))
    # # plot3($eoe_m(1,:),$eoe_m(2,:),$eoe_m(3,:))
    # # hold off
    # # "
    # #
    # #
    # # @show norm(X_eci[end] - X_eoe_eci[end])
    # # @show norm(X_eci[end] - X_ks_eci[end])
end

driver()
