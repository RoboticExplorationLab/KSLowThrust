using SatelliteDynamics, LinearAlgebra

function kron_delta(a::Real, b::Real)
    return a == b ? 1 : 0
end

"""
Internal helper function to aid in denormalization of gravity field coefficients.
Provides method co compute the factorial ratio: (n-m)!/(n+m)!
Arguments:
- `n::Real`: Gravity model degree, n.
- `m::Real`: Gravity model order, m.
Returns:
- `p::Real`: Factorial product
"""
function facprod(n::Int, m::Int)
    p = 1.0

    for i in (n-m+1):(n+m+1) # Is this correct? Should it be (n-m+1):(n+m)
        p = p/i
    end

    return p
end
"""
Compute the gravitational acceleration at a model given a spherical harmonic
gravity field model.
Arguments:
- `r::Array{<:Real, 1}`: Position of the point in the body (field) fixed frame. [m]
- `coef::Array{<:Real, 2}`: Gravity coefficients stored in dense matrix form. C_nm terns are stored along rows indexed from C_00 in coef[1, 1] to C_nm in coef[n+1, m+1]. S_nm terms are stored along matrix columns with S_n0 stored in coef[0, n+1]
- `n_max::Integer`: Maximum degree coefficient to use in expansion
- `m_max::Integer`: Maximum order coefficient to use in the expansion. Must be less than the degree.
- `r_ref::Real`: Reference distance of the gravity field.
- `GM::Real`: Gravitational constant of central body
- `normralized::Bool`: Whether the input gravity field coefficients are normalized coefficients (Default: true)
Returns:
- `a::Array{<:Real, 1}`: Acceleration in X, Y, and Z inertial directions [m/s^2]
"""
function spherical_harmonic_gravity2(r::Array{<:Real, 1}, coef::Array{<:Real, 2}, n_max::Integer, m_max::Integer, r_ref::Real, GM::Real; normalized::Bool=true)
    # Intermediate computations
    r_sqr = dot(r, r)
    rho   = r_ref^2/r_sqr
    x0    = r_ref * r[1] / r_sqr
    y0    = r_ref * r[2] / r_sqr
    z0    = r_ref * r[3] / r_sqr

    # Initialize Intermediary Matrices
    V = zeros(eltype(r), n_max+2, n_max+2)
    W = zeros(eltype(r), n_max+2, n_max+2)

    # Calculate zonal terms V(n, 0). Set W(n,0)=0.0
    V[0+1, 0+1] = r_ref /sqrt(r_sqr)
    W[0+1, 0+1] = 0.0

    V[1+1, 0+1] = z0 * V[0+1, 0+1]
    W[1+1, 0+1] = 0.0

    for n in 2:(n_max+1)
        V[n+1, 0+1] = ((2*n-1)*z0*V[n+1-1, 0+1] - (n-1)*rho*V[n+1-2,0+1])/n
        W[n+1, 0+1] = 0.0
    end

    # Calculate tesseral and sectoral terms
    for m in 1:m_max+1
        # Calculate V(m,m) to V(n_max+1,m)
        V[m+1, m+1] = (2*m-1)*(x0*V[m+1-1, m+1-1] - y0*W[m+1-1, m+1-1])
        W[m+1, m+1] = (2*m-1)*(x0*W[m+1-1, m+1-1] + y0*V[m+1-1, m+1-1])

        if m <= m_max
            V[m+1+1, m+1] = (2*m+1)*z0*V[m+1,m+1]
            W[m+1+1, m+1] = (2*m+1)*z0*W[m+1,m+1]
        end

        for n in (m+2):(n_max+1)
            V[n+1,m+1] = ((2*n-1)*z0*V[n+1-1,m+1]-(n+m-1)*rho*V[n+1-2,m+1])/(n-m)
            W[n+1,m+1] = ((2*n-1)*z0*W[n+1-1,m+1]-(n+m-1)*rho*W[n+1-2,m+1])/(n-m)
        end
    end

    # Calculate accelerations
    ax = 0.0
    ay = 0.0
    az = 0.0

    for m in 0:m_max
        for n in m:n_max
            C = 0.0
            S = 0.0
            if m == 0
                # Denormalize Coeeficients
                if normalized
                    N = sqrt(2*n+1)
                    C = N*coef[n+1, 0+1]
                else
                    C = coef[n+1, 0+1]
                end

                ax -= C*V[n+1+1, 1+1]
                ay -= C*W[n+1+1, 1+1]
                az -= (n+1)*C*V[n+1+1, 0+1]

            else
                if normalized
                    N = sqrt((2 - kron_delta(0,m))*(2*n+1)*facprod(n,m))
                    C = N*coef[n+1,   m+1]
                    S = N*coef[m+1-1, n+1]
                else
                    C = coef[n+1,   m+1]
                    S = coef[m+1-1, n+1]
                end

                fac =  0.5 * (n-m+1)*(n-m+2)
                ax += +0.5 * (-C * V[n+1+1, m+1+1] - S * W[n+1+1, m+1+1])
                      +fac * (+C * V[n+1+1, m+1-1] + S * W[n+1+1, m+1-1])
                ay += +0.5 * (-C * W[n+1+1, m+1+1] + S * V[n+1+1, m+1+1])
                      +fac * (-C * W[n+1+1, m+1-1] + S * V[n+1+1, m+1-1])
                az += (n-m+1)*(-C*V[n+1+1, m+1] - S*W[n+1+1, m+1])
            end
        end
    end

    a = (GM / (r_ref^2)) * [ax, ay, az]

    return a
end

# export accel_gravity
"""
Computes the accleration caused by Earth gravity as modeled by a spherical
harmonic gravity field.
Arguments:
- `r_sat::Array{<:Real, 1}`: Satellite position in the inertial frame [m]
- `R_eci_ecef::Array{<:Real, 2}`: Rotation matrix transforming a vector from the inertial to body-fixed reference frames.
- `n_max::Integer`: Maximum degree coefficient to use in expansion
- `m_max::Integer`: Maximum order coefficient to use in the expansion. Must be less than the degree.

Return:
- `a::Array{<:Real, 1}`: Gravitational acceleration in X, Y, and Z inertial directions [m/s^2]
References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.56-68.
"""
function accel_gravity2(x, R_eci_ecef, n_max::Int=5, m_max::Int=5)

    # Gravitational parameters of primary body
    # GM    = GRAVITY_MODEL.GM
    # R_ref = GRAVITY_MODEL.R

    # Body-fixed position
    r_bf = R_eci_ecef * x[1:3]

    # Compute spherical harmonic acceleration
    a_ecef = spherical_harmonic_gravity2(r_bf, GRAVITY_MODEL.data, n_max, m_max,
     GRAVITY_MODEL.R, GRAVITY_MODEL.GM, normalized=GRAVITY_MODEL.normalized)

    # Inertial acceleration
    a_eci = R_eci_ecef' * a_ecef

    # Finished
    return a_eci
end

function accel_point_mass2(x, gm_body::Real=GM_EARTH)
    # Restrict inputs to position only. Considered in body frame
    r  = x[1:3]

    # Acceleration
    a = -gm_body * r/norm(r)^3

    return a
end

function myRz(θ)
    return [cos(θ) sin(θ) 0;
           -sin(θ) cos(θ) 0;
            0      0      1]
end
function spherical_harmonic_gravity_wrapper(r_eci,R_eci_ecef)

    # R_eci_ecef = rECItoECEF(epc)

    # spherical harmonic gravity (ECI) - point mass gravity (ECI)
    return accel_gravity2(r_eci,R_eci_ecef) - accel_point_mass2(r_eci)
end
