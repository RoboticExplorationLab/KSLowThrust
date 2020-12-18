using LinearAlgebra, StaticArrays, ModelingToolkit








function LQR_cost(Q,R,x,u,xf)
    return .5*(x-xf)'*Q*(x - xf) + .5*u'*R*u
end
function c_fx(x,u)
    # u_max = 1.0
    u_max = 2.0
    return [x[12:14];u_max]
end
function Π(x)
    v = SVector(x[1],x[2],x[3])
    s = x[4]
    nv = norm(v)

    # if     nv <= -s
        # return zeros(4)
    # elseif nv <= s
        return x
    # elseif nv > abs(s)
        # return SVector{4}(.5*(1 + s/nv)*[v;nv])
    # else
    #     @infiltrate
    # end
end

function Lag(Q,R,x,u,xf,λ,μ)

    proj_gap = Π(λ - μ*c_fx(x,u))

    return (LQR_cost(Q,R,x,u,xf) +
            (1/(2*μ))*(  sum(proj_gap.^2) - sum(λ.^2))
    # return LQR_cost(Q,R,x,u,xf)
end
