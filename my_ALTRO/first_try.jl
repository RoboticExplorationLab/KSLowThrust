using ForwardDiff, LinearAlgebra, MATLAB, Infiltrator
using Attitude, StaticArrays



function dynamics(x,u,t)

    p = SVector(x[1],x[2],x[3])
    ω = SVector(x[4],x[5],x[6])

    J = @views params.J
    invJ = @views params.invJ

    return SVector{6}([pdot_from_w(p,ω);
           invJ*(u - ω × (J*ω)]
end
function rk4(f, x_n, u, t,dt)

    x_n = SVector{nx}(x_n)

    k1 = h*f(x_n,u,t)
    k2 = h*f(x_n+k1/2, u,t+dt/2)
    k3 = h*f(x_n+k2/2, u,t+dt/2)
    k4 = h*f(x_n+k3, u+dt)


    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

end
function discrete_dynamics(x,u,t,dt)
    return rk4(dynamics,x,u,t,dt)
end

struct AL_struct
    λ::Array{Array{Float64,1},1}
    μ::Float64
end



function iLQR(x0,xg,Q,R,Qf,N,dt,Ad,Bd,utraj,λ_vec,μ)

    # state and control dimensions
    # Nx, Nu = size(Bd)
    Nx = length(x0)
    Nu = size(R,1)

    # initial trajectory is initial conditions the whole time
    xtraj = fill(x0,N)
    # utraj = fill(zeros(Nu),N-1)
    J = 0.0
    for k = 1:N-1
        xtraj[k+1] = discrete_dynamics(xtraj[k],utraj[k],(i-1)*dt,dt)
        J += L(Q,R,xtraj[k],utraj[k],λ_vec[k],μ)
    end
    J += quad(Qf,xtraj[N])


    @show J
    # allocate K and l
    K = fill(zeros(Nu,Nx),N-1)
    l = fill(zeros(Nu),N-1)
    λ_vec = fill(zeros(Nu +1),N-1)
    μ = 1.0

    # allocate the new states and controls
    xnew = fill(zeros(Nx),N)
    unew = fill(zeros(Nu),N-1)

    # main loop
    for iter = 1:50

        # cost to go matrices at the end of the trajectory
        S = Qf
        s = Qf*(xtraj[N]-xg)

        # backwards pass
        for k = (N-1):-1:1

            # calculate cost gradients for this time step
            q = Q*(xtraj[k]-xg)
            # r = R*utraj[k]
            x = xtraj[k]
            u = utraj[k]
            λ = λ_vec[k]
            _Lx(x) = L(Q,R,x,u,λ,μ)
            _Lu(u) = L(Q,R,x,u,λ,μ)

            # Q = ForwardDiff.hessian(_Lx, x)
            # q = ForwardDiff.gradient(_Lx, x)
            R = ForwardDiff.hessian(_Lu, u)
            r = ForwardDiff.gradient(_Lu, u)

            # jacobians
            Ak, Bk = Ad, Bd

            # linear solve
            l[k] = (R + Bk'*S*Bk)\(r + Bk'*s)
            K[k] = (R + Bk'*S*Bk)\(Bk'*S*Ak)
            # update
            Snew = Q + K[k]'*R*K[k] + quad(S,(Ak-Bk*K[k]))
            snew = q - K[k]'*r + K[k]'*R*l[k] + (Ak-Bk*K[k])'*(s - S*Bk*l[k])

            # update S's
            S = copy(Snew)
            s = copy(snew)
            @show k

        end
        # @infiltrate
        # error()
        # initial conditions
        xnew[1] = copy(x0)

        # learning rate
        alpha = 1.0

        # forward pass
        Jnew = 0.0
        for inner_i = 1:6
            @show inner_i
            Jnew = 0.0
            # rollout the dynamics
            for k = 1:N-1
                unew[k] = utraj[k] - alpha*l[k] - K[k]*(xnew[k]-xtraj[k])
                xnew[k+1] = Ad*xnew[k] + Bd*unew[k]
                Jnew += L(Q,R,xnew[k],unew[k],λ_vec[k],μ)
            end
            Jnew += quad(Qf,xnew[N])
            @infiltrate
            error()
            # if the new cost is lower, we keep it
            if Jnew<J
                break
            else# this also pulls the linesearch back if hasnan(xnew)
                alpha = (1/2)*alpha
            end
        end

        # update trajectory and control history
        xtraj = copy(xnew)
        utraj = copy(unew)


    end

        return xtraj, utraj, K
end


x0 = [10;15;20.0;0;0;0]
xg = zeros(6)
Q = 1*I(6)
Qf = copy(Q)
N = 50
R = 50*I(3)
utraj = fill(zeros(3),N-1)
λ_vec = fill(zeros(4),N-1)
μ = 1.0
x, u, K = iLQR(x0,xg,Q,R,Qf,N,dt,Ad,Bd,utraj,λ_vec,μ)




# x = mat_from_vec(x)
# u = mat_from_vec(u)
#
# mat"
# figure
# hold on
# plot($x')
# "
#
# mat"
# figure
# hold on
# plot($u')
# "
