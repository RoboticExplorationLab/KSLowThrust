using ForwardDiff, LinearAlgebra, MATLAB, Infiltrator
using Attitude, StaticArrays
const FD = ForwardDiff
function cfill(nx,N)
    return [zeros(nx) for i = 1:N]
end
function cfill(nu,nx,N)
    return [zeros(nu,nx) for i = 1:N]
end

function dynamics(x, u_dot,t)
    """ Here is the state
            p  ∈ R⁴ [1:4]
            p′ ∈ R⁴ [5:8]
            h  ∈ R  [9]
           |r| ∈ R  [10]
            r₃ ∈ R  [11]
            u  ∈ R³ [12:14]
    """

    # unpack state
    p = SVector(x[1],x[2],x[3],x[4])
    p_prime = SVector(x[5],x[6],x[7],x[8])
    h = x[9]



    P_j2 = SVector(x[12],x[13],x[14],0)*(params.tscale^2/params.dscale)/params.uscale
    # Levi-Civita transformation matrix
    L = L_fx(p)

    # append 0 to end of P_j2 acceleration vector
    # P_j2 = SVector{4}([P_j2;0]*params.tscale^2/params.dscale)

    p_dp = -(h/2)*p + ((dot(p,p))/(2))*L'*P_j2

    return SVector(p_prime[1],p_prime[2],p_prime[3],p_prime[4],
                   p_dp[1],p_dp[2],p_dp[3],p_dp[4],
                   -2*dot(p_prime,L'*P_j2),
                   2*p'*p_prime,
                   2*(p_prime[1]*p[3] + p[1]*p_prime[3] + p_prime[2]*p[4] + p[2]*p_prime[4]),
                   u_dot[1],
                   u_dot[2],
                   u_dot[3])

end
function rk4(f, x_n, u, t,dt)

    # x_n = SVector{6}(x_n)

    k1 = dt*f(x_n,u,t)
    k2 = dt*f(x_n+k1/2, u,t+dt/2)
    k3 = dt*f(x_n+k2/2, u,t+dt/2)
    k4 = dt*f(x_n+k3, u,t+dt)


    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

end
function discrete_dynamics(x,u,t,dt)
    return rk4(dynamics,x,u,t,dt)
end
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

    if     nv <= -s
        return zeros(4)
    elseif nv <= s
        return x
    elseif nv > abs(s)
        return SVector{4}(.5*(1 + s/nv)*[v;nv])
    else
        @infiltrate
    end
end

function Lag(Q,R,x,u,xf,λ,μ)

    return (LQR_cost(Q,R,x,u,xf) +
            (1/(2*μ))*(  norm(Π(λ - μ*c_fx(x,u)))^2 - dot(λ,λ)))
    # return LQR_cost(Q,R,x,u,xf)
end

function ilqr(x0,utraj,xf,Q_lqr,Qf_lqr,R_lqr,N,dt,μ,λ)

    # state and control dimensions
    Nx = length(x0)
    Nu = size(R_lqr,1)

    # initial trajectory is initial conditions the whole time
    xtraj = cfill(Nx,N)
    xtraj[1] = copy(x0)

    J = 0.0
    for k = 1:N-1
        xtraj[k+1] = discrete_dynamics(xtraj[k],utraj[k],(k-1)*dt,dt)
        J += Lag(Q_lqr,R_lqr,xtraj[k],utraj[k],xf,λ[k],μ)
    end
    J += .5*(xtraj[N]-xf)'*Qf_lqr*(xtraj[N] - xf)

    # allocate K and l
    K = cfill(Nu,Nx,N-1)
    l = cfill(Nu,N-1)

    # initial stuff
    xtraj_initial = copy(xtraj)
    utraj_initial = copy(utraj)
    solver_failed = false

    # allocate the new states and controls
    xnew = cfill(Nx,N)
    unew = cfill(Nu,N-1)

    # allocate gradients and hessians
    Q = zeros(Nx,Nx)
    q = zeros(Nx)
    R = zeros(Nu,Nu)
    r = zeros(Nu)
    S = zeros(Nx,Nx)
    Snew = zeros(Nx,Nx)
    s = zeros(Nx)
    snew = zeros(Nx)


    # allocate jacobians
    Ak = zeros(Nx,Nx)
    Bk = zeros(Nx,Nu)

    # main loop
    for iter = 1:3000

        # cost to go matrices at the end of the trajectory
        S .= Qf_lqr
        s .= Qf_lqr*(xtraj[N]-xf)

        # backwards pass
        for k = (N-1):-1:1

            # calculate cost gradients for this time step
            _L_fx(xl) = Lag(Q_lqr,R_lqr,xl,utraj[k],xf,λ[k],μ)
            FD.hessian!(Q,_L_fx,xtraj[k])
            FD.gradient!(q,_L_fx,xtraj[k])

            _L_fx2(ul) = Lag(Q_lqr,R_lqr,xtraj[k],ul,xf,λ[k],μ)
            FD.hessian!(R,_L_fx2,utraj[k])
            FD.gradient!(r,_L_fx2,utraj[k])

            # discrete dynamics jacobians
            _A_fx(x) = discrete_dynamics(x,utraj[k],(k-1)*dt,dt)
            FD.jacobian!(Ak,_A_fx,xtraj[k])
            _B_fx(u) =  discrete_dynamics(xtraj[k],u,(k-1)*dt,dt)
            FD.jacobian!(Bk,_B_fx,utraj[k])

            # linear solve
            F = factorize(R + Bk'*S*Bk)
            l[k] .= F\(r + Bk'*s)
            K[k] .= F\(Bk'*S*Ak)

            # update
            Snew .= Q + K[k]'*R*K[k] + (Ak-Bk*K[k])'*S*(Ak-Bk*K[k])
            snew .= q - K[k]'*r + K[k]'*R*l[k] + (Ak-Bk*K[k])'*(s - S*Bk*l[k])

            # update S's
            S .= Snew
            s .= snew
        end

        # initial conditions
        xnew[1] .= copy(x0)

        # learning rate
        alpha = 1.0

        # forward pass
        Jnew = 0.0
        for inner_i = 1:10

            Jnew = 0.0
            # rollout the dynamics
            for k = 1:N-1
                unew[k] = utraj[k] - alpha*l[k] - K[k]*(xnew[k]-xtraj[k])
                xnew[k+1] = discrete_dynamics(xnew[k],unew[k],(k-1)*dt,dt)
                if hasnan(xnew[k+1])
                    Jnew += Inf
                else
                    Jnew += Lag(Q_lqr,R_lqr,xnew[k],unew[k],xf,λ[k],μ)
                end
            end
            Jnew += .5*(xnew[N]-xf)'*Q_lqr*(xnew[N] - xf)

            # if the new cost is lower, we keep it
            if Jnew<J
                break
            else# this also pulls the linesearch back if hasnan(xnew)
                alpha = (1/2)*alpha
            end
            if inner_i == 10
                # error("Line Search Failed")
                print("line search failed")

                solver_failed = true
                # we do this to get out of it
                Jnew = copy(J)
                break
            end

        end

        # update trajectory and control history
        xtraj .= xnew
        utraj .= unew

        # termination criteria
        if abs(J - Jnew)<params.dJ_tol
            break
        end
        # else
        #     # update trajectory and control history
        #     xtraj .= xnew
        #     utraj .= unew
        # end


        # ----------------------------output stuff-----------------------------
        if rem((iter-1),4)==0
            println("iter       alpha      maxL    Cost")
        end
        maxL = round(maximum(vec(mat_from_vec(l))),digits = 3)
        J = Jnew
        J_display = round(J,digits = 3)
        alpha_display = round(alpha,digits = 3)
        println("$iter          $alpha_display      $maxL    $J_display")


    end

    if solver_failed
        return xtraj_initial, utraj_initial, K
    else
        return xtraj, utraj, K
    end
end


function runit()
x0 = vec([0.1, 0.0, 0.0, 0.8048687470637682, 0.0, 4.647638402647582,
-6.869088938897239, 0.0, 33.24699010602865, 0.6578137000000001, 0.0,0,0,0])
xf = SVector{14}([zeros(9);4.2;zeros(4)])
scaling_num = 10
R_scale = .4
udot_scale = 1e-4
Q = Diagonal(SVector{14}([1e-8*ones(9);scaling_num;10;R_scale*ones(3)]))
R = Diagonal(SVector(udot_scale,udot_scale,udot_scale))
Qf = 10*Q
dscale = 1e7 # m
tscale = 20000 # s
uscale = 1000.0
global params = (tscale = tscale, uscale = uscale, dscale = dscale,
dJ_tol = 1e-1)
dt = 5.0e-2
N = 500


μ = 100
ϕ = 3
λ = cfill(4,N-1)
utraj = cfill(3,N-1)
xtraj = cfill(12,N-1)
# constraint_violation = zeros(6,N-1)
for i = 1:10
    xtraj, utraj, K = ilqr(x0,utraj,xf,Q,Qf,R,N,dt,μ,λ)
    for k = 1:length(λ)
        λ[k] = Π( λ[k] - μ*c_fx(xtraj[k],utraj[k]) )
    end
    @info "here is new out loop iteration"
    @show i
    μ*=ϕ

    # xm = mat_from_vec(xtraj)
    u_norm = [norm(xtraj[i][12:14]) for i = 1:length(xtraj)]
    con_violation = max(0,maximum(u_norm) - 2)
    @info "max constraint violation is $con_violation"
    #
    # if con_violation < 1e-4
    #     break
    # end
end

r_eci_hist = mat_from_vec([x_from_u(xtraj[i]) for i = 1:length(xtraj)])


mat"
figure
hold on
plot3($r_eci_hist(1,:),$r_eci_hist(2,:),$r_eci_hist(3,:))
%axis equal
hold off
"
xm = mat_from_vec(xtraj)
um = mat_from_vec(utraj)
# mat"
# figure
# hold on
# plot($xm')
# hold off
# "
mat"
figure
title('U dot')
hold on
plot($um')
hold off
"
mat"
figure
title('U')
hold on
plot($xm(12:14,:)')
hold off
"
u_norm = [norm(xtraj[i][12:14]) for i = 1:length(xtraj)]

mat"
figure
hold on
plot($u_norm)
hold off
"



    return xm, um
end

xm,um = runit()

# using JLD2
# @save "sixty_5_days_8.jld2" xm um

function get_time_transfer(xm,dt,tscale)
    """get the time in days for a run. This allows us to remove t from state"""
    X = vec_from_mat(xm)
    t = 0.0
    for i = 1:(length(X)-1)

        p0 = X[i][1:4]
        R0 = dot(p0,p0)

        p1 = X[i+1][1:4]
        R1 = dot(p1,p1)

        t += (R0 + R1)*dt/2
    end
    t_days = t*tscale/(24*3600)
    return t_days
end

t_days = get_time_transfer(xm,5e-2,params.tscale)
