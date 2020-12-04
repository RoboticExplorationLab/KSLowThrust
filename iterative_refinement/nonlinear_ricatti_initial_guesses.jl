# using Pkg; cd(joinpath(dirname(@__DIR__),"my_ALTRO"))
# Pkg.activate(".")

using LinearAlgebra, SparseArrays, SuiteSparse, Attitude
using ForwardDiff
using StaticArrays
const FD = ForwardDiff

# random utilities
function dynamics(x,u,t)

    p = SVector(x[1],x[2],x[3])
    ω = SVector(x[4],x[5],x[6])

    return SVector{6}([pdot_from_w(p,ω);
           invJ*(u - ω × (J*ω))])
end
function rk4(f, x_n, u, t,dt)
    x_n = SVector{6}(x_n)

    k1 = dt*f(x_n,u,t)
    k2 = dt*f(x_n+k1/2, u,t+dt/2)
    k3 = dt*f(x_n+k2/2, u,t+dt/2)
    k4 = dt*f(x_n+k3, u,t+dt)
    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

end
function discrete_dynamics(x,u,t,dt)
    return rk4(dynamics,x,u,t,dt)
end

# forward rollout
x0 = [p_from_phi(deg2rad(10)*normalize(randn(3)));
      deg2rad(10)*normalize(randn(3))]
nx = 6
nu = 3
dt = 0.1
N = 100
J = Diagonal(SVector(1,2,3))
invJ = inv(J)

function rollout(x0,utraj,dt,N)
    u = zeros(3)
    X = cfill(nx,N)
    X[1] = x0
    for i = 1:N-1
        X[i+1] = discrete_dynamics(X[i],utraj[i],0,dt)
    end
    return X
end

# initial controls and states
utraj = [randn(nu) for i = 1:N-1]
xtraj = rollout(x0,utraj,dt,N)

# cost function
Q = sparse(Array(float(I(nx))))
Qf = 10*Q
R = 2*sparse(Array(float(I(nu))))

function gen_AB(xtraj,utraj,dt,N)
    """Create the A and B matrices for the trajectory"""
    A = cfill(nx,nx,N)
    B = cfill(nx,nu,N)
    for i = 1:N-1
        # jacobians
        _A_fx(x) = discrete_dynamics(x,utraj[i],(i-1)*dt,dt)
        A[i] = FD.jacobian(_A_fx,xtraj[i])
        _B_fx(u) =  discrete_dynamics(xtraj[i],u,(i-1)*dt,dt)
        B[i] = FD.jacobian(_B_fx,utraj[i])
    end
    return A,B
end

A,B = gen_AB(xtraj,utraj,dt,N)


function kkt_stuff(A,B,nx,nu,x0,N,Q,Qf,R)
    """Solve for the updated dU and dX using direct collocation"""

    # number of constraints
    nc = (N-1)*nx

    # length of the [x0 u0 x1 u1 ... xN] vector
    nz = (N)*nx + (N-1)*nu
    idx_x = [(t-1)*(nx+nu) .+ (1:nx) for t = 1:N]
    idx_u = [(t-1)*(nx+nu) .+ ((nx+1):(nx+nu)) for t = 1:N-1]

    # indexing for the constraint vector
    idx_cons = [(t-1)*(nx) .+ (1:nx) for t = 1:(N-1)]

    # constraint jacobian
    C = spzeros(nc,nz)
    for i = 1:N-1
        # Aᵢδxᵢ + Bᵢδuᵢ - δx\_{i+1} = 0
        C[idx_cons[i],idx_x[i]] = A[i]
        C[idx_cons[i],idx_u[i]] = B[i]
        C[idx_cons[i],idx_x[i+1]] = -I(nx)
    end
    C_lower = spzeros(nx,nz)
    C_lower[1:nx,1:nx] = sparse(I(nx))
    C = [C;C_lower]
    d = zeros(nc + nx)

    # initial condition
    d[(end-nx+1):end] = zeros(nx)

    # cost hessian
    P = sparse(blockdiag(kron(sparse(I(N-1)),blockdiag(Q,R)),Qf))

    # cost gradient
    q = zeros(nz)
    for i = 1:N-1
        q[idx_x[i]] = Q*xtraj[i]
        q[idx_u[i]] = R*utraj[i]
    end
    q[idx_x[N]] = Qf*xtraj[N]

    # KKT system
    kkt_A = [P C';
             C zeros(size(C,1),size(C,1))]

    kkt_b = [-q;d]

    # solve the KKT system
    z_sol = kkt_A\kkt_b

    # unpack the dX and dU from the solution
    dX = [z_sol[idx_x[i]] for i = 1:length(idx_x)]
    dU = [z_sol[idx_u[i]] for i = 1:length(idx_u)]

    return dX, dU
end

dX_kkt, dU_kkt = kkt_stuff(A,B,nx,nu,x0,N,Q,Qf,R)

function ilqr_stuff(A,B,Nx,Nu,x0,N,Q,Qf,R)
    """Solve for an iteration with iLQR instead of direct with KKT"""

    # allocate K and l
    K = cfill(Nu,Nx,N-1)
    l = cfill(Nu,N-1)

    # cost to go matrices at the end of the trajectory
    S = Qf
    s = Qf*(xtraj[N])

    # backwards pass
    for k = (N-1):-1:1

        # cost gradients at this time step
        q = Q*xtraj[k]
        r = R*utraj[k]

        # jacobians
        Ak = @views A[k]
        Bk = @views B[k]

        # linear solve
        l[k] = (R + Bk'*S*Bk)\(r + Bk'*s)
        K[k] = (R + Bk'*S*Bk)\(Bk'*S*Ak)

        # update
        Snew = Q + K[k]'*R*K[k] + (Ak-Bk*K[k])'*S*(Ak-Bk*K[k])
        snew = q - K[k]'*r + K[k]'*R*l[k] + (Ak-Bk*K[k])'*(s - S*Bk*l[k])

        # update S's
        S = copy(Snew)
        s = copy(snew)
    end

    # copy over old trajectory
    xnew = copy(xtraj)
    unew = copy(utraj)
    dX_ilqr = cfill(nx,N)
    dU_ilqr = cfill(nu,N-1)

    α = 1.0

    # forward rollout to get the correct dU_ilqr
    for k = 1:N-1

        # δx
        dX_ilqr[k] = xnew[k]-xtraj[k]

        # δu = -α*l - K*δx
        dU_ilqr[k] = -α*l[k] - K[k]*dX_ilqr[k]

        # unew = uold + δu
        unew[k] = utraj[k] + dU_ilqr[k]

        # rollout with xnew = xold + A*δx + B*δu
        xnew[k+1] = xtraj[k+1] + A[k]*dX_ilqr[k] + B[k]*dU_ilqr[k]
    end

    # δx for last knot point
    dX_ilqr[N] = xnew[N] - xtraj[N]

    return dX_ilqr, dU_ilqr

end

dX_ilqr, dU_ilqr = ilqr_stuff(A,B,nx,nu,x0,N,Q,Qf,R)


#control checking
u_error = [norm(dU_kkt[i] - dU_ilqr[i]) for i = 1:length(dU_ilqr)]
x_error = [norm(dX_kkt[i] - dX_ilqr[i]) for i = 1:length(dX_ilqr)]
mat"
figure
title('Control Error between KKT and iLQR')
hold on
plot($u_error)
hold off
%saveas(gcf,'ilqrvskkt_u.png')
"

mat"
figure
title('State Error between KKT and iLQR')
hold on
plot($x_error)
hold off
%saveas(gcf,'ilqrvskkt_x.png')
"

# NOTE: this commented out part is just to test the KKT solution if you want

# # now i'll verify KKT solution with convex.jl
# using Convex, Mosek, MosekTools
#
# function cvx_stuff()
# # delta x and delta u
# δx = Variable(nx,N)
# δu = Variable(nu,N-1)
#
# # initial condition constraint
# cons = Constraint[δx[:,1] == zeros(nx)]
#
# # linearized dynamics constraints
# for i = 1:N-1
#     push!(cons,δx[:,i+1] == A[i]*δx[:,i] + B[i]*δu[:,i])
# end
#
# # cost function
# J = 0
# for i = 1:N-1
#     J+= quadform((xtraj[i] + δx[:,i]),Q)
#     J+= quadform((utraj[i] + δu[:,i]),R)
# end
# J+= quadform((xtraj[N] + δx[:,N]),Qf)
#
# #solve
# problem = minimize(J, cons)
# solve!(problem, () -> Mosek.Optimizer(MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-15))
#
#
# u_cvx = vec_from_mat(evaluate(δu))
#
# return u_cvx
# end
# u_cvx = cvx_stuff()
# u_cvx_error = [norm(dU_kkt[i] - u_cvx[i]) for i = 1:length(l)]
#
# mat"
# figure
# title('Error between KKT and Mosek solution')
# hold on
# plot($u_cvx_error)
# hold off
# "
