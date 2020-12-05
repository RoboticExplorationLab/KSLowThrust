using LinearAlgebra, SparseArrays, SuiteSparse
using Attitude

# function ctrb(A,B)
#     nx,nu = size(B)
#     C = zeros(nx,nu*nx)
#
#     for i = 1:nx
#         C[:,(i-1)*nu .+ (1:nu)] = A^(i-1)*B
#     end
#     if rank(C) != nx
#         error("not controllable")
#     end
# end
# nx = 5
# nu = 3
# A = randn(nx,nx)
# B = randn(nx,nu)
# ctrb(A,B)
As = [zeros(3,3) I;
       zeros(3,6)]
Bs = [zeros(3,3);I]
dt = 0.1
exp_discrete = exp([As Bs; zeros(3,9)]*dt)
A = exp_discrete[1:6,1:6]
B = exp_discrete[1:6,7:9]
nx,nu = size(B)

Q = randn(nx,nx);Q = sparse(Q'*Q);
Qf = 10*Q
# q = randn(nx)
R = randn(nu,nu);R = sparse(R'*R);
# r = randn(nu)

x0 = randn(nx)
N = 20

function kkt_stuff(A,B,nx,nu,x0,N,Q,Qf,R)
    nc = (N-1)*nx
    nz = (N)*nx + (N-1)*nu
    idx_x = [(t-1)*(nx+nu) .+ (1:nx) for t = 1:N]
    idx_u = [(t-1)*(nx+nu) .+ ((nx+1):(nx+nu)) for t = 1:N-1]
    idx_cons = [(t-1)*(nx) .+ (1:nx) for t = 1:(N-1)]
    C = spzeros(nc,nz)
    for i = 1:N-1
        C[idx_cons[i],idx_x[i]] = A
        C[idx_cons[i],idx_u[i]] = B
        C[idx_cons[i],idx_x[i+1]] = -I(nx)
    end
    C_lower = spzeros(nx,nz)
    C_lower[1:nx,1:nx] = sparse(I(nx))
    C = [C;C_lower]
    d = spzeros(nc + nx)
    d[(end-nx+1):end] = x0


    P = sparse(blockdiag(kron(sparse(I(N-1)),blockdiag(Q,R)),Qf))

    # KKT System
    kkt_A = [P C';
             C zeros(size(C,1),size(C,1))]

    kkt_b = [spzeros(nz);spzeros(nc);x0]


    z_sol = kkt_A\Array(kkt_b)

    X = [z_sol[idx_x[i]] for i = 1:length(idx_x)]
    U = [z_sol[idx_u[i]] for i = 1:length(idx_u)]
    λ = [z_sol[idx_cons[i] .+ (N*nx + (N-1)*nu)] for i = 1:length(idx_cons)]

    # @infiltrate
    # error()

    return X, U, λ
end

X_kkt, U_kkt, λ = kkt_stuff(A,B,nx,nu,x0,N,Q,Qf,R)


dyn_errors_kkt = [norm( A*X_kkt[i] + B*U_kkt[i] - X_kkt[i+1]) for i = 1:length(X_kkt)-1]


function ilqr_stuff(Ak,Bk,Nx,Nu,x0,N,Q,Qf,R)
    K = cfill(Nu,Nx,N-1)
    l = cfill(Nu,N-1)

    # allocate the new states and controls
    xtraj = cfill(Nx,N)
    utraj = cfill(Nu,N-1)

    # allocate gradients and hessians
    S = zeros(Nx,Nx)
    Snew = zeros(Nx,Nx)
    s = zeros(Nx)
    snew = zeros(Nx)


    # cost to go matrices at the end of the trajectory
    S .= Qf
    s .= zeros(Nx)
    q = zeros(Nx)
    r = zeros(Nu)

    # backwards pass
    for k = (N-1):-1:1

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
    xtraj[1] .= copy(x0)

    for i = 1:N-1
        utraj[i] = -K[i]*xtraj[i]
        # xtraj[i+1] = (A-B*K[i])*xtraj[i]
        xtraj[i+1] = A*xtraj[i] + B*utraj[i]
    end

    return xtraj, utraj

end


xtraj, utraj = ilqr_stuff(A,B,nx,nu,x0,N,Q,Qf,R)

x_error = [norm(xtraj[i]-X_kkt[i]) for i = 1:length(xtraj)]
u_error = [norm(utraj[i]-U_kkt[i]) for i = 1:length(utraj)]

# verify optimiality conditions
R*U_kkt[5] + B'*λ[5]


test_λ = cfill(nx,N-1)
test_λ[N-1] = Qf*xtraj[N]


# kkk = N-2
# kkk =5
# Q*xtraj[kkk] + A'*λ[kkk+1] - λ[kkk]


# @show norm(Q*xtraj[5] + A'*λ[6] - λ[5])
# @show norm(Q*xtraj[5] + A'*λ[6] - λ[5])
# @show norm(Q*xtraj[6] + A'*λ[6] - λ[5])
#
for i = N-2:-1:1
    test_λ[i] = Q*xtraj[i+1] + A'*test_λ[i+1]
end

λ_error = [norm(λ[i] - test_λ[i]) for i = 1:length(λ)]
mat"
figure
hold on
plot($λ_error)
hold off
"
mat"
figure
hold on
plot($x_error)
hold off
"
mat"
figure
hold on
plot($u_error)
hold off
"
