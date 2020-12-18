using ForwardDiff, StaticArrays
using LinearAlgebra
const FD = ForwardDiff
# const MTK = ModelingToolkit


function dynamics(x,u)
    """ Here is the state
            p  ∈ R⁴ [1:4]
            p′ ∈ R⁴ [5:8]
            h  ∈ R  [9]
           |r| ∈ R  [10]
            r₃ ∈ R  [11]
            u  ∈ R³ [12:14]
            u̇  ∈ R³ [15:17]
            t  ∈ R  [18]
    """
    return params.invJ*(tanh.(u)*x[1] -(x×(params.J*x)))
end

function rk4(x_in)

    # x_n = SVector{6}(x_n)
    # x_n = @views x[1:14]
    # u = @views x[15:17]
    # t = @views x[18]
    # dt = @views x[19]
    x_n = x_in[1:3]
    u = x_in[4:6]
    dt = x_in[7]
    # u = 0
    # t = 0
    k1 = dt*dynamics(x_n,u)
    k2 = dt*dynamics(x_n+k1/2,u)
    k3 = dt*dynamics(x_n+k2/2,u)
    k4 = dt*dynamics(x_n+k3,u)

    return (x_n + (1/6)*(k1+2*k2+2*k3+k4))

end
function bdynamics(xt,ut,tt)
    _Afx(x) = dynamics(x,ut)
    _Bfx(u) = dynamics(xt,u)
    return dynamics(xt,ut), FD.jacobian(_Afx,xt), FD.jacobian(_Bfx,ut)
    # return dynamics(xt), FD.jacobian(dynamics,xt), zeros(3,3)
end
function rk4step_jacobians(x1,u0,dt)
    t = 0
    # u0 = 0
    xdot1, A1, B1 = bdynamics(x1,u0,t)
    k1 = xdot1*dt

    x2 = x1 + .5*k1
    xdot2, A2, B2 = bdynamics(x2,u0,t+dt*.5)
    k2 = dt*xdot2;


    x3 = x1 + .5*k2
    xdot3, A3, B3 = bdynamics(x3,u0,t+dt*.5)
    k3 = dt*xdot3

    x4 = x1 + k3
    xdot4, A4, B4 = bdynamics(x4,u0,t+dt)
    k4 = dt*xdot4

    x_tp1 = x1 + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

    # # RK4 method
    # dk1_dx1 = dt*A1
    # # dx2_dx1 = I + .5*dk1_dx1
    # dk2_dx1 = dt*A2*(I + .5*dk1_dx1)
    # # dx3_dx1 = I + .5*dk2_dx1
    # dk3_dx1 = dt*A3*(I + .5*dk2_dx1)
    # # dx4_dx1 = I + dk3_dx1
    # dk4_dx1 = dt*A4*(I + dk3_dx1)
    # A_d = I + (1/6)*(dk1_dx1 + 2*dk2_dx1 + 2*dk3_dx1 + dk4_dx1)
    #
    # # B_d
    # dk1_du = dt*B1
    # # dx2_du = .5*dk1_du
    # dk2_du = dt*A2*(.5*dk1_du) + dt*B2
    # # dx3_du = .5*dk2_du
    # dk3_du = dt*A3*(.5*dk2_du) + dt*B3
    # # dx4_du = dk3_du
    # dk4_du = dt*A4*dk3_du + dt*B4
    # B_d = (1/6)*(dk1_du + 2*dk2_du + 2*dk3_du + dk4_du)
    # RK4 method
    dk1_dx1 = dt*A1
    dx2_dx1 = I + .5*dk1_dx1
    dk2_dx1 = dt*A2*dx2_dx1
    dx3_dx1 = I + .5*dk2_dx1
    dk3_dx1 = dt*A3*dx3_dx1
    dx4_dx1 = I + dk3_dx1
    dk4_dx1 = dt*A4*dx4_dx1
    A_d = I + (1/6)*(dk1_dx1 + 2*dk2_dx1 + 2*dk3_dx1 + dk4_dx1)

    # B_d
    dk1_du = dt*B1
    dx2_du = .5*dk1_du
    dk2_du = dt*A2*dx2_du + dt*B2
    dx3_du = .5*dk2_du
    dk3_du = dt*A3*dx3_du + dt*B3
    dx4_du = dk3_du
    dk4_du = dt*A4*dx4_du + dt*B4
    B_d = (1/6)*(dk1_du + 2*dk2_du + 2*dk3_du + dk4_du)

    return A_d, B_d, x_tp1
    # return A_d, x_tp1
end

function testall()
xt = SVector{3}([.2;7;.2])
ut = SVector{3}(randn(3))
# tt = randn()
dtt = 0.12
J = Diagonal(SVector{3}([1;2;3]))
global params = (J = J, invJ = inv(J))
xtp1_1 = rk4([xt;ut;dtt])


A_d, B_d, x_tp1 = rk4step_jacobians(xt,ut,dtt)
@info "first time"
@btime A_d, B_d, x_tp1 = rk4step_jacobians($xt,$ut,$dtt)
@show norm(xtp1_1 - x_tp1)

# @time AA = FD.jacobian(rk4,[xt;ut;dtt])
# @time AA = FD.jacobian(rk4,[xt;ut;dtt])
AA = FD.jacobian(rk4,[xt;ut;dtt])
@info "second time"
@btime AA = FD.jacobian(rk4,[$xt;$ut;$dtt])
A_d2 = AA[:,1:3]
B_d2 = AA[:,4:6]

@show norm(A_d - A_d2)
@show norm(B_d - B_d2)
@show B_d

end

testall()
