using LinearAlgebra

function dyn(t,x,u)
    J = Diagonal([1;2;3])
    ω = x[1:3]
    return J\(u - ω × (J*ω))
end


function  RK8(func, t, y,u, h )

    k_1 = func(t         ,y            ,u)
    k_2 = func(t+h*(4/27),y+(h*4/27)*k_1                        ,u)
    k_3 = func(t+h*(2/9) ,y+  (h/18)*(k_1+3*k_2)    ,u)
    k_4 = func(t+h*(1/3) ,y+  (h/12)*(k_1+3*k_3)               ,u)
    k_5 = func(t+h*(1/2) ,y+   (h/8)*(k_1+3*k_4)         ,u)
    k_6 = func(t+h*(2/3) ,y+  (h/54)*(13*k_1-27*k_3+42*k_4+8*k_5) ,u)
    k_7 = func(t+h*(1/6) ,y+(h/4320)*(389*k_1-54*k_3+966*k_4-824*k_5+243*k_6),u)
    k_8 = func(t+h       ,y+  (h/20)*(-234*k_1+81*k_3-1164*k_4+656*k_5-
    122*k_6+800*k_7)               ,u)
    k_9 = func(t+h*(5/6) ,y+ (h/288)*(-127*k_1+18*k_3-678*k_4+456*k_5-9*k_6+
    576*k_7+4*k_8)            ,u)
    k_10= func(t+h       ,y+(h/820)*(1481*k_1-81*k_3+7104*k_4-3376*k_5+72*k_6
    -5040*k_7-60*k_8+720*k_9),u)

    y = y + h/840*(41*k_1+27*k_4+272*k_5+27*k_6+216*k_7+216*k_9+41*k_10)

    return y

end

function testrk8()
    x0 = [.2;7;.2]

    N = 100
    dt = .1
    X = cfill(3,N)
    X[1] = x0
    for i = 1:N-1
        X[i+1] = RK8(dyn,0,X[i],zeros(3),dt)
    end


    Xm = mat_from_vec(X)

    mat"
    figure
    hold on
    plot($Xm')
    hold off
    "
end

testrk8()
