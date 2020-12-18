clear


syms x1 x2 x3 lambda1 lambda2 lambda3 lambda4 mu real 


x = [x1;x2;x3];
lambda = [lambda1; lambda2; lambda3; lambda4];

L = Lag(x,0,lambda,mu)

simplify(gradient(L,x))
simplify(hessian(L,x))


function out = c_fx(x,u)
    % u_max = 1.0
    syms u_max real
    out = [x(1:3);u_max];
end
function out = proj(x)
    v = x(1:3);
    s = x(4);
    nv = norm(v);

  
%         out = zeros(4,1);
%         out = x;
        out = (.5*(1 + s/nv)*[v;nv]);
end

function out = Lag(x,u,lambda,mu)

    
        out =  (1/(2*mu))*(  norm(proj(lambda - mu*c_fx(x,u)))^2 - dot(lambda,lambda));
%     # return LQR_cost(Q,R,x,u,xf)
end
