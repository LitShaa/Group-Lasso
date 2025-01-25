function [x, iter, out] = gl_FProxGD_primal(x0, A, b, mu, opts)
%FISTA, fixed stepsize+ continuation
[m,n]=size(A);
l=size(b,2);
p=eigs(A'*A,1,'LM');
c=p;
count=1;
k=1;
v=mu*10^5;

x1=prox(x0-1/c*(A'*(A*x0-b)),1/c*v);
while v>=mu
    fopt=fun(x0);
    count=count+1;
    k=k+1;
    y=x1+(k-2)/(k+1)*(x1-x0); %FISTA
    x0=x1;
    x1=prox(y-1/c*(A'*(A*y-b)),1/c*v);
    if norm(x0-x1,'fro')<10^-5 && v>mu %continuation
        v=v/10;k=1;
    end
    if v==mu && norm(x0-x1,'fro')<10^-8
        break
    end
end
x=x0;
iter=count;
out=struct('fval',fopt);
function prox=prox(x,t)%proximal operator
    prox=repmat(max(zeros(n,1),ones(n,1)-t./sqrt(sum(x.^2,2))),1,l).*x;
end
    function y=fun(x)%original function
        y=1/2*norm(A*x-b,'fro')^2+mu*sum(sqrt(sum(x.^2,2)));      
end
end
 
