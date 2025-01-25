function [x, iter, out] = gl_GD_primal(x0, A, b, mu, opts)
%smoothing technique + backward linesearch + BB step size + continuation
[m,n]=size(A);
l=size(b,2);
p=eigs(A'*A,1,'LM');
count=1;
u=10^-3;
v=1;
L=p+2*v/u;
x1=x0-grad_smoothedfun(x0,u,10^3*u)/p;
while u>1.1*10^-8   %or 1.1*10^-7 faster but with lower accuracy
    fopt=fun(x0,mu);
    s=x1-x0;
    y=grad_smoothedfun(x1,u,v)-grad_smoothedfun(x0,u,v);
    g=grad_smoothedfun(x1,u,v);
    a=norm(s,'fro')^2/sum(sum(s.*y));
    a=max(a,1/p);
    r=norm(g,'fro')^2;
    while a>1/L && smoothedfun(x1,u,v)-smoothedfun(x1-a*g,u,v)<1/(2*L)*r 
        a=a/2;
    end
    x0=x1;
    x1=x1-a*g;
    count=count+1;
    if norm(grad_smoothedfun(x1,u,v),'fro')<=10^-5
        u=u/10;
        v=max(mu,v/10);
        L=p+2*v/u;
    end
end
x=x0;
iter=count;
out=struct('fval',fopt);
    function y=smoothedfun(x,u,mu)%smoothed function
        h=sqrt(sum(x.^2,2));
        uu=u*ones(n,1);
        y=1/2*norm(A*x-b,'fro')^2+mu*sum(min(max(h-uu./2,uu./2),(h.^2./uu)/2));
    end
    function g=grad_smoothedfun(x,u,mu)%gradient of smoothed function
        h=repmat(sqrt(sum(x.^2,2)),1,l);
        uu=u*ones(n,1);
        g=A'*(A*x-b)+mu*x./max(uu,h);
    end
    function y=fun(x,mu)%original function
        y=1/2*norm(A*x-b,'fro')^2+mu*sum(sqrt(sum(x.^2,2)));      
    end

end
