function [x, iter, out] = gl_FGD_primal(x0, A, b, mu, opts7)
%光滑化+Nesterov加速+连续化技巧
[m,n]=size(A);
l=size(b,2);
p=eigs(A'*A,1,'LM');
count=1;
k=1;
u=10^-3;
v=1;
L=p+2*v/u;
x1=x0-grad_smoothedfun(x0,u,v)/p;
while u>10^-9       %也可以更小,求得的解更精确
    count=count+1;
    k=k+1;
    fopt=fun(x0,mu);
    y=x1+(k-2)/(k+1)*(x1-x0);
    x0=x1;
    x1=y-1/L*grad_smoothedfun(y,u,v);
    if norm(grad_smoothedfun(x1,u,v),'fro')<=10^-5
        u=u/10;
        v=max(mu,v/10);
        L=p+2*v/u;
        k=1;
    end
end
x=x0;
iter=count;
out=struct('fval',fopt);
    function g=grad_smoothedfun(x,u,mu)%光滑化函数梯度
        h=repmat(sqrt(sum(x.^2,2)),1,2);
        uu=u*ones(n,l);
        g=A'*(A*x-b)+mu*x./max(uu,h);
    end
    function y=fun(x,mu)%原函数
        y=1/2*norm(A*x-b,'fro')^2+mu*sum(sqrt(sum(x.^2,2)));      
    end
end