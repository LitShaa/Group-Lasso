function [x, iter, out] = gl_SGD_primal(x0, A, b, mu, opts)
%backward linesearch + bb step size + continuation
n=size(A,2);
l=size(b,2);
p=eigs(A'*A,1,'LM');
v=mu*5^5;
c=p;                   %Lipchitz constant
count=0;
epis=(10^-6)*ones(n,1);%regard x as 0 when ||x||<=10^-6
x1=x0-subgrad_fun(x0,v)/c;
opt=fun(x0,mu);
while  v>mu || abs(fun(x0,mu)-fun(x1,mu))>10^-10 
   s=x1-x0;
   y=subgrad_fun(x1,v)-subgrad_fun(x0,v);
   g=subgrad_fun(x1,v);
   a=norm(s,'fro')^2/sum(sum(s.*y));%bb step size
   a=max(a,1/c);
   r=norm(g,'fro')^2;
   %line search
   while fun(x1,v)-fun(x1-a*g,v)<1/(2*c)*r && a>1/c
        a=a/2;
   end
   while fun(x1,v)-fun(x1-a*g,v)<0
       a=a/2;
   end
        x0=x1;
        x1=x1-a*g;
        count=count+1;
        opt=fun(x0,mu);
        %continuation
        if fun(x0,v)-fun(x1,v)<=10^-5*v && v>mu
            v=v/5;
        end
end
x=x0;
iter=count;
out=struct('fval',opt);
    function y=fun(x,mu)
        y=1/2*norm(A*x-b,'fro')^2+mu*sum(sqrt(sum(x.^2,2)));
    end
    function g=subgrad_fun(x,mu)%subgradient of fun
        h=max(epis,sqrt(sum(x.^2,2)));
        g=A'*(A*x-b)+mu*x./repmat(h,1,l);
    end
end