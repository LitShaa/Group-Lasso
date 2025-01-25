function [x, iter, out] = gl_ProxGD_primal(x0, A, b, mu, opts)
%近似点梯度法+bb步长选取+连续化技巧 %proximal+bb stepsize+continuation
[m,n]=size(A);
l=size(b,2);
p=eigs(A'*A,1,'LM');
c=p;
count=1;
x1=prox(x0-1/c*(A'*(A*x0-b)),1/c*mu);
v=mu*5^5;
while v>=mu 
      fopt=fun(x0,mu);
      s=x0-x1;
      if norm(s,'fro')<10^-5 && v>mu %continuation
          v=v/5;
      end
      a=max((norm(s,'fro')/norm(A*s,'fro'))^2,1/c);%bb stepsize
      x0=x1;
      g=A'*(A*x1-b);
      while sqrt(a)*norm(A*g,'fro')>norm(g,'fro') % backward linesearch with armoji criteria
          a=a/2;
      end
      x1=prox(x1-a*g,a*v);
      count=count+1;
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
function y=fun(x,mu) %original function
        y=1/2*norm(A*x-b,'fro')^2+mu*sum(sqrt(sum(x.^2,2)));      
end
end