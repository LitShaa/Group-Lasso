function [x, iter, out] = gl_ADMM_primal(x0, A, b, mu, opts)
[m,n]=size(A);
L=eigs(A'*A,1,'LM');
c=A'*b;
l=size(b,2);    
x1=x0+ones(n,l);
x2=x0;
p=randn(m,l);
veczero=zeros(n,1);
vecmu=mu*ones(n,1);
iter=0;             %num of total iterations
sigma=0.0001;       %initial penalty factor
a1=10;a2=10;        %adjust rate of penalty factor
u=10;               %coordinate ratio of KKT condition tolerance
while norm(A*(x1-x2),'fro')>10^-10
    r=norm((x1-x2),'fro')^2/sigma/norm(A*(x1-x2),'fro')^2; %BB step size
    x1=x2;iter=iter+1;
    s=sigma*(A*x2-b+p)/(1+sigma); %update first variable
    %linesearch
    while  fun(x2,s,p,sigma)-fun(x2-r*gradfun(x2,s,p,sigma),s,p,sigma)<r/2*norm(gradfun(x2,s,p,sigma),'fro')^2
        r=r/2;
    end
    r=max(r,1/L/sigma);
    %update second variable using linearization
    g=x2-r*sigma*(A'*(A*x2-b-s+p/sigma));
    x2=repmat(max(veczero,(1-vecmu./sqrt(sum(g.^2,2))*r)),1,l).*g;
    p=p+sigma*(A*x2-b-s); %update multiplers
    %adjust penalty function
    if u*norm(A*(x1-x2),'fro')<norm(A*x2-b-s,'fro')
        sigma=sigma*a1;
    elseif u*norm(A*x2-b-s,'fro')<norm(A*(x1-x2),'fro')
        sigma=sigma\a2;
    end
end
x=x1;
fopt=1/2*norm(A*x-b,'fro')^2+mu*sum(sqrt(sum(x.^2,2)));
out=struct('fval',fopt);
    function a=fun(x,s,p,sigma)
        a=sigma/2*norm(A*x-b-s+p/sigma,'fro')^2;
    end
    function g=gradfun(x,s,p,sigma) 
        g=sigma*(A'*(A*x-b-s+p/sigma));
    end
end


