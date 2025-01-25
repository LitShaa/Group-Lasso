function [x, iter, out] = gl_ADMM_dual(x0, A, b, mu, opts)
%ADMM for dual problem with dynamic adjustment for penalty factor
[m,n]=size(A);
l=size(b,2);
iter=0;
p=zeros(n,l);
z1=ones(n,l);
z2=zeros(n,l);
x1=zeros(m,l);
Q=A*A';
sigma=100;          %penalty factor
u=10;a1=2;a2=2;     %parameters for adjusting penalty factor
vecmu=mu*ones(n,l);
R=chol(eye(m)+sigma*Q);
vecsigma=sigma*ones(n,l);
while norm(A*(z1-z2),'fro')>10^-10
    iter=iter+1;
    z1=z2;
    x1=R'\(b-A*p+sigma*(A*z1));
    x1=R\x1;
    r=p+sigma*(A'*x1);
    s=repmat(sqrt(sum(r.^2,2)),1,l);
    z2=r./s.*min(vecmu,s./vecsigma);
    p=p+sigma*(A'*x1-z2);
    if norm(z2-A'*x1,'fro')>u*norm(A*(z1-z2),'fro')
        sigma=sigma*a1; R=chol(eye(m)+sigma*Q);vecsigma=a1*vecsigma;
    elseif norm(A*(z1-z2),'fro')>u*norm(z2-A'*x1,'fro')
        sigma=sigma/a2; R=chol(eye(m)+sigma*Q);vecsigma=vecsigma/a2;
    end
end
x=p;
fopt=1/2*norm(A*x-b,'fro')^2+mu*sum(sqrt(sum(x.^2,2)));
out=struct('fval',fopt);
end