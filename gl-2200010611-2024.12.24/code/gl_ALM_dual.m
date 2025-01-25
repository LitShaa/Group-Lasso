function [t, iter, out] = gl_ALM_dual(x0, A, b, mu, opts)
%Augmented Lagrangian method for dual problem
[m,n]=size(A);
l=size(b,2);
err0=10^-8;     %violation tolerance
accy0=10^-5;    %accuracy tolerance
sigma=100000;   %initial penalty factor
rho=10;         %growth rate of penalty factor
accy=10^-2;     %initial accuracy
err=1;          %initial violation tolerance
count=0;        %num of outer cycle
iter=0;         %num of total iterations
p=zeros(n,1);
x=zeros(m,l);
while vio(x,p,sigma)>err0 || norm(grad_AL_fun(x,p,sigma),'fro')>accy0
    if sigma>10^7
        break
    end
    count=count+1;
    x=argmin_AL_fun(x,p,sigma,accy);
    if vio(x,p,sigma)<=err
        if vio(x,p,sigma)<=err0 && norm(grad_AL_fun(x,p,sigma),'fro')<=accy0
            break
        end
        p=update_fun(x,p,sigma);
        err=err/rho;
    else
        sigma=rho*sigma;
        accy=accy/5;
    end
end
Q=A'*x;
t=zeros(n,l);
for i=1:n
    t(i,:)=Q(i,:)/norm(Q(i,:),2)*2*p(i)*mu;
end
fopt=fun(t);
out=struct('fval',fopt);
    function y=fun(x)%original function
        y=1/2*norm(A*x-b,'fro')^2+mu*sum(sqrt(sum(x.^2,2)));
    end
    function x1=argmin_AL_fun(x,p,sigma,accy)%minimize AL function
        y1=x;y0=x; k=1 ;a=1;
        g=grad_AL_fun(y1,p,sigma);
        while norm(g,'fro')>accy && a>10^-8
            z=y1+(k-1)/(k+2)*(y1-y0);
            g=grad_AL_fun(z,p,sigma);
            while AL_fun(z,p,sigma)-AL_fun(z-a*g,p,sigma)<a/2*norm(g,'fro')^2
                if a<10^-8 
                    break
                end
                a=a/2;              
            end
            y0=y1;y1=z-a*g;
            k=k+1;
            iter=iter+1;
            g=grad_AL_fun(y1,p,sigma);
        end
        x1=y1;
    end
    function output = update_fun(x, p, sigma)%multiplers update function
        output = max(zeros(n,1), sigma * (sum((A'*x).^2,2) - mu^2) + p);
    end
    function output = dual_fun(x)%objective function of dual problem
        output = 1/2*norm(x,'fro')^2 - sum(x(:).*b(:));
    end
    function output = AL_fun(x, p, sigma)%Augmented Lagrangian function
        output = dual_fun(x) + 1/2*(sum(update_fun(x, p, sigma).^2) - sum(p.^2))/sigma;
    end
    function output = vio(x, p, sigma)%violation tolerance
        output = norm((update_fun(x, p, sigma) - p),2)/sigma;
    end
    function output = grad_AL_fun(x, p, sigma)%gradient of AL function
        output = x - b + 2*(A*((repmat(update_fun(x,p,sigma),1,l).*(A'*x))));
    end
end









