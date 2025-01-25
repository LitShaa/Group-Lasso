function [x, iter, out] = gl_gurobi(x0, A, b, mu, opts)
[m,n]=size(A);
l=size(b,2);
% addpath('C:\gurobi1200\win64\matlab')
% savepath
gurobi_setup;
%variable u(i,j),x(k,l),s(i) omitted
F=zeros(m*l+n*(l+1));%linear term of objective function
model.Q=sparse((1:m*l),(1:m*l),1/2*ones(1,m*l),m*l+n*(1+l),m*l+n*(1+l));% quadratic term of objective function
model.lb=[-inf*ones((m+n)*l,1);zeros(n,1)];%lower bound
model.ub=inf*ones((m+n)*l+n,1);%upper bound
model.vtype=repmat('C',1,n*(l+1)+m*l);
model.obj=[zeros(1,(m+n)*l),mu*ones(1,n)];%objective function
model.modelsense='min';
%equivalent form of Ax-u=-b 
F(1:m*l,1:m*l)=-eye(m*l);
for j=1:l
    F(j:l:(m-1)*l+j,m*l+j:l:(m+n-1)*l+j)=A;
end
model.A=sparse(F);
model.rhs=[reshape(b',m*l,1);zeros(n*l+n,1)];
model.sense=repmat('=',1,n*(l+1)+m*l);
%quadratic constraints for each s(i)
for i=1:n
    model.quadcon(i).Qrow=[m*l+(i-1)*l+1:m*l+i*l,m*l+n*l+i];
    model.quadcon(i).Qcol=[m*l+(i-1)*l+1:m*l+i*l,m*l+n*l+i];
    model.quadcon(i).Qval=[ones(1,l),-1];
    model.quadcon(i).rhs=0;
    model.quadcon(i).q=sparse((m+n)*l+n,1);
    model.quadcon(i).sense='<';
end
params.feasibilityTol=1e-9;
params.OptimalityTol=1e-9;
result=gurobi(model,params);
if strcmp(result.status, 'OPTIMAL')
    disp('Optimal solution found.');
else
    disp(result);
end
x=reshape(result.x(m*l+1:(m+n)*l,1)',l,n)';
iter=result.baritercount;
out=struct('fval',result.objval);