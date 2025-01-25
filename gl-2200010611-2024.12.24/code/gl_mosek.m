function [x, iter, out] = gl_mosek(x0, A, b, mu, opts)
[m,n]=size(A);
l=size(b,2);
% addpath('C:\Users\沙柯岑\mosek\10.2\toolbox\R2017a');
% savepath;
[rcode,res]=mosekopt('symbcon echo(0)');
prob=[];
prob.a = sparse(0, n*l+n+1);
prob.c=[zeros(1,n*l),mu*ones(1,n),1/2];%variable:x(i,j),s(i),t
%write all constraints in the form Fx+g \in SOC
FQ=[];%corresponds to F
gQ=zeros(n*(l+1),1);%corresponds to g
cQ=[res.symbcon.MSK_DOMAIN_QUADRATIC_CONE l+1];%show that we're writing quadratic cones 
s=repmat(cQ,1,n);%there are n cones for each s(i)
for k=1:n
    %write FQ for constraints ||x(i,1:l)||_{2}<=s(i)
    FQ=[FQ;sparse([zeros(1,n*l+k-1),1,zeros(1,n+1-k);zeros(l,l*(k-1)),eye(l),zeros(l,n+1+(n-k)*l)])];
end
%write FQ for constraint:norm(Ax-b,'fro')^2<=t
%need to transform R-SOC into SOC
FFQ=zeros(2+m*l,n*l+n+1);
FFQ(1,:)=[zeros(1,n*l+n),1/2];
FFQ(2,:)=[zeros(1,n*l+n),1/2];
for i=1:m
    for j=1:l
        for k=1:n
            FFQ(2+l*(i-1)+j,(k-1)*l+j)=A(i,k);
        end
    end
end
ggQ=[1/2;-1/2;-reshape(b',m*l,1)];%corresponds to g
ccQ=[res.symbcon.MSK_DOMAIN_QUADRATIC_CONE m*l+2];%QC for the above constraint about t
%combine all constraints together
prob.f=[reshape(FQ,n*(l+1),n*l+n+1);FFQ];
prob.g=[gQ;ggQ];
prob.accs=[s ccQ];
[r,res]=mosekopt('minimize info',prob);
x=reshape(res.sol.itr.xx(1:n*l,1),l,n)';
iter=res.info.MSK_IINF_INTPNT_ITER;
out=struct('fval',res.sol.itr.pobjval);
