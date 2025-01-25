function [x, iter, out] =gl_cvx_gurobi(x0, A, b, mu, opts)
n=size(A,2);
l=size(b,2);
cvx_begin
    cvx_solver gurobi
    variable t(n,l);
    minimize (1/2*square_pos(norm(A*t-b,"fro"))+mu*sum(norms(t,2,2)));
cvx_end
x=t;
iter=-1;
out=struct('fval',cvx_optval);

end