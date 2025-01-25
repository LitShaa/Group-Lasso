function plot_results(v, titleText, savePath, u, x1,x2)
% plot the solution of each method
% This function assumes l=2
[n,l]=size(v);
fig = figure(1);
subplot(2,2,1); 
plot(1:n, v(:,1), '*', 1:n, v(:,2), 'o');
title(titleText)
xlim([1 n])
subplot(2,2,2);
plot(1:n, v(:,1)-u(:,1), '*', 1:n, v(:,2)-u(:,2), 'o');
xlim([1 n])
title('Error to exact');
subplot(2,2,3);
plot(1:n, v(:,1)-x1(:,1), '*', 1:n, v(:,2)-x1(:,2), 'o');
xlim([1 n])
title('Error to cvx-mosek');
subplot(2,2,4);
plot(1:n, v(:,1)-x2(:,1), '*', 1:n, v(:,2)-x2(:,2), 'o');
xlim([1 n])
title('Error to cvx-gurobi');
saveas(fig,savePath);
end