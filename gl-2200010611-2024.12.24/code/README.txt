Matlab使用的是R2023a版本
cvx使用的是2.2.2版本
mosek使用的是10.2.7版本
gurobi使用的是12.0.0版本
所有程序均使用MATLAB自带工具箱

以随机种子 24847563（已经设置成默认种子）运行Test_group_lasso.m可以得到报告中的结果
报告中最后的可视化图可通过运行Test_group_lasso.m的plot_results部分得到
在本地运行时，gl_mosek.m和gl_gurobi.m程序中有添加路径这一步骤，提交时已将其注释化