function handle_fig=plot_parametric_1Ylog(Y1,name_Y1,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    loglog(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,end),'-g','LineWidth',1.5);
else
    semilogy(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,end),'-g','LineWidth',1.5);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']); %\it斜体 \rm 正体 ^上标 _下标 {}组合
ylabel(name_Y1);
L1=legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
set(L1,'FontSize',10);
title([name_Y1 ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end