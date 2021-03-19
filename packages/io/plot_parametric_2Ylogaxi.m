function handle_fig=plot_parametric_2Ylogaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1与Y2双纵轴随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    yyaxis left;
    loglog(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    loglog(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
else
    yyaxis left;
    semilogy(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    semilogy(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
set(L1,'FontSize',10);
title([name_Y ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end