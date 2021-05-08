function handle_fig=plot_2Yaxis(Y1, name_Y1, Y2, name_Y2)
% CHARLIE_raza_sweep case
% plot Y for specified X(size specified in this file: p×f)
p=[0.3, 0.5, 1, 3, 5, 10]'; 
handle_fig=figure;
if max(Y1(:))-min(Y1(:))<100 && max(Y1(:))/min(Y1(:))<100 ...
        && max(Y2(:))-min(Y2(:))<100 && max(Y2(:))/min(Y2(:))<100
    yyaxis left
    semilogx(p,Y1(:,1),'-.c');
    ylabel(name_Y1);
    axis([0.3,10,-inf,inf])
    yyaxis right
    semilogx(p,Y2(:,1),'--y');
    ylabel(name_Y2);
    axis([0.3,10,-inf,inf])
    hold on
    yyaxis left
    semilogx(p,Y1(:,2),'-.m');
    yyaxis right
    semilogx(p,Y2(:,2),'--g');
else
    yyaxis left
    loglog(p,Y1(:,1),'-.c');
    ylabel(name_Y1);
    axis([0.3,10,-inf,inf])
    yyaxis right
    loglog(p,Y2(:,1),'--y');
    ylabel(name_Y2);
    axis([0.3,10,-inf,inf])
    hold on
    yyaxis left
    loglog(p,Y1(:,2),'-.m');
    yyaxis right
    loglog(p,Y2(:,2),'--g');
end
xticks(p)
xlabel('{\itp} [Pa]');
L1=legend([name_Y1 ',f=1MHz'],[name_Y1 ',f=4MHz'],...
    [name_Y2 ',f=1MHz'],[name_Y2 ',f=4MHz']);
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格
end