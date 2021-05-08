function handle_fig=plot_2Ylog(Y1, name_Y1, Y2, name_Y2, name_Y)
% CHARLIE_raza_sweep case
% plot Y for specified X(size specified in this file: p×f)
p=[0.3, 0.5, 1, 3, 5, 10]';
handle_fig=figure;
if max(Y1(:))-min(Y1(:))<100 && max(Y1(:))/min(Y1(:))<100 ...
        && max(Y2(:))-min(Y2(:))<100 && max(Y2(:))/min(Y2(:))<100
    semilogx(p,Y1(:,1),'-.c');
    axis([0.3,10,-inf,inf])
    hold on
    semilogx(p,Y1(:,2),'-.m');
    semilogx(p,Y2(:,1),'--y');
    semilogx(p,Y2(:,2),'--g');
else
    loglog(p,Y1(:,1),'-.c');
    axis([0.3,10,-inf,inf])
    hold on
    loglog(p,Y1(:,2),'-.m');
    loglog(p,Y2(:,1),'--y');
    loglog(p,Y2(:,2),'--g');
end
ylabel(name_Y);
xticks(p)
xlabel('{\itp} [Pa]');
L1=legend([name_Y1 ',f=1MHz'],[name_Y1 ',f=4MHz'],...
    [name_Y2 ',f=1MHz'],[name_Y2 ',f=4MHz']);
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格
end