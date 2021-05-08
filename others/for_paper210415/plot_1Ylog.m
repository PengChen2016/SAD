function handle_fig=plot_1Ylog(Y, name_Y)
% CHARLIE_raza_sweep case
% plot Y for specified X(size specified in this file: p×f)
p=[0.3, 0.5, 1, 3, 5, 10]'; 
handle_fig=figure;
loglog(p,Y(:,1),'-.c');
hold on
loglog(p,Y(:,2),'-.m');
axis([0.3,10,-inf,inf])
xticks(p)
xlabel('{\itp} [Pa]');
ylabel(name_Y);
L1=legend('f=1MHz','f=4MHz');
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格
end