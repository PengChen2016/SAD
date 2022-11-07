% BUG fig for paper
% 说明有仿真聚变负源的能力，显示其优点
clear
close all

path='d:\School\DoctorProgram\eP-项目笔记文件夹\eP-190821-01激励器FEM模型\paper_v2-210427期刊论文\绘图\';
% ft.print_fed('vec_emf',[path 'Fig4'])

%% input
load('BUG_results.mat')
load('paper_BUG_edge_t210517.mat')
plasma=input.plasma;
source_t=source;
load('paper_BUG_edge_a210517.mat')
source_a=source;
plasma_a=input_m.plasma;
load('paper_BUG_edge_tj210517.mat')
plasma_j=input.plasma;
source_j=source;

constants=get_constants();
w1MHz=2*pi*1e6;
co=my.get_color_order('');
ft=my.figtool;
ft.init_fed();

d1=@(mat,ip) reshape(mat(ip,1,:),1,3);
d2=@(mat,iPin) reshape(mat(:,1,iPin),1,3);
o=@(mat) mat([1,4,7,2,5,8,3,6,9]');
p=[0.2,0.3,0.5]';
P_RF=[18,36,55];

%% 选取代表性参数组
% % PER-p
% figure
% % experiment
% subplot(1,2,1);
% legend_text={};
% plot(p,d2(experiment.PER,2),'-o')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,d2(experiment.PER,3),'-.o')
% legend_text{end+1}='{\itP}_{S}=55kW';
% plot(p,d2(experiment.PER,1),'--o')
% legend_text{end+1}='{\itP}_{S}=18kW';
% ylabel('PER [\Omega]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% L1=legend(legend_text);
% set(L1,'location','best');
% set(L1,'box','off')
% % model
% subplot(1,2,2);
% legend_text={};
% plot(p,d2(fem.PER,2),'-s')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,d2(fem.PER,3),'-.s')
% legend_text{end+1}='{\itP}_{S}=55kW';
% plot(p,d2(fem.PER,1),'--s')
% legend_text{end+1}='{\itP}_{S}=18kW';
% ylabel('PER [\Omega]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% L1=legend(legend_text);
% set(L1,'location','best');
% set(L1,'box','off')
% 
% % PER-P_RF
% figure
% % experiment
% subplot(1,2,1);
% legend_text={};
% plot(P_RF,d1(experiment.PER,2),'-o')
% legend_text{end+1}='p=0.3Pa';
% hold on
% plot(P_RF,d1(experiment.PER,3),'-.o')
% legend_text{end+1}='p=0.5Pa';
% plot(P_RF,d1(experiment.PER,1),'--o')
% legend_text{end+1}='p=0.2Pa';
% ylabel('PER [\Omega]')
% xlabel('{\itP}_{S} [kW]')
% xticks(P_RF)
% grid on
% L1=legend(legend_text);
% set(L1,'location','best');
% set(L1,'box','off')
% % model
% subplot(1,2,2);
% legend_text={};
% plot(P_RF,d2(fem.PER,2),'-s')
% legend_text{end+1}='p=0.3Pa';
% hold on
% plot(P_RF,d2(fem.PER,3),'-.s')
% legend_text{end+1}='p=0.5Pa';
% plot(P_RF,d2(fem.PER,1),'--s')
% legend_text{end+1}='p=0.2Pa';
% ylabel('PER [\Omega]')
% xlabel('{\itP}_{S} [kW]')
% xticks(P_RF)
% grid on
% L1=legend(legend_text);
% set(L1,'location','best');
% set(L1,'box','off')
% 
% % 决定：PER-p, 36/55kW两条线
y36kW=@(mat) d2(mat,2);
y55kW=@(mat) d2(mat,3);

y18kW=@(mat) d2(mat,1);

%% 等离子体模型
% ne, Te
figure
% ne
subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
legend_text={};
temp=co(1,:);
plot(p,y36kW(plasma.ne),'-o','Color',temp)
legend_text{end+1}='{\itP}_{S}=36kW';
hold on
plot(p,y55kW(plasma.ne),'-.o','Color',temp)
legend_text{end+1}='{\itP}_{S}=55kW';
ylabel('{\itn}_e [m^{-3}]')
xlabel('{\itp}_{fill} [Pa]')
xticks(p)
grid on
L1=legend(legend_text);
set(L1,'location','south');
set(L1,'box','off')
max_y=6e17;
axis([0.2,0.5,0,max_y]) 
text(0.25,0.95*max_y,'(a)')
% Te
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
legend_text={};
temp=co(2,:);
plot(p,y36kW(plasma.Te),'-s','Color',temp)
legend_text{end+1}='{\itP}_{S}=36kW';
hold on
plot(p,y55kW(plasma.Te),'-.s','Color',temp)
legend_text{end+1}='{\itP}_{S}=55kW';
ylabel('{\itT}_e [eV]')
xlabel('{\itp}_{fill} [Pa]')
xticks(p)
grid on
L1=legend(legend_text);
set(L1,'location','south');
set(L1,'box','off')
max_y=10;
axis([0.2,0.5,0,max_y]) 
text(0.25,0.95*max_y,'(b)')

% 手动调整图例位置
fig_name='Fig13';
ft.print_fed('vec_emf', [path fig_name])
ft.print_fed('vec_eps', [path fig_name])


% % frequency
% figure
% % 36kW
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% plot(p,y36kW(plasma.nu_m),'-o')
% legend_text{end+1}='{\it\nu}_{m}';
% hold on
% plot(p,y36kW(plasma.nu_st),'-s')
% legend_text{end+1}='{\it\nu}_{st}';
% ylabel('{frequency [s^{-1}]}')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% line([0.2,0.5],[w1MHz, w1MHz],'linestyle','--','color','k');
% legend_text{end+1}='{\it\omega}_{RF}';
% L1=legend(legend_text);
% set(L1,'location','east');
% set(L1,'box','off')
% max_y=2e7;
% axis([0.2,0.5,0,max_y]) 
% text(0.25,0.95*max_y,'(a) {\itP}_{S}=36kW')
% % 55kW
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% plot(p,y55kW(plasma.nu_m),'-.o')
% legend_text{end+1}='{\it\nu}_{m}';
% hold on
% plot(p,y55kW(plasma.nu_st),'-.s')
% legend_text{end+1}='{\it\nu}_{st}';
% ylabel('{frequency [s^{-1}]}')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% line([0.2,0.5],[w1MHz, w1MHz],'linestyle','--','color','k');
% legend_text{end+1}='{\it\omega}_{RF}';
% L1=legend(legend_text);
% set(L1,'location','east');
% set(L1,'box','off')
% max_y=2e7;
% axis([0.2,0.5,0,max_y]) 
% text(0.25,0.95*max_y,'(b) {\itP}_{S}=55kW')

% % sigma, eps
% figure
% % sigma
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% plot(p,y36kW(plasma.sigma),'-o')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,y55kW(plasma.sigma),'-.o')
% legend_text{end+1}='{\itP}_{S}=55kW';
% ylabel('{\it\sigma} [S/m]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% L1=legend(legend_text);
% set(L1,'location','south');
% set(L1,'box','off')
% max_y=600;
% axis([0.2,0.5,0,max_y]) 
% text(0.25,0.95*max_y,'(a)')
% % -omega*eps
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% plot(p,y36kW(-w1MHz*plasma.eps_prime),'-s')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,y55kW(-w1MHz*plasma.eps_prime),'-.s')
% legend_text{end+1}='{\itP}_{S}=55kW';
% ylabel('-{\it\omega\epsilon} [S/m]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% L1=legend(legend_text);
% set(L1,'location','northeast');
% set(L1,'box','off')
% max_y=600;
% axis([0.2,0.5,0,max_y]) 
% text(0.25,0.95*max_y,'(b)')
% % eps_r
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% plot(p,y36kW(plasma.eps_r),'-s')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,y55kW(plasma.eps_r),'-.s')
% legend_text{end+1}='{\itP}_{S}=55kW';
% ylabel('{\it\epsilon}_r')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% L1=legend(legend_text);
% set(L1,'location','south');
% set(L1,'box','off')

delta_m_geo=get_plasma_skin_depth('as-medium-simplified-finite-radius',...
    plasma.f, plasma.nu_eff, plasma.wpe, plasma.r);
r_plasma=1e3*plasma.r;
% % skindepth, wavelength
% figure
% % skindepth
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% plot(p,1e3*y36kW(plasma.skin_depth),'-o')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,1e3*y55kW(plasma.skin_depth),'-.o')
% legend_text{end+1}='{\itP}_{S}=55kW';
% % plot(p,1e3*y55kW(delta_m_geo),'--d')
% % legend_text{end+1}='{\itP}_{S}=55kW, finite-radius';
% ylabel('{\it\delta} [mm]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% line([0.2,0.5],[r_plasma, r_plasma],'linestyle','--','color','k');
% legend_text{end+1}='{\itr}_{plasma}';
% L1=legend(legend_text);
% set(L1,'location','south');
% set(L1,'box','off')
% max_y=r_plasma;
% axis([0.2,0.5,0,max_y]) 
% text(0.25,0.95*max_y,'(a)')
% % wavelength
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% plot(p,1e3*y36kW(plasma.wavelength),'-s')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,1e3*y55kW(plasma.wavelength),'-.s')
% legend_text{end+1}='{\itP}_{S}=55kW';
% ylabel('\lambda [mm]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% line([0.2,0.5],[r_plasma, r_plasma],'linestyle','--','color','k');
% legend_text{end+1}='{\itr}_{plasma}';
% L1=legend(legend_text);
% set(L1,'location','south');
% set(L1,'box','off')
% max_y=240;
% axis([0.2,0.5,r_plasma,max_y]) 
% text(0.25,0.95*(max_y-r_plasma)+r_plasma,'(b)')




%% basic
figure
set(gca, 'Position',[0.25+0.12 0.19 0.355 0.73])
legend_text={};
plot(nan,nan,'-k')
legend_text{end+1}='{\itP}_{S}=36kW';
hold on
plot(nan,nan,'-.k')
legend_text{end+1}='{\itP}_{S}=55kW';
set(gca,'ColorOrderIndex',1);
plot(p,y36kW(experiment.PER),'-o')
plot(p,y55kW(experiment.PER),'-.o')
text(0.22,2,'experiment')
set(gca,'ColorOrderIndex',1);
plot(p,y36kW(fem.PER),'-s')
plot(p,y55kW(fem.PER),'-.s')
text(0.22,7.5,'FEM')
set(gca,'ColorOrderIndex',1);
plot(p,y36kW(source_t.PER),'-d')
plot(p,y55kW(source_t.PER),'-.d')
text(0.22,14,'transformer')
% plot(p,d2(fem.PER,1),'--s')
% legend_text{end+1}='{\itP}_{S}=18kW';
% plot(p,y36kW(source_j.PER),'-x')
% plot(p,y55kW(source_j.PER),'-.x')
ylabel('PER [\Omega]')
xlabel('{\itp}_{fill} [Pa]')
xticks(p)
grid on
L1=legend(legend_text);
set(L1,'location','west');
set(L1,'box','off')
min_y=0;
max_y=17;
axis([0.2,0.5,min_y,max_y]) 


figure
set(gca, 'Position',[0.25+0.12 0.19 0.355 0.73])
legend_text={};
plot(nan,nan,'-k')
legend_text{end+1}='{\itP}_{S}=36kW';
hold on
plot(nan,nan,'-.k')
legend_text{end+1}='{\itP}_{S}=55kW';
set(gca,'ColorOrderIndex',1);
plot(p,y36kW(experiment.PER),'-o')
plot(p,y55kW(experiment.PER),'-.o')
text(0.22,2,'experiment')
set(gca,'ColorOrderIndex',1);
plot(p,y36kW(fem.PER),'-s')
plot(p,y55kW(fem.PER),'-.s')
text(0.22,7.5,'FEM')
set(gca,'ColorOrderIndex',1);
plot(p,y36kW(source_t.PER),'-d')
plot(p,y55kW(source_t.PER),'-.d')
text(0.22,14,'transformer')
% plot(p,d2(fem.PER,1),'--s')
% legend_text{end+1}='{\itP}_{S}=18kW';
% plot(p,y36kW(source_j.PER),'-x')
% plot(p,y55kW(source_j.PER),'-.x')
ylabel('PER [\Omega]')
xlabel('{\itp}_{fill} [Pa]')
xticks(p)
grid on
L1=legend(legend_text);
set(L1,'location','west');
set(L1,'box','off')
min_y=0;
max_y=17;
axis([0.2,0.5,min_y,max_y]) 

% 手动调整图例位置
fig_name='Fig14';
ft.print_fed('vec_emf', [path fig_name])
ft.print_fed('vec_eps', [path fig_name])


% 
% % PER
% figure
% % experiment
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% plot(p,y36kW(experiment.PER),'-o')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,y55kW(experiment.PER),'-.o')
% legend_text{end+1}='{\itP}_{S}=55kW';
% ylabel('PER [\Omega]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% L1=legend(legend_text);
% set(L1,'location','best');
% set(L1,'box','off')
% max_y=0.9;
% axis([0.2,0.5,0.5,max_y]) 
% text(0.25,0.95*(max_y-0.5)+0.5,'(a)')
% % model
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% plot(nan,nan,'-')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(nan,nan,'-.')
% legend_text{end+1}='{\itP}_{S}=55kW';
% set(gca,'ColorOrderIndex',1);
% plot(p,y36kW(fem.PER),'-s')
% plot(p,y55kW(fem.PER),'-.s')
% text(0.22,7.5,'FEM')
% set(gca,'ColorOrderIndex',1);
% plot(p,y36kW(source_t.PER),'-d')
% plot(p,y55kW(source_t.PER),'-.d')
% text(0.22,14,'transformer')
% % plot(p,d2(fem.PER,1),'--s')
% % legend_text{end+1}='{\itP}_{S}=18kW';
% ylabel('PER [\Omega]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% L1=legend(legend_text);
% set(L1,'location','best');
% set(L1,'box','off')
% min_y=5;
% max_y=17;
% axis([0.2,0.5,min_y,max_y]) 
% text(0.3,0.95*(max_y-min_y)+min_y,'(b)')
% 

% % PER
% figure
% % all
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% plot(nan,nan,'-k')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(nan,nan,'-.k')
% legend_text{end+1}='{\itP}_{S}=55kW';
% set(gca,'ColorOrderIndex',1);
% plot(p,y36kW(experiment.PER),'-o')
% plot(p,y55kW(experiment.PER),'-.o')
% text(0.22,2,'experiment')
% set(gca,'ColorOrderIndex',1);
% plot(p,y36kW(fem.PER),'-s')
% plot(p,y55kW(fem.PER),'-.s')
% text(0.22,7.5,'FEM')
% set(gca,'ColorOrderIndex',1);
% plot(p,y36kW(source_t.PER),'-d')
% plot(p,y55kW(source_t.PER),'-.d')
% text(0.22,14,'transformer')
% % plot(p,d2(fem.PER,1),'--s')
% % legend_text{end+1}='{\itP}_{S}=18kW';
% % plot(p,y36kW(source_j.PER),'-x')
% % plot(p,y55kW(source_j.PER),'-.x')
% ylabel('PER [\Omega]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% L1=legend(legend_text);
% set(L1,'location','west');
% set(L1,'box','off')
% min_y=0;
% max_y=17;
% axis([0.2,0.5,min_y,max_y]) 
% text(0.25,0.95*(max_y-0.5)+0.5,'(a)')
% % FEM and experiment
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% plot(nan,nan,'-')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(nan,nan,'-.')
% legend_text{end+1}='{\itP}_{S}=55kW';
% set(gca,'ColorOrderIndex',1);
% plot(p,y36kW(fem.PER),'-s')
% plot(p,y55kW(fem.PER),'-.s')
% text(0.22,7.5,'FEM')
% set(gca,'ColorOrderIndex',1);
% plot(p,y36kW(experiment.PER),'-o')
% plot(p,y55kW(experiment.PER),'-.o')
% text(0.22,14,'transformer')
% % plot(p,d2(fem.PER,1),'--s')
% % legend_text{end+1}='{\itP}_{S}=18kW';
% % ylabel('PER [\Omega]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% % L1=legend(legend_text);
% % set(L1,'location','best');
% % set(L1,'box','off')
% min_y=0;
% max_y=1.5;
% axis([0.2,0.5,min_y,max_y]) 
% text(0.3,0.95*(max_y-min_y)+min_y,'(b)')

% Rs, Ls
figure
% Rs
subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
legend_text={};
plot(nan,nan,'-k')
legend_text{end+1}='{\itP}_{S}=36kW';
hold on
plot(nan,nan,'-.k')
legend_text{end+1}='{\itP}_{S}=55kW';
set(gca,'ColorOrderIndex',1);
plot(p,y36kW(experiment.PER+input.external.Rmetal),'-o')
plot(p,y55kW(experiment.PER+input.external.Rmetal),'-.o')
text(0.22,2,'experiment')
set(gca,'ColorOrderIndex',1);
plot(p,y36kW(fem.PER+fem.Rmetal),'-s')
plot(p,y55kW(fem.PER+fem.Rmetal),'-.s')
text(0.22,6.5,'FEM')
ylabel('{\itR}_s [\Omega]')
xlabel('{\itp}_{fill} [Pa]')
xticks(p)
grid on
L1=legend(legend_text);
set(L1,'location','best');
set(L1,'box','off')
min_y=0;
max_y=8.5;
axis([0.2,0.5,min_y,max_y]) 
text(0.25,0.95*(max_y-min_y)+min_y,'(a)')
% Ls
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
legend_text={};
plot(nan,nan,'-k')
legend_text{end+1}='{\itP}_{S}=36kW';
hold on
plot(nan,nan,'-.k')
legend_text{end+1}='{\itP}_{S}=55kW';
set(gca,'ColorOrderIndex',1);
plot(p,y36kW(fem.Ls),'-s')
plot(p,y55kW(fem.Ls),'-.s')
% text(0.4,7.5,'FEM')
set(gca,'ColorOrderIndex',1);
plot(p,1e6*y36kW(source_t.Lsys),'-d')
plot(p,1e6*y55kW(source_t.Lsys),'-.d')
text(0.3,6,'transformer')
% plot(p,1e6*y55kW(source_j.Lsys),'--d')
% ylabel('{\itL}_s [\muH]')
xlabel('{\itp}_{fill} [Pa]')
xticks(p)
grid on
L1=legend(legend_text);
set(L1,'location','best');
set(L1,'box','off')
min_y=6.5;
max_y=7.6;
axis([0.2,0.5,min_y,max_y]) 
text(0.25,0.95*(max_y-min_y)+min_y,'(b)')

% 手动调整图例位置
% visio处理，导出emf和tif
fig_name='Fig16_matlab';
ft.print_fed('vec_emf', [path fig_name])
ft.print_fed('vec_eps', [path fig_name])

% % Rs
% figure
% % experiment
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% plot(p,y36kW(experiment.PER+input.external.Rmetal),'-o')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,y55kW(experiment.PER+input.external.Rmetal),'-.o')
% legend_text{end+1}='{\itP}_{S}=55kW';
% ylabel('{\itR}_S [\Omega]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% L1=legend(legend_text);
% set(L1,'location','best');
% set(L1,'box','off')
% min_value=0.5+0.6;
% max_value=0.9+0.6;
% axis([0.2,0.5,min_value,max_value]) 
% text(0.25,0.95*(max_value-min_value)+min_value,'(a)')
% % model
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% plot(p,y36kW(fem.PER+fem.Rmetal),'-s')
% legend_text{end+1}='{\itP}_{S}=36kW';
% hold on
% plot(p,y55kW(fem.PER+fem.Rmetal),'-.s')
% legend_text{end+1}='{\itP}_{S}=55kW';
% ylabel('{\itR}_S [\Omega]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% yticks([6.9,7.2,7.6,8])
% grid on
% L1=legend(legend_text);
% set(L1,'location','south');
% set(L1,'box','off')
% min_value=6.7;
% max_value=8;
% axis([0.2,0.5,min_value,max_value]) 
% text(0.25,0.95*(max_value-min_value)+min_value,'(b)')



%% magnetized
% w_ce
Bsurface=source_a.emf.Bzm_r(plasma.r);
wce_max=constants.e*Bsurface/constants.me;
plasma.vte=sqrt(plasma.Te*constants.e/constants.me);
r_L_min=constants.me*plasma.vte./(constants.e*Bsurface);
wce_skindepth=zeros(3,1,3);
for i_p=1:3
    for i_Pin=1:3
    Bskindepth=source_a.emf.Bzm_r(plasma.r-plasma.skin_depth(i_p,1,i_Pin));
    wce_skindepth(i_p,1,i_Pin)=constants.e*Bskindepth(i_p,1,i_Pin)/constants.me;
    end
end

figure
% skindepth and wavelength
% subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
legend_text={};
plot(p,1e3*y36kW(plasma.skin_depth),'-o')
legend_text{end+1}='{\itP}_{S}=36kW';
hold on
plot(p,1e3*y55kW(plasma.skin_depth),'-.o')
legend_text{end+1}='{\itP}_{S}=55kW';
line([0.2,0.5],1e3*[plasma.r, plasma.r],'linestyle','--','color','k');
legend_text{end+1}='{\itr}_{plasma}';
plot(p,1e3*y36kW(r_L_min),'-s')
plot(p,1e3*y55kW(r_L_min),'-.s')
ylabel('{length} [mm]')
xlabel('{\itp}_{fill} [Pa]')
xticks(p)
grid on
L1=legend(legend_text);
set(L1,'location','south');
set(L1,'box','off')
max_value=240;
axis([0.2,0.5,0,max_value]) 

% % frequency
% figure
% % 36kW
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% semilogy(p,y36kW(plasma.nu_eff),'-o')
% legend_text{end+1}='{\it\nu}_{eff}';
% hold on
% semilogy(p,y36kW(wce_skindepth),'-s')
% legend_text{end+1}='{\it\omega}_{ce}';
% ylabel('{frequency [s^{-1}]}')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% line([0.2,0.5],[w1MHz, w1MHz],'linestyle','--','color','k');
% legend_text{end+1}='{\it\omega}_{RF}';
% L1=legend(legend_text);
% set(L1,'location','east');
% set(L1,'box','off')
% max_y=4e9;
% axis([0.2,0.5,5e6,max_y]) 
% text(0.25,2.3e9,'(a) {\itP}_{S}=36kW')
% % 55kW
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% semilogy(p,y55kW(plasma.nu_eff),'-.o')
% legend_text{end+1}='{\it\nu}_{eff}';
% hold on
% semilogy(p,y55kW(wce_skindepth),'-.s')
% legend_text{end+1}='{\it\omega}_{ce}';
% ylabel('{frequency [s^{-1}]}')
% xlabel('{\itp}_{fill} [Pa]')
% xticks(p)
% grid on
% line([0.2,0.5],[w1MHz, w1MHz],'linestyle','--','color','k');
% legend_text{end+1}='{\it\omega}_{RF}';
% L1=legend(legend_text);
% set(L1,'location','east');
% set(L1,'box','off')
% max_y=4e9;
% axis([0.2,0.5,5e6,max_y]) 
% text(0.25,2.3e9,'(b) {\itP}_{S}=55kW')

% frequency
figure
% nu_m,nu_st
subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
legend_text={};
semilogy(nan,nan,'-k')
legend_text{end+1}='{\itP}_{S}=36kW';
hold on
semilogy(nan,nan,'-.k')
legend_text{end+1}='{\itP}_{S}=55kW';
set(gca,'ColorOrderIndex',1);
semilogy(p,y36kW(plasma.nu_m),'-o')
semilogy(p,y36kW(plasma.nu_st),'-s')
set(gca,'ColorOrderIndex',1);
semilogy(p,y55kW(plasma.nu_m),'-.o')
semilogy(p,y55kW(plasma.nu_st),'-.s')
ylabel('{frequency [s^{-1}]}')
xlabel('{\itp}_{fill} [Pa]')
xticks(p)
grid on
line([0.2,0.5],[w1MHz, w1MHz],'linestyle','--','color','k');
text(0.4,0.4e7,'{\it\nu}_{m}')
text(0.4,1.4e7,'{\it\nu}_{st}')
text(0.23,0.55e7,'{\it\omega}');
L1=legend(legend_text);
set(L1,'location','east');
set(L1,'box','off')
max_y=2e7;
axis([0.2,0.5,0,max_y]) 
text(0.25,1.8e7,'(a)')
% w_ce, nu_eff
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
legend_text={};
semilogy(nan,nan,'-k')
legend_text{end+1}='{\itP}_{S}=36kW';
hold on
semilogy(nan,nan,'-.k')
legend_text{end+1}='{\itP}_{S}=55kW';
set(gca,'ColorOrderIndex',3);
semilogy(p,y36kW(plasma.nu_eff),'-d')
semilogy(p,y36kW(wce_skindepth),'->')
set(gca,'ColorOrderIndex',3);
semilogy(p,y55kW(plasma.nu_eff),'-.d')
semilogy(p,y55kW(wce_skindepth),'-.>')
ylabel('{frequency [s^{-1}]}')
xlabel('{\itp}_{fill} [Pa]')
xticks(p)
grid on
line([0.2,0.5],[w1MHz, w1MHz],'linestyle','--','color','k');
text(0.4,4e7,'{\it\nu}_{eff}')
text(0.4,2.3e9,'{\it\omega}_{ce}')
text(0.23,9e6,'{\it\omega}');
L1=legend(legend_text);
set(L1,'location','east');
set(L1,'box','off')
max_y=4e9;
axis([0.2,0.5,5e6,max_y]) 
text(0.25,2.8e9,'(b)')

% 手动调整图例位置
% visio处理，导出emf和tif
fig_name='Fig15';
ft.print_fed('vec_emf', [path fig_name])
ft.print_fed('vec_eps', [path fig_name])


% 0.3Pa, 55kW
idx=(plasma.p==0.3 & plasma.Pin==55e3);
% B(r)
size_mat=size(source_a.PER);
r1=0;
r2=input.geometry.r_plasma_eff;
r3=input.geometry.r_coil;
r_a=[r1:(r2-r1)/100:r2 r2+(r3-r2)/10:(r3-r2)/10:r3];
% 前101个在r_plasma内，后10个在>r_plasma
len_r=length(r_a);

B_a=zeros(1,len_r);
B_ma=zeros(1,len_r);
for i_r=1:len_r
    temp=source_a.emf.Bzm_r(r_a(i_r));
    B_a(i_r)=temp(idx);
    temp=source_ma.emf.Bzm_r(r_a(i_r));
    B_ma(i_r)=temp(idx);
end

% B_rcoil_a=B_a(end);
% B_rplasma_a=B_a(101);
% B_mean_plasma=input_m.plasma.Bz_mean(idx);

% FEM model中Bz，变换到同一Im下
Im=sqrt(2)*input.external.Icoil_rms(idx);

radialHz= csvread('cHz_0.3Pa.csv',1,0);
fem.radial_r=radialHz(:,1);
radialHz(radialHz(:,2)==0,2)=nan;
radialHz(radialHz(:,3)==0,3)=nan;
fem.radial_Bz=Im/25.71*constants.mu0*radialHz(:,2);
fem.radial_Bz_no_plasma=Im/21.65*constants.mu0*radialHz(:,3);
radialHz= csvread('cHz_0.3Pa_no_FS.csv',1,0);
% fem.radial_r2=radialHz(:,1);
% find(~(fem.radial_r==fem.radial_r2))
radialHz(radialHz(:,2)==0,2)=nan;
radialHz(radialHz(:,3)==0,3)=nan;
fem.radial_Bz_no_FS=Im/23.21*constants.mu0*radialHz(:,2);
fem.radial_Bz_no_FS_no_plasma=Im/16.92*constants.mu0*radialHz(:,3);



% B(r)
figure
legend_text={};
max_y=23;
rectangle('Position',[116,16,3,max_y],'Curvature', [0 0], 'FaceColor',[196 203 207]/255)
hold on
plot(fem.radial_r,1e3*fem.radial_Bz,'-')
legend_text{end+1}='FEM';
plot(fem.radial_r,1e3*fem.radial_Bz_no_FS,'--')
legend_text{end+1}='FEM, no FS';
plot(fem.radial_r,1e3*fem.radial_Bz_no_plasma,'--')
legend_text{end+1}='FEM, no plasma';
plot(fem.radial_r,1e3*fem.radial_Bz_no_FS_no_plasma,'--')
legend_text{end+1}='FEM, no FS no plasma';
plot(1e3*r_a,1e3*B_a,'-.k')
legend_text{end+1}='analytical, no FS';
ylabel('{\itB}_z [mT]')
xlabel('{\itr} [mm]')
grid on
L1=legend(legend_text);
set(L1,'location','south');
set(L1,'box','off')
min_y=0;
% max_x=119;
max_x=1e3*input.geometry.r_coil;
axis([0,max_x,min_y,max_y]) 
text(max_x/6,0.95*max_y,'(a)')
text(max_x/6,0.95*(max_y-min_y)+min_y,'{\itp}_{fill}=0.3Pa, {\itP}_{S}=55kW')


% 
Bfun_e=@(r) interp1(fem.radial_r,1e3*fem.radial_Bz,r);
B_at_r_plasma=Bfun_e(r_plasma);
B_by_e=B_at_r_plasma/exp(1);
B_at_delta=Bfun_e(r_plasma-1e3*plasma.skin_depth(2,1,3));

% 指数衰减？
figure
semilogy(fem.radial_r,1e3*fem.radial_Bz,'-')
hold on 



