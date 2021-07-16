% CHARLIE fig for paper

% 图形美化-位置、大小合适等要根据论文来，而且别浪费时间

% 整体上，突出重点；单张图上，不要混杂
% 重点：本文FEM model优于现有阻抗和效率估算方法
% 优于文献常用的解析电磁和变压器模型

% 次要-讨论
% 文献常用的良导体-FEM model，与有损介质FEM model
% 不均匀分布的影响

% 输入：ne，Te

% 输出
% 论文只对比关键参数，不对比过程值
% 效率：PER
% 阻抗：Rs，Ls
path='d:\School\DoctorProgram\eP-项目笔记文件夹\eP-190821-01激励器FEM模型\paper_v2-210427期刊论文\绘图\';
% ft.print_fed('vec_emf',[path 'Fig4'])

% - 1MHz -. 4MHz 线型 频率  
% or / sb / dc 颜色与标记 物理量
% --k 辅助线

clear
solution_name='for_paper210415';
% pwd == the root of SAD code
addpath(genpath(['./others/' solution_name '/']))
addpath(genpath('./packages'))

close all

%% input
% get_CHARLIE_resultsmat();
load('CHARLIE_results.mat')
load('paper_CHARLIE_raza_t210517.mat')
plasma=input.plasma;
source_t=source;
load('paper_CHARLIE_raza_a210517.mat')
source_a=source;
plasma_ma=input_m.plasma;
load('paper_CHARLIE_raza_tj210517.mat')
plasma_j=input.plasma;
source_j=source;

constants=get_constants();
ft=my.figtool;
ft.init_fed();
p=[0.3, 0.5, 1, 3, 5, 10]';
w1MHz=2*pi*1e6;
w4MHz=2*pi*4e6;
co=my.get_color_order('');

%% input
% 考虑ne不均匀分布，估算平均值/边缘值
norm_ne_r10=nonuniform_dist.get_ne_r_CHARLIE(10); % origin data from experiments
for i=1:5
    dist_rp(i)=nonuniform_dist.get_nonuniform_dist_CHARLIE(['rp' num2str(i)]);
    ratio_origin2goal(i).ne_r=dist_rp(i).ne_r/norm_ne_r10;
end
data_ref=nonuniform_dist.get_ref_CHARLIE();
r_line=(0:45.5/100:45.5)';
ne_r_fit=nonuniform_dist.get_ne_r_CHARLIE(r_line);
Te_r_fit=nonuniform_dist.get_Te_r_CHARLIE(r_line);

z_line=(-200:10:200)';
ne_z_fit=nonuniform_dist.get_ne_z_CHARLIE(z_line);
Te_z_fit=nonuniform_dist.get_Te_z_CHARLIE(z_line);
mean_ne_z=nonuniform_dist.get_ne_z_CHARLIE([-200,200]);
mean_coil_ne_z=nonuniform_dist.get_ne_z_CHARLIE([-50,50]);

% ne 分布
figure
% 径向
subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
plot(r_line,ne_r_fit) 
hold on
temp_color=co(2,:);
line([0,45.5],[dist_rp(1).ne_r(1) dist_rp(1).ne_r(1)],'linestyle','--','color',temp_color);
axis([0,45.5,0,1]) 
xlabel('{\itr} [mm]')
xticks([0,20,45.5])
ylabel('normalized {\itn}_e')
grid on
line([10,10],[0.8, 1],'linestyle','--','color','k','LineWidth',1);
text(2,0.7,'diagnostic\newlineport')
L1=legend('PIC/MCC','averaged');
set(L1,'location','southwest');
set(L1,'box','off')
text(45.5/3,0.95,'(a)')
% 轴向
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
rectangle('Position',[0,0,50,1],'LineStyle', 'none', 'FaceColor',[196 203 207]/255)
hold on 
idx=z_line>=0;
plot(z_line(idx),ne_z_fit(idx))
line([0,200],[mean_ne_z, mean_ne_z],'linestyle','--','color',temp_color);
axis([0,200,0,1]) 
xlabel('{\itz} [mm]')
xticks(0:100:200)
ylabel('normalized {\itn}_e')
grid on
text(1,0.6,'coil')
L1=legend('experiment','averaged');
set(L1,'location','north');
set(L1,'box','off')
text(200/3,0.95,'(b)')
% % 轴向
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% rectangle('Position',[-50,0,100,1],'LineStyle', 'none', 'FaceColor',[196 203 207]/255)
% hold on 
% plot(z_line,ne_z_fit)
% line([-200,200],[mean_ne_z, mean_ne_z],'color',[0.8500    0.3250    0.0980]);
% axis([-200,200,0,1]) 
% xlabel('{\itz} [mm]')
% xticks(-200:200:200)
% ylabel('normalized {\itn}_e')
% grid on
% text(-30,0.6,'coil')
% L1=legend('experiment','mean');
% set(L1,'location','south');
% set(L1,'box','off')
% text(-200+400/6,0.95,'(b)')

% ne, Te
figure
% ne
subplot(1,2,1);
% set(gca, 'Position',[0.1 0.19 0.355 0.73])
legend_text={};
temp=co(1,:);
semilogx(p,plasma.ne(:,1),'-o','Color',temp)
legend_text{end+1}='{\itf}=1MHz';
hold on
semilogx(p,plasma.ne(:,2),'-.o','Color',temp)
legend_text{end+1}='{\itf}=4MHz';
ylabel('{\itn}_e [m^{-3}]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southeast');
set(L1,'box','off')
max_y=2e17;
axis([0.3,10,0,max_y]) 
text(0.4,0.95*max_y,'(a)')
% Te
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
legend_text={};
temp=co(2,:);
semilogx(p,plasma.Te(:,1),'-s','Color',temp)
legend_text{end+1}='{\itf}=1MHz';
hold on
semilogx(p,plasma.Te(:,2),'-.s','Color',temp)
legend_text{end+1}='{\itf}=4MHz';
ylabel('{\itT}_e [eV]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southeast');
set(L1,'box','off')
max_y=6;
axis([0.3,10,0,max_y]) 
text(0.4,0.95*max_y,'(b)')

%% 等离子体模型输出
% % sigma, eps_r
% figure
% % sigma
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% semilogx(p,plasma.sigma(:,1),'-o')
% legend_text{end+1}='{\itf}=1MHz';
% hold on
% semilogx(p,plasma.sigma(:,2),'-.o')
% legend_text{end+1}='{\itf}=4MHz';
% ylabel('\sigma [S/m]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% L1=legend(legend_text);
% set(L1,'location','south');
% set(L1,'box','off')
% min_y=0;
% max_y=1e2;
% axis([0.3,10,min_y,max_y]) 
% text(2,93,'(a)')
% % 4MHz
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% semilogx(p,-w1MHz*plasma.eps_prime(:,1),'-s')
% legend_text{end+1}='{\itf}=1MHz';
% hold on
% semilogx(p,-w4MHz*plasma.eps_prime(:,2),'-.s')
% legend_text{end+1}='{\itf}=4MHz';
% ylabel('-{\it\omega\epsilon} [S/m]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% L1=legend(legend_text);
% set(L1,'location','north');
% set(L1,'box','off')
% min_y=0;
% max_y=1e2;
% axis([0.3,10,min_y,max_y]) 
% text(2,93,'(b)')

r_plasma=1e3*plasma.r;


% % length
% figure
% % skindepth
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% semilogx(p,1e3*plasma.skin_depth(:,1),'-o')
% legend_text{end+1}='{\itf}=1MHz';
% hold on
% semilogx(p,1e3*plasma.skin_depth(:,2),'-.o')
% legend_text{end+1}='{\itf}=4MHz';
% ylabel('\delta [mm]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% line([0.3,10],[r_plasma, r_plasma],'linestyle','--','color','k');
% legend_text{end+1}='{\itr}_{plasma}';
% L1=legend(legend_text);
% set(L1,'location','north');
% set(L1,'box','off')
% min_y=0;
% max_y=1e2;
% axis([0.3,10,min_y,max_y]) 
% text(0.7,0.95*max_y,'(a)')
% % wavelength
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% semilogx(p,1e3*plasma.wavelength(:,1),'-o')
% legend_text{end+1}='{\itf}=1MHz';
% hold on
% semilogx(p,1e3*plasma.wavelength(:,2),'-.o')
% legend_text{end+1}='{\itf}=4MHz';
% ylabel('\lambda [mm]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% line([0.3,10],[r_plasma, r_plasma],'linestyle','--','color','k');
% legend_text{end+1}='{\itr}_{plasma}';
% L1=legend(legend_text);
% set(L1,'location','north');
% set(L1,'box','off')
% min_y=0;
% max_y=7e2;
% axis([0.3,10,min_y,max_y]) 
% text(0.7,0.95*max_y,'(b)')

%% main: different electric models
% FEM model，解析电磁，变压器model，实验
% PER
figure
%1MHz
subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
legend_text={};
semilogx(p,experiment.PER(:,1),'-o')
legend_text{end+1}='experiment';
hold on
semilogx(p,fem.dielectric_PER(:,1),'-s')
legend_text{end+1}='FEM';
semilogx(p,source_a.PER(:,1),'-d')
legend_text{end+1}='analytical';
% semilogx(p,source_t.PER(:,1),'->')
% legend_text{end+1}='transformer-base';
% semilogx(p,source_j.PER(:,1),'-x')
% legend_text{end+1}='transformer-Jain';
ylabel('PER [\Omega]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','north');
set(L1,'box','off')
min_y=0;
max_y=4;
axis([0.3,10,min_y,max_y]) 
text(0.4,0.95*max_y,'(a) {\itf}=1MHz')
% 4MHz
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
legend_text={};
semilogx(p,experiment.PER(:,2),'-.o')
legend_text{end+1}='experiment';
hold on
semilogx(p,fem.dielectric_PER(:,2),'-.s')
legend_text{end+1}='FEM';
semilogx(p,source_a.PER(:,2),'-.d')
legend_text{end+1}='analytical';
% semilogx(p,source_t.PER(:,2),'->')
% legend_text{end+1}='transformer-base';
% semilogx(p,source_j.PER(:,2),'-x')
% legend_text{end+1}='transformer-Jain';
ylabel('PER [\Omega]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','north');
set(L1,'box','off')
min_y=0;
max_y=24;
axis([0.3,10,min_y,max_y]) 
text(0.4,0.95*max_y,'(b) {\itf}=4MHz')

% Rs, Ls
figure
% Rs
subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
legend_text={};
semilogx(nan,nan,'-o')
legend_text{end+1}='experiment';
hold on
semilogx(nan,nan,'-s')
legend_text{end+1}='FEM';
set(gca,'ColorOrderIndex',1);
semilogx(p,experiment.PER(:,1)+input.external.Rmetal(:,1),'-o')
semilogx(p,fem.dielectric_PER(:,1)+fem.dielectric_Rmetal(:,1),'-s')
text(3,1.8,'{\itf}=1MHz')
set(gca,'ColorOrderIndex',1);
semilogx(p,experiment.PER(:,2)+input.external.Rmetal(:,2),'-.o')
semilogx(p,fem.dielectric_PER(:,2)+fem.dielectric_Rmetal(:,2),'-.s')
text(1,6.1,'{\itf}=4MHz')
ylabel('{\itR}_s [\Omega]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','north');
set(L1,'box','off')
min_y=0;
max_y=10;
axis([0.3,10,min_y,max_y]) 
text(0.4,0.95*max_y,'(a)')
% Ls
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
legend_text={};
set(gca,'ColorOrderIndex',2);
semilogx(nan,nan,'-s')
legend_text{end+1}='FEM';
hold on
semilogx(nan,nan,'-d')
legend_text{end+1}='analytical';
set(gca,'ColorOrderIndex',2);
semilogx(p,fem.dielectric_Ls(:,1),'-s')
semilogx(p,1e6*source_a.Lsys(:,1),'-d')
text(1,3.3,'{\itf}=1MHz')
set(gca,'ColorOrderIndex',2);
semilogx(p,fem.dielectric_Ls(:,2),'-.s')
semilogx(p,1e6*source_a.Lsys(:,2),'-.d')
text(1,1.8,'{\itf}=4MHz')
ylabel('{\itL}_s [\muH]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','south');
set(L1,'box','off')
min_y=0;
max_y=4;
axis([0.3,10,min_y,max_y]) 
text(0.4,0.95*max_y,'(b)')

% % Rs, Ls
% figure
% % Rs
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% semilogx(nan,nan,'-o')
% legend_text{end+1}='experiment';
% hold on
% semilogx(nan,nan,'-s')
% legend_text{end+1}='FEM';
% set(gca,'ColorOrderIndex',1);
% semilogx(p,input.external.Rmetal(:,1),'-o')
% semilogx(p,fem.dielectric_Rmetal(:,1),'-s')
% % text(1,2.1,'{\itf}=1MHz')
% set(gca,'ColorOrderIndex',1);
% semilogx(p,input.external.Rmetal(:,2),'-.o')
% semilogx(p,fem.dielectric_Rmetal(:,2),'-.s')
% % text(1,6.1,'{\itf}=4MHz')
% ylabel('{\itR}_{metal} [\Omega]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% L1=legend(legend_text);
% set(L1,'location','north');
% set(L1,'box','off')
% min_y=0;
% max_y=10;
% axis([0.3,10,min_y,max_y]) 
% text(0.4,0.95*max_y,'(a)')
% % Ls
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% semilogx(nan,nan,'-s')
% legend_text{end+1}='FEM';
% hold on
% semilogx(nan,nan,'-d')
% legend_text{end+1}='analytical';
% set(gca,'ColorOrderIndex',1);
% semilogx(p,fem.dielectric_Ls(:,1),'-s')
% semilogx(p,1e6*source_a.Lsys(:,1),'-d')
% text(1,3.3,'{\itf}=1MHz')
% set(gca,'ColorOrderIndex',1);
% semilogx(p,fem.dielectric_Ls(:,2),'-.s')
% semilogx(p,1e6*source_a.Lsys(:,2),'-.d')
% text(1,1.8,'{\itf}=4MHz')
% ylabel('{\itL}_s [\muH]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% L1=legend(legend_text);
% set(L1,'location','south');
% set(L1,'box','off')
% min_y=0;
% max_y=4;
% axis([0.3,10,min_y,max_y]) 
% text(0.4,0.95*max_y,'(b)')

% ft.print_fed('vec_emf','test')

%% 加热机制
% nu_c
figure
% 1MHz
subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.75])
legend_text={};
temp=co(1,:);
loglog(nan,nan,'-','Color',temp)
legend_text{end+1}='{\itf}=1MHz';
hold on
semilogx(nan,nan,'-.','Color',temp)
legend_text{end+1}='{\itf}=1MHz';
temp=co(1,:);
loglog(p,plasma.nu_m(:,1),'-o','Color',temp)
loglog(p,plasma.nu_m(:,2),'-.o','Color',temp)
text(9,1.5e8,'{\it\nu}_{m}')
temp=co(2,:);
loglog(p,plasma.nu_st(:,1),'-s','Color',temp)
loglog(p,plasma.nu_st(:,2),'-.s','Color',temp)
text(1,4e6,'{\it\nu}_{st}')
temp=co(3,:);
line([0.3,10],[w1MHz, w1MHz],'linestyle','-','color',temp);
line([0.3,10],[w4MHz, w4MHz],'linestyle','-.','color',temp);
text(3,4e6,'{\it\omega}')
ylabel('frequency [s^{-1}]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','northwest');
set(L1,'box','off')
min_y=3e6;
max_y=2e8;
axis([0.3,10,min_y,max_y]) 
text(1.2,1.6e8,'(a)')
% error-ratio
subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.75])
legend_text={};
temp=co(4,:);
semilogx(plasma.nu_m(:,1)./plasma.nu_st(:,1),...
    fem.dielectric_PER(:,1)./experiment.PER(:,1),'->','Color',temp)
legend_text{end+1}='{\itf}=1MHz';
hold on
semilogx(plasma.nu_m(:,2)./plasma.nu_st(:,2),...
    fem.dielectric_PER(:,2)./experiment.PER(:,2),'-.>','Color',temp)
legend_text{end+1}='{\itf}=4MHz';
ylabel('PER_{FEM}/PER_{experiment}')
xlabel('\nu_{m}/\nu_{st}')
xticks([0.3,1,10])
grid on
min_y=0;
max_y=10;
axis([3e-1,2e1,min_y,max_y]) 
line([1,1],[min_y,max_y],'linestyle','--','color','k');
L1=legend(legend_text);
set(L1,'location','northwest');
set(L1,'box','off')
text(1.4,9.5,'(b)')

% % nu_c
% % loglog 与其他图semilogx不一致，且难以突出不同物理量数量对比
% % 但有利于体现nu_m与p成指数关系
% figure
% % 1MHz
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.81])
% legend_text={};
% loglog(p,plasma.nu_m(:,1),'-o')
% legend_text{end+1}='{\it\nu}_{m}';
% hold on
% loglog(p,plasma.nu_st(:,1),'-s')
% legend_text{end+1}='{\it\nu}_{st}';
% ylabel('frequency [s^{-1}]')
% % yticks([3e6,1e7,1e8,2e8])
% % yticklabels({'3\times10^6','10^7','10^8','2\times10^8'})
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% line([0.3,10],[w1MHz, w1MHz],'linestyle','--','color','k');
% legend_text{end+1}='\omega';
% L1=legend(legend_text);
% set(L1,'location','northwest');
% set(L1,'box','off')
% min_y=3e6;
% max_y=2e8;
% axis([0.3,10,min_y,max_y]) 
% text(0.7,1.6e8,'(a) {\itf}=1MHz')
% % 4MHz
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.81])
% loglog(p,plasma.nu_m(:,2),'-.o')
% hold on
% loglog(p,plasma.nu_st(:,2),'-.s')
% ylabel('frequency [s^{-1}]')
% % yticks([3e6,1e7,1e8,2e8])
% % yticklabels({'3\times10^6','10^7','10^8','2\times10^8'})
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% line([0.3,10],[w4MHz, w4MHz],'linestyle','--','color','k');
% L1=legend(legend_text);
% set(L1,'location','northwest');
% set(L1,'box','off')
% min_y=3e6;
% max_y=2e8;
% axis([0.3,10,min_y,max_y]) 
% text(0.7,1.6e8,'(b) {\itf}=4MHz')
% nu_m中nu_enp占绝大多数，没必要画图

% % frequency ratio
% % loglog 与其他图semilogx不一致，且难以突出不同物理量数量对比
% % 但有利于体现nu_m与p成指数关系
% figure
% % nu_m/nu_st
% subplot(1,2,1);
% % set(gca, 'Position',[0.12 0.19 0.355 0.81])
% legend_text={};
% temp=co(1,:);
% loglog(p,plasma.nu_m(:,1)./plasma.nu_st(:,1),'-o','Color',temp)
% legend_text{end+1}='{\itf}=1MHz';
% hold on
% loglog(p,plasma.nu_m(:,2)./plasma.nu_st(:,2),'-.o','Color',temp)
% legend_text{end+1}='{\itf}=4MHz';
% ylabel('\nu_{m}/\nu_{st}')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% line([0.3,10],[1, 1],'linestyle','--','color','k');
% L1=legend(legend_text);
% set(L1,'location','northwest');
% set(L1,'box','off')
% min_y=1e-1;
% max_y=1e2;
% axis([0.3,10,min_y,max_y]) 
% text(0.4,1.8e-1,'(a)')
% % nu_eff/omega
% subplot(1,2,2);
% % set(gca, 'Position',[0.61 0.19 0.355 0.81])
% temp=co(2,:);
% loglog(p,plasma.nu_eff(:,1)./plasma.w_RF(:,1),'-s','Color',temp)
% hold on
% loglog(p,plasma.nu_eff(:,2)./plasma.w_RF(:,2),'-.s','Color',temp)
% ylabel('\nu_{eff}/\omega')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% line([0.3,10],[1, 1],'linestyle','--','color','k');
% L1=legend(legend_text);
% set(L1,'location','northwest');
% set(L1,'box','off')
% min_y=1e-1;
% max_y=1e2;
% axis([0.3,10,min_y,max_y]) 
% text(0.4,1.8e-1,'(b)')

%% 磁化
% 磁场分布
% 0.3Pa, 1MHz
idx1=(plasma.p==0.3 & plasma.f==1e6);
idx2=(plasma.p==0.3 & plasma.f==4e6);
% B(r)
size_mat=size(source_a.PER);
r1=0;
r2=input.geometry.r_plasma_eff;
r3=input.geometry.r_coil;
r_a=[r1:(r2-r1)/100:r2 r2+(r3-r2)/10:(r3-r2)/10:r3];

% 前101个在r_plasma内，后10个在>r_plasma
len_r=length(r_a);
B_a1=zeros(1,len_r);
B_a2=B_a1;
for i_r=1:len_r
    temp=source_a.emf.Bzm_r(r_a(i_r));
    B_a1(i_r)=temp(idx1);
    B_a2(i_r)=temp(idx2);
end

figure
legend_text={};
plot(1e3*r_a,1e3*B_a1,'-.')
legend_text{end+1}='analytical, 1MHz0.3Pa';
hold on
plot(1e3*r_a,1e3*B_a2,'-.')
legend_text{end+1}='analytical, 4MHz0.3Pa';
ylabel('{\itB}_z [mT]')
xlabel('{\itr} [mm]')
grid on
L1=legend(legend_text);
set(L1,'location','south');
set(L1,'box','off')

B_max=source_a.emf.Bzm_r(plasma.r);
B_min=source_a.emf.Bzm_r(0); % 与流体模型结果差太多
plasma.vte=sqrt(plasma.Te*constants.e/constants.me);
wce_max=constants.e*B_max/constants.me;
r_L_min=constants.me*plasma.vte./(constants.e*B_max);

% % w_ce at the skindepth
% Bsurface=source_a.emf.Bzm_r(plasma.r);
% wce_max=constants.e*Bsurface/constants.me;
% wce_skindepth=zeros(6,2);
% for i_p=1:6
%     for i_f=1:2
%     Bskindepth=source_a.emf.Bzm_r(plasma.r-plasma.skin_depth(i_p,i_f));
%     wce_skindepth(i_p,i_f)=constants.e*Bskindepth(i_p,i_f)/constants.me;
%     end
% end

figure
% skindepth,r_L
subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
legend_text={};
temp=co(1,:);
semilogx(nan,nan,'-','Color',temp)
legend_text{end+1}='{\itf}=1MHz';
hold on
semilogx(nan,nan,'-.','Color',temp)
legend_text{end+1}='{\itf}=4MHz';
temp=co(1,:);
semilogx(p,1e3*plasma.skin_depth(:,1),'-o','Color',temp)
semilogx(p,1e3*plasma.skin_depth(:,2),'-.o','Color',temp)
text(0.4,30,'\delta')
temp=co(2,:);
semilogx(p,1e3*r_L_min(:,1),'-s','Color',temp)
semilogx(p,1e3*r_L_min(:,2),'-.s','Color',temp)
text(0.6,30,'r_{Larmor}')
ylabel('length [mm]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
line([0.3,10],[r_plasma, r_plasma],'linestyle','--','color','k');
text(3,40,'r_{plasma}')
L1=legend(legend_text);
set(L1,'location','north');
set(L1,'box','off')
min_y=1e0;
max_y=1e2;
axis([0.3,10,min_y,max_y]) 
text(0.4,0.95*max_y,'(a)')
% w_ce, nu_eff
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
legend_text={};
temp=co(3,:);
loglog(p,wce_max(:,1)./plasma.nu_eff(:,1),'-s','Color',temp)
legend_text{end+1}='{\itf}=1MHz';
hold on
loglog(p,wce_max(:,2)./plasma.nu_eff(:,2),'-.s','Color',temp)
legend_text{end+1}='{\itf}=4MHz';
ylabel('\omega_{ce}/\nu_{eff}')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
line([0.3,10],[1, 1],'linestyle','--','color','k');
L1=legend(legend_text);
set(L1,'location','north');
set(L1,'box','off')
min_y=5e-1;
max_y=1e2;
axis([0.3,10,min_y,max_y]) 
text(0.4,8e1,'(b)')

% % length
% figure
% % skindepth
% subplot(1,2,1);
% set(gca, 'Position',[0.12 0.19 0.355 0.73])
% legend_text={};
% semilogx(p,1e3*plasma.skin_depth(:,1),'-o')
% legend_text{end+1}='{\itf}=1MHz';
% hold on
% semilogx(p,1e3*plasma.skin_depth(:,2),'-.o')
% legend_text{end+1}='{\itf}=4MHz';
% ylabel('\delta [mm]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% line([0.3,10],[r_plasma, r_plasma],'linestyle','--','color','k');
% legend_text{end+1}='{\itr}_{plasma}';
% L1=legend(legend_text);
% set(L1,'location','north');
% set(L1,'box','off')
% min_y=0;
% max_y=1e2;
% axis([0.3,10,min_y,max_y]) 
% text(0.7,0.95*max_y,'(a)')
% % wavelength
% subplot(1,2,2);
% set(gca, 'Position',[0.61 0.19 0.355 0.73])
% legend_text={};
% semilogx(p,1e3*plasma.wavelength(:,1),'-o')
% legend_text{end+1}='{\itf}=1MHz';
% hold on
% semilogx(p,1e3*plasma.wavelength(:,2),'-.o')
% legend_text{end+1}='{\itf}=4MHz';
% ylabel('\lambda [mm]')
% xlabel('{\itp}_{fill} [Pa]')
% xticks([0.3,1,3,10])
% grid on
% line([0.3,10],[r_plasma, r_plasma],'linestyle','--','color','k');
% legend_text{end+1}='{\itr}_{plasma}';
% L1=legend(legend_text);
% set(L1,'location','north');
% set(L1,'box','off')
% min_y=0;
% max_y=7e2;
% axis([0.3,10,min_y,max_y]) 
% text(0.7,0.95*max_y,'(b)')


% % 加热区域分析
% % local regime
% delta_l=constants.c./plasma.wpe;
% % anomalous
% % ka=2.8;
% ka=sqrt(pi);
% 
% % nonlinear regime
% % 2009Froese - Nonlinear skin effect in a collisionless plasma
% knl=0.4; % assumed
% B0=Bmax;
% idx1=(plasma.p==0.3 & plasma.f==1e6);
% % delta_nl=
% 
% % 非线性和反常的分界
% temp_right=constants.eps0*constants.c^2*constants.me^5*ka^5/(constants.e^5*knl^6);
% temp_left=constants.e*plasma.Te.*plasma.ne.*B0.^3./plasma.w_RF.^5;
% temp_left-temp_right
% % 20210530 失败 


%% discussion: 不均匀ne的影响
% 先简单展示一下，说明有必要考虑，以后再定量对比
% FEM model 分层
% 以1MHz, 10Pa为例
idx=(plasma.p(:,1)==10 & plasma.f(1,:)==1e6);
figure
% 径向
subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
plot(r_line,ne_r_fit) 
hold on
temp_color=co(2,:);
line([0,45.5],[dist_rp(1).ne_r(1) dist_rp(1).ne_r(1)],'linestyle','--','color',temp_color);
axis([0,45.5,0,1]) 
xlabel('{\itr} [mm]')
xticks([0,20,45.5])
ylabel('normalized {\itn}_e')
grid on

temp_color=co(3,:);
line([0,45.5],[1,1],'linestyle','-.','color',temp_color, 'LineWidth',3);

temp_color=co(4,:);
for i=3:3
    r=0:45.5/i:45.5;
    ne_r=[dist_rp(i).ne_r; 0];
    for j=1:i
        line([r(j),r(j+1)],[ne_r(j),ne_r(j)],'Color',temp_color,'linestyle',':')
        line([r(j+1),r(j+1)],[ne_r(j),ne_r(j+1)],'Color',temp_color,'linestyle',':')
    end
end
L1=legend('origin','uniform averaged','uniform central','nonuniform');
set(L1,'location','southwest');
set(L1,'box','off')
text(45.5/3,0.95,'(a)')
% PER
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
legend_text={};
temp_color=co(4,:);
plot(1:5,fem.nonuniform_in4_PER(2:end),'--o','Color',temp_color)
legend_text{end+1}='nonuniform';
hold on
scatter(1,fem.nonuniform_in4_PER(2),'s','MarkerEdgeColor',co(2,:),'LineWidth',4)
legend_text{end+1}='uniform averaged';
scatter(1,fem.nonuniform_in4_PER(1),'d','MarkerEdgeColor',co(3,:),'LineWidth',2)
legend_text{end+1}='uniform central';
ylabel('PER [\Omega]')
xlabel('discrete parts')
grid on
line([1,5],[experiment.PER(idx) experiment.PER(idx)],'linestyle','--','color','k')
L1=legend(legend_text);
set(L1,'location','northeast');
set(L1,'box','off')
min_y=0;
max_y=1;
axis([1,5,min_y,max_y]) 
text(2.5,0.95*max_y,'(b)')
text(1.3,0.1,'experiment measured')


%% 讨论直流电导率应用条件
% sigma_dc的适用条件: nu_c>>omega
% PER
figure
%1MHz
subplot(1,2,1);
set(gca, 'Position',[0.12 0.19 0.355 0.73])
legend_text={};
semilogx(p,experiment.PER(:,1),'-o')
legend_text{end+1}='experiment';
hold on
semilogx(p,fem.dielectric_PER(:,1),'-s')
legend_text{end+1}='as dielectric';
set(gca,'ColorOrderIndex',4);
semilogx(p,fem.conductor_PER(:,1),'-^')
legend_text{end+1}='as conductor';
ylabel('PER [\Omega]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','north');
set(L1,'box','off')
min_y=0;
max_y=3;
axis([0.3,10,min_y,max_y]) 
text(0.4,0.95*max_y,'(a) {\itf}=1MHz')
% 4MHz
subplot(1,2,2);
set(gca, 'Position',[0.61 0.19 0.355 0.73])
legend_text={};
semilogx(p,experiment.PER(:,2),'-.o')
legend_text{end+1}='experiment';
hold on
semilogx(p,fem.dielectric_PER(:,2),'-.s')
legend_text{end+1}='as dielectric';
set(gca,'ColorOrderIndex',4);
semilogx(p,fem.conductor_PER(:,2),'-^')
legend_text{end+1}='as conductor';
ylabel('PER [\Omega]')
xlabel('{\itp}_{fill} [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','north');
set(L1,'box','off')
min_y=0;
max_y=24;
axis([0.3,10,min_y,max_y]) 
text(0.4,0.95*max_y,'(b) {\itf}=4MHz')

%% discussion: ne、Te耦合的影响
% 以后再说
% 解析电磁、FEM model
% 以1MHz, 10Pa为例，假设ne、Te解耦且与其他参数解耦
% PER

% Rs, Ls
