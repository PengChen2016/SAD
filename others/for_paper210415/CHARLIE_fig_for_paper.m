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
% 效率：PER，PTE
% 阻抗：Rs，Ls

clear
close all
constants=get_constants();
solution_name='for_paper210415';
% pwd == the root of SAD code
addpath(genpath(['./others/' solution_name '/']))
addpath(genpath('./packages'))

% - 1MHz -. 4MHz  --k 辅助线
% or / sb / dc 物理量
% TODO: 批量修改线型
lstyle={'-','-.','--'};
lmarker={'o','s','d','x','+'};

q1='or';


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

p=[0.3, 0.5, 1, 3, 5, 10]';
w1MHz=2*pi*1e6;
w4MHz=2*pi*4e6;

%% input
% 考虑ne不均匀分布，估算平均值/边缘值
norm_ne_r10=nonuniform_dist.get_ne_r(10); % origin data from experiments
for i=1:5
    dist_rp(i)=nonuniform_dist.get_nonuniform_dist_CHARLIE(['rp' num2str(i)]);
    ratio_origin2goal(i).ne_r=dist_rp(i).ne_r/norm_ne_r10;
end
data_ref=nonuniform_dist.get_ref_CHARLIE();
r_line=(0:45.5/100:45.5)';
ne_r_fit=nonuniform_dist.get_ne_r(r_line);
Te_r_fit=nonuniform_dist.get_Te_r(r_line);

z_line=(-200:10:200)';
ne_z_fit=nonuniform_dist.get_ne_z(z_line);
Te_z_fit=nonuniform_dist.get_Te_z(z_line);
mean_ne_z=nonuniform_dist.get_ne_z([-200,200]);
mean_coil_ne_z=nonuniform_dist.get_ne_z([-50,50]);

% ne 分布
figure
% 径向
subplot(1,2,1);
% scatter(data_ref.ne_r(:,1),data_ref.ne_r(:,2),'o','MarkerEdgeColor','k')
plot(r_line,ne_r_fit,'-r')
hold on 
line([0,45.5],[dist_rp(1).ne_r(1) dist_rp(1).ne_r(1)],'linestyle','-','color','b');
% scatter(data_ref.ne_r(:,1),data_ref.ne_r(:,2),'o','MarkerEdgeColor','k')
axis([0,45.5,0,1]) 
xlabel('{radial position} [mm]')
xticks([0,20,45.5])
ylabel('normalized n_{\rme}')
grid on
line([10,10],[0.8, 1],'linestyle','--','color','k','LineWidth',1);
text(5,1.01,'diagnostic port')
L1=legend('PIC/MCC','mean');
set(L1,'location','southwest');
set(L1,'box','off')
text(-20,0,'(a)')
% 轴向
subplot(1,2,2);
% scatter(data_ref.ne_z(:,1),data_ref.ne_z(:,2),'o','MarkerEdgeColor','k')
plot(z_line,ne_z_fit,'-r')
hold on 
line([-200,200],[mean_ne_z, mean_ne_z],'linestyle','-','color','b');
line([-200,200],[mean_coil_ne_z, mean_coil_ne_z],'linestyle','-.','color','m');
% scatter(data_ref.ne_z(:,1),data_ref.ne_z(:,2),'o','MarkerEdgeColor','k')
axis([-200,200,0,1]) 
xlabel('{axial position} [mm]')
xticks(-200:200:200)
ylabel('normalized n_{\rme}')
grid on
line([-50,-50],[0.5, 1],'linestyle','--','color','k','LineWidth',1);
line([50,50],[0.5, 1],'linestyle','--','color','k','LineWidth',1);
text(-50,0.6,'coil')
L1=legend('experiment','mean1','mean2');
set(L1,'location','southwest');
set(L1,'box','off')
text(-360,0,'(b)')
% set(gca,'yaxislocation','right');

% ne, Te
figure
% ne
subplot(1,2,1);
legend_text={};
semilogx(p,plasma.ne(:,1),'-or')
legend_text{end+1}='1MHz';
hold on
semilogx(p,plasma.ne(:,2),'-.or')
legend_text{end+1}='4MHz';
axis([0.3,10,0,2e17]) 
ylabel('{n}_e [m^{-3}]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southeast');
set(L1,'box','off')
text(0.4,1.9e17,'(a)')
% Te
subplot(1,2,2);
legend_text={};
semilogx(p,plasma.Te(:,1),'-sb')
legend_text{end+1}='1MHz';
hold on
semilogx(p,plasma.Te(:,2),'-.sb')
legend_text{end+1}='4MHz';
axis([0.3,10,0,6]) 
ylabel('{T}_e [eV]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southeast');
set(L1,'box','off')
text(0.4,5.7,'(b)')

% 等离子体模型输出
% nu_c
% semilogx
figure
% 1MHz
subplot(1,2,1);
legend_text={};
semilogx(p,plasma.nu_m(:,1),'-or')
legend_text{end+1}='\nu_m';
hold on
semilogx(p,plasma.nu_st(:,1),'-sb')
legend_text{end+1}='\nu_{st}';
axis([0.3,10,0,1.3e8]) 
ylabel('frequency [s^{-1}]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
line([0.3,10],[w1MHz, w1MHz],'linestyle','--','color','k');
legend_text{end+1}='\omega';
L1=legend(legend_text);
set(L1,'location','northwest');
set(L1,'box','off')
text(0.35,1.21e8,'(a) f=1MHz')
% 4MHz
subplot(1,2,2);
legend_text={};
semilogx(p,plasma.nu_m(:,2),'-.or')
legend_text{end+1}='\nu_m';
hold on
semilogx(p,plasma.nu_st(:,2),'-.sb')
legend_text{end+1}='\nu_{st}';
axis([0.3,10,0,1.3e8]) 
ylabel('frequency [s^{-1}]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
line([0.3,10],[w4MHz, w4MHz],'linestyle','--','color','k');
legend_text{end+1}='\omega';
L1=legend(legend_text);
set(L1,'location','northwest');
set(L1,'box','off')
text(0.35,1.21e8,'(b) f=4MHz')
% loglog 与其他图semilogx不一致，且难以突出不同物理量数量对比
% 但有利于体现nu_m与p成指数关系
figure
% 1MHz
subplot(1,2,1);
legend_text={};
loglog(p,plasma.nu_m(:,1),'-or')
legend_text{end+1}='\nu_m';
hold on
loglog(p,plasma.nu_st(:,1),'-sb')
legend_text{end+1}='\nu_{st}';
axis([0.3,10,3e6,2e8]) 
ylabel('frequency [s^{-1}]')
% yticks([3e6,1e7,1e8,2e8])
% yticklabels({'3\times10^6','10^7','10^8','2\times10^8'})
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
line([0.3,10],[w1MHz, w1MHz],'linestyle','--','color','k');
legend_text{end+1}='\omega';
L1=legend(legend_text);
set(L1,'location','northwest');
set(L1,'box','off')
text(0.35,1.6e8,'(a) f=1MHz')
% 4MHz
subplot(1,2,2);
legend_text={};
loglog(p,plasma.nu_m(:,2),'-.or')
legend_text{end+1}='\nu_m';
hold on
loglog(p,plasma.nu_st(:,2),'-.sb')
legend_text{end+1}='\nu_{st}';
axis([0.3,10,3e6,2e8]) 
ylabel('frequency [s^{-1}]')
% yticks([3e6,1e7,1e8,2e8])
% yticklabels({'3\times10^6','10^7','10^8','2\times10^8'})
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
line([0.3,10],[w4MHz, w4MHz],'linestyle','--','color','k');
legend_text{end+1}='\omega';
L1=legend(legend_text);
set(L1,'location','northwest');
set(L1,'box','off')
text(0.35,1.6e8,'(b) f=4MHz')

% nu_m中nu_enp占绝大多数，没必要画图

% sigma, eps_r
figure
% sigma
subplot(1,2,1);
legend_text={};
semilogx(p,plasma.sigma(:,1),'-or')
legend_text{end+1}='1MHz';
hold on
semilogx(p,plasma.sigma(:,2),'-.or')
legend_text{end+1}='4MHz';
axis([0.3,10,0,1e2]) 
ylabel('\sigma [S/m]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southwest');
set(L1,'box','off')
text(5,93,'(a)')
% 4MHz
subplot(1,2,2);
legend_text={};
semilogx(p,plasma.eps_r(:,1),'-sb')
legend_text{end+1}='1MHz';
hold on
semilogx(p,plasma.eps_r(:,2),'-.sb')
legend_text{end+1}='4MHz';
axis([0.3,10,-7e5,0]) 
ylabel('{\bf\epsilon}_r')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southeast');
set(L1,'box','off')
text(0.35,-0.5e5,'(b)')


%%%%%%%%%% main: different electric models
% FEM model，解析电磁，变压器model，实验
% PER
%1MHz
h_fig=plot_nY(p,'{p} [Pa]',...
    experiment.PER(:,1),'subtractive method',...
    fem.dielectric_PER(:,1),'FEM model',...
    source_t.PER(:,1),'transformer model',...
    source_a.PER(:,1),'analytical EM model',...
    'PER [\Omega]','plot');
for i=1:length(h_fig)
    figure(h_fig{i})
    title('1MHz')
%     yticks([0.1,0.2,0.5,1,2,5,10])
%     axis([0.3,10,0.1,10])
end
% yticks([0.1,0.2,0.5,1,2,3]) % for loglog
xticks('auto') % for plot

% 4MHz
h_fig=plot_nY(p,'{p} [Pa]',...
    experiment.PER(:,2),'subtractive method',...
    fem.dielectric_PER(:,2),'FEM model',...
    source_t.PER(:,2),'transformer model',...
    source_a.PER(:,2),'analytical EM model',...
    'PER [\Omega]','plot');
for i=1:length(h_fig)
    figure(h_fig{i})
    title('4MHz')
%     yticks([0.1,0.2,0.5,1,2,5,10])
%     axis([0.3,10,0.1,10])
end
xticks('auto') % for plot

% Rs, Ls
% 1MHz

% 4MHz

%%%%% 误差分析
% 磁化与离子动力学的影响
% 回旋频率
constants=get_constants();
% H
Hz_edge=source_a.emf.Hz_plasma(45.5e-3);
Hz_half=source_a.emf.Hz_plasma(0.5*45.5e-3);
Hz_center=source_a.emf.Hz_plasma(0);
% Bm
Bz_edge=constants.mu0*abs(Hz_edge);
Bz_half=constants.mu0*abs(Hz_half);
Bz_center=constants.mu0*abs(Hz_center);
w_ce_min=constants.e*Bz_center/constants.me;

figure
% wpe,wpi,w
subplot(1,3,1);
legend_text={};
loglog(p,plasma.wpe(:,1),'-or')
legend_text{end+1}='\omega_{pe}';
hold on
loglog(p,plasma.nu_st(:,1),'-sb')
legend_text{end+1}='\omega_{pi}';
axis([0.3,10,5e6,3e10]) 
ylabel('frequency [s^{-1}]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
line([0.3,10],[w1MHz, w1MHz],'linestyle','-','color','k');
legend_text{end+1}='\omega, 1MHz';
line([0.3,10],[w4MHz, w4MHz],'linestyle','-.','color','k');
legend_text{end+1}='\omega, 4MHz';
L1=legend(legend_text);
set(L1,'location','west');
set(L1,'box','off')
text(0.35,2.2e10,'(a)')
% nu_eff, wce
% 1MHz
subplot(1,3,2);
legend_text={};
loglog(p,plasma.nu_eff(:,1),'-or')
legend_text{end+1}='\nu_{c}';
hold on
loglog(p,w_ce_min(:,1),'-sb')
legend_text{end+1}='\omega_{ce}';
axis([0.3,10,1e7,5e8]) 
ylabel('frequency [s^{-1}]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','best');
set(L1,'box','off')
text(0.8,1.6e8,'(b) f=1MHz')
subplot(1,3,3);
legend_text={};
loglog(p,plasma.nu_eff(:,2),'-.or')
legend_text{end+1}='\nu_{c}';
hold on
loglog(p,w_ce_min(:,2),'-.sb')
legend_text{end+1}='\omega_{ce}';
% axis([0.3,10,3e6,2e8]) 
ylabel('frequency [s^{-1}]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','northwest');
set(L1,'box','off')
text(0.35,1.6e8,'(c) f=4MHz')

% 变压器模型不适用-密度不够高 delta
figure
% sigma
subplot(1,2,1);
legend_text={};
semilogx(p,plasma.sigma(:,1),'-or')
legend_text{end+1}='1MHz';
hold on
semilogx(p,plasma.sigma(:,2),'-.or')
legend_text{end+1}='4MHz';
axis([0.3,10,0,1e2]) 
ylabel('\sigma [S/m]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southwest');
set(L1,'box','off')
text(5,93,'(a)')
% 4MHz
subplot(1,2,2);
legend_text={};
semilogx(p,plasma.eps_r(:,1),'-sb')
legend_text{end+1}='1MHz';
hold on
semilogx(p,plasma.eps_r(:,2),'-.sb')
legend_text{end+1}='4MHz';
axis([0.3,10,-7e5,0]) 
ylabel('{\bf\epsilon}_r')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southeast');
set(L1,'box','off')
text(0.35,-0.5e5,'(b)')

%%%%%%%%%% discussion: 有损介质的影响
% sigma/(-w*eps), sigma,sigma_dc的适用条件是什么？
figure
% sigma
subplot(1,2,1);
legend_text={};
plot(p,-w1MHz*plasma.eps_prime(:,1),'-.ob')
legend_text{end+1}='-\omega\eps';
hold on
semilogx(p,plasma.sigma(:,2),'-.or')
legend_text{end+1}='4MHz';
axis([0.3,10,0,1e2]) 
ylabel('\sigma [S/m]')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southwest');
set(L1,'box','off')
text(5,93,'(a)')
% 4MHz
subplot(1,2,2);
legend_text={};
semilogx(p,plasma.eps_r(:,1),'-sb')
legend_text{end+1}='1MHz';
hold on
semilogx(p,plasma.eps_r(:,2),'-.sb')
legend_text{end+1}='4MHz';
axis([0.3,10,-7e5,0]) 
ylabel('{\bf\epsilon}_r')
xlabel('filling pressure [Pa]')
xticks([0.3,1,3,10])
grid on
L1=legend(legend_text);
set(L1,'location','southeast');
set(L1,'box','off')
text(0.35,-0.5e5,'(b)')

% FEM model，实验 结果
% PER
% 1MHz
h_fig=plot_nY(p,'{p} [Pa]',...
    experiment.PER(:,1),'subtractive method',...
    fem.dielectric_PER(:,1),'FEM model',...
    fem.conductor_PER(:,1),'FEM model',...
    'PER [\Omega]','plot');
title('1MHz')

% 4MHz
h_fig=plot_nY(p,'{p} [Pa]',...
    experiment.PER(:,2),'subtractive method',...
    fem.dielectric_PER(:,2),'FEM model',...
    fem.conductor_PER(:,2),'FEM model',...
    'PER [\Omega]','plot');
title('4MHz')
% Rs, Ls
% 1MHz

% 4MHz

%%%%%%%%%% discussion: 不均匀ne的影响
% 先简单展示一下，说明有必要考虑，以后再定量对比
figure
scatter(data_ref.ne_r(:,1),data_ref.ne_r(:,2),'o','MarkerEdgeColor','k')
hold on 
plot(r_line,ne_r_fit,'-y')
line_color={'r','g','b','c','m','y','k'};
for i=1:5
    r=0:45.5/i:45.5;
    ne_r=[dist_rp(i).ne_r; 0];
    for j=1:i
        line([r(j),r(j+1)],[ne_r(j),ne_r(j)],'Color',line_color{i},'linestyle','-.')
        line([r(j+1),r(j+1)],[ne_r(j),ne_r(j+1)],'Color',line_color{i},'linestyle','-.')
    end
end
axis([0,45.5,0,1]) 
xlabel('{r} [mm]')
ylabel('Normalized n_{\rme}')
grid on
L1=legend('PIC/MCC','Fit','Uniform case');
set(L1,'location','southwest');
set(L1,'box','off')

% FEM model 分层
% 以1MHz, 10Pa为例
% PER

% Rs, Ls

%%%%%%%%%% discussion: ne、Te耦合的影响
% 以后再说
% 解析电磁、FEM model
% 以1MHz, 10Pa为例，假设ne、Te解耦且与其他参数解耦
% PER

% Rs, Ls


