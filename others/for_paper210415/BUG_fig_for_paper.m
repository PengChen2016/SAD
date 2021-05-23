% BUG fig for paper
% 说明有仿真聚变负源的能力，显示其优点
clear
close all

constants=get_constants();

load('BUG_results.mat')
load('paper_BUG_edge_t210517.mat')
plasma=input.plasma;
source_t=source;
load('paper_BUG_edge_a210517.mat')
source_a=source;
plasma_ma=input_m.plasma;
load('paper_BUG_edge_tj210517.mat')
plasma_j=input.plasma;
source_j=source;

p=[0.3, 0.5, 1, 3, 5, 10]';
w1MHz=2*pi*1e6;

% PER

% Rs, Ls

%% magnetized
% intial
size_mat=size(source_ma.PER);
r1=0;
r2=input.geometry.r_plasma_eff;
r3=input.geometry.r_coil;
r_a=[r1:(r2-r1)/100:r2 r2+(r3-r2)/10:(r3-r2)/10:r3];
% 前101个在r_plasma内，后10个在>r_plasma
len_r=length(r_a);

B_a=zeros(1,len_r);
B_ma=zeros(1,len_r);
idx=(plasma.p==0.3 & plasma.Pin==18e3);
for i_r=1:len_r
    temp=source_a.emf.Bzm_r(r_a(i_r));
    B_a(i_r)=temp(idx);
    temp=source_ma.emf.Bzm_r(r_a(i_r));
    B_ma(i_r)=temp(idx);
end

% B_rcoil_a=B_a(end);
% B_rplasma_a=B_a(101);
% B_mean_plasma=input_m.plasma.Bz_mean(idx);

figure
plot(1e3*fem.radialHztypical_r,constants.mu0*fem.radialHztypical_Hz,'-')
hold on
plot(1e3*r_a,B_a,'-.')
plot(1e3*r_a,B_ma,'--')
grid on
xlabel('{radial position} [mm]')
ylabel('B_{\rmz}')
grid on
% line([-50,-50],[0.5, 1],'linestyle','--','color','k','LineWidth',1);
% line([50,50],[0.5, 1],'linestyle','--','color','k','LineWidth',1);
% text(-50,0.6,'coil')
L1=legend('FEM','analytical','analytical-magnetized');
set(L1,'location','southwest');
set(L1,'box','off')
% text(-360,0,'(b)')
