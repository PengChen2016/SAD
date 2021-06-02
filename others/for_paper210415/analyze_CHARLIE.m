% analyze the result of CHARLIE_sweep20210415.m

clear
solution_name='for_paper210415';
% pwd == the root of SAD code
addpath(genpath(['./others/' solution_name '/']))
addpath(genpath('./packages'))

variable_Y_plasma={'ne','nu_m','nu_st',...
    'skin_depth','sigma_dc','eps_r'};
name_Y_plasma={'n_e [m^{-3}]','\nu_m [Hz]','\nu_{st} [Hz]',...
    '\delta_{skin} [m]','\sigma_{dc} [S/m]','\epsilon_r'};
variable_Y_source={'PER','PTE','Xplasma',...
    'transform.Rp'};
name_Y_source={'R_{plasma} [\Omega]','\eta','X_{plasma} [\Omega]'...
    ,'R_p'};

% %% 对比ne
% % compare CHARLIE_raza_sweep210423 and CHARLIE_r10za_sweep210423
% close all
% load('CHARLIE_raza_sweep210423.mat')
% plasma_raza=input.plasma;
% source_raza=source;
% 
% load('CHARLIE_r10za_sweep210423.mat')
% plasma_r10za=input.plasma;
% source_r10za=source;
% 
% for i=1:length(variable_Y_plasma)
% plot_2Ylog(eval(['plasma_raza.' variable_Y_plasma{i}]), 'raza',...
%     eval(['plasma_r10za.' variable_Y_plasma{i}]), 'r10za',...
%     name_Y_plasma{i});
% end
% for i=1:length(variable_Y_source)
% plot_2Ylog(eval(['source_raza.' variable_Y_source{i}]), 'raza',...
%     eval(['source_r10za.' variable_Y_source{i}]), 'r10za',...
%     name_Y_source{i});
% end
% 
% %% 对比集肤深度公式
% % compare CHARLIE_raza_sweep210423 and CHARLIE_raza_vahedi_skin_depth210424
% close all
% load('CHARLIE_raza_vahedi_skin_depth210424.mat')
% plasma_raza_smaller_delta=input.plasma;
% source_raza_smaller_delta=source;
% 
% for i=1:length(variable_Y_plasma)
% plot_2Ylog(eval(['plasma_raza.' variable_Y_plasma{i}]), 'semi-infinite',...
%     eval(['plasma_raza_smaller_delta.' variable_Y_plasma{i}]), 'finite-size',...
%     name_Y_plasma{i});
% end
% for i=1:length(variable_Y_source)
% plot_2Ylog(eval(['source_raza.' variable_Y_source{i}]), 'semi-infinite',...
%     eval(['source_raza_smaller_delta.' variable_Y_source{i}]), 'finite-size',...
%     name_Y_source{i});
% end

%% FEM sweep
% results
% get_CHARLIE_resultsmat();
load('CHARLIE_results.mat')
load('CHARLIE_raza_sweep210423.mat')
plasma_raza=input.plasma;
source_raza=source;
plasma_t=input.plasma;
source_t=source;
load('CHARLIE_raza_sweep_analyticalEM210508.mat')
plasma_a=input.plasma;
source_a=source;
load('CHARLIE_raza_a_long210529.mat')
source_a_lplasma=source;

% 变压器模型的细节，不在本文范围之内。作为一个test去对比更合适
% load('CHARLIE_raza_vahedi_skin_depth210424.mat')
% plasma_raza_smaller_delta=input.plasma;
% source_raza_smaller_delta=source;
% load('CHARLIE_r10za_sweep210423.mat')
% plasma_r10za=input.plasma;
% source_r10za=source;
% load('CHARLIE_razcoil_sweep210428.mat')
% plasma_razcoil=input.plasma;
% source_razcoil=source;

p=[0.3, 0.5, 1, 3, 5, 10]';

h_fig=plot_nY(p,'{\itp} [Pa]',...
    plasma_raza.nu_enp(:,1),'\nu_{enp}, 1MHz',...
    plasma_raza.nu_eniz(:,1),'\nu_{eniz}, 1MHz',...
    plasma_raza.nu_eip(:,1),'\nu_{eip}, 1MHz',...
        plasma_raza.nu_enp(:,2),'\nu_{enp}, 4MHz',...
    plasma_raza.nu_eniz(:,2),'\nu_{eniz}, 4MHz',...
    plasma_raza.nu_eip(:,2),'\nu_{eip}, 4MHz',...
    'n_e [m^{-3}]','loglog');


h_fig=plot_nY(p,'{\itp} [Pa]',...
    plasma_raza.skin_depth(:,1),'\delta, 1MHz',...
    plasma_raza.wavelength(:,1),'\lambda, 1MHz',...
    plasma_raza.skin_depth(:,2),'\delta, 4MHz',...
    plasma_raza.wavelength(:,2),'\lambda, 4MHz',...
    'length [m]','semilogx');


% % 对比全部PER
% handle_fig=figure;
% legend_text={};
% 
% loglog(p,experiment.PER(:,1),'-xk');
% legend_text{end+1}='subtractive method, 1MHz';
% 
% hold on
% xticks(p)
% xlabel('{\itp} [Pa]');
% axis([0.3,10,-inf,inf])
% grid on
% 
% loglog(p,experiment.PER(:,2),'-.xk');
% legend_text{end+1}='subtractive method, 4MHz';
% loglog(p,fem.dielectric_PER(:,1),'-or','MarkerSize',8);
% legend_text{end+1}='FEM model, 1MHz';
% loglog(p,fem.dielectric_PER(:,2),'-.or','MarkerSize',8);
% legend_text{end+1}='FEM model, 4MHz';
% loglog(p,source_raza.PER(:,1),'-sb');
% legend_text{end+1}='transformer model, 1MHz';
% loglog(p,source_raza.PER(:,2),'-.sb');
% legend_text{end+1}='transformer model, 4MHz';
% ylabel('PER');
% yticks([0.1,0.2,0.5,1,2,5,10])
% 
% L1=legend(legend_text);
% set(L1,'Location','best');
% set(L1,'AutoUpdate','off');

h_fig=plot_nY(p,'{\itp} [Pa]',...
    plasma_raza.ne(:,1),'1MHz',...
    plasma_raza.ne(:,2),'4MHz',...
    'n_e [m^{-3}]','semilogy');
xticks('auto') % for plot
axis([0.3,10,-inf,inf])

%1MHz 对比实验与模型（标准）的PER
h_fig=plot_nY(p,'{\itp} [Pa]',...
    experiment.PER(:,1),'subtractive method',...
    fem.dielectric_PER(:,1),'FEM model',...
    source_t.PER(:,1),'transformer model',...
    source_a.PER(:,1),'analytical EM model',...
    source_a_lplasma.PER(:,1),'analytical EM model,long',...
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
h_fig=plot_nY(p,'{\itp} [Pa]',...
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

h_fig=plot_nY(p,'{\itp} [Pa]',...
    fem.dielectric_Rmetal(:,1),'FEM model-1MHz',...
    input.external.Rmetal_ex(:,1),'subtractive method-1MHz',...
    fem.dielectric_Rmetal(:,2),'FEM model-4MHz',...
    input.external.Rmetal_ex(:,2),'subtractive method-4MHz',...
    'R_{metal} [\Omega]','plot');
hold on
line([0.3,10],[input.external.Rcoil_th(1,1),input.external.Rcoil_th(1,1)],'linestyle','--','color','k');
line([0.3,10],[input.external.Rcoil_th(1,2),input.external.Rcoil_th(1,2)],'linestyle','--','color','k');
text(2,0.035,'Rcoil-analytical-1MHz')
text(2,0.06,'Rcoil-analytical-4MHz')

h_fig=plot_nY(p,'{\itp} [Pa]',...
    fem.dielectric_Ls(:,1),'FEM model-1MHz',...
    fem.dielectric_Ls(:,2),'FEM model-4MHz',...
    'L_{s} [\muH]','plot');
hold on
line([0.3,10],1e6*[input.external.Lcoil_ex,input.external.Lcoil_ex],'linestyle','--','color','k');
line([0.3,10],1e6*[input.external.Lcoil_th,input.external.Lcoil_th],'linestyle','--','color','k');
text(2,2.95,'Lcoil-experiment')
text(2,2.25,'Lcoil-experience')
xticks('auto') % for plot

%% lossy dielectric vs. good conductor
h_fig=plot_nY(p,'{\itp} [Pa]',...
    experiment.PER(:,1),'subtractive method',...
    fem.dielectric_PER(:,1),'FEM model',...
    fem.conductor_PER(:,1),'FEM model',...
    'PER [\Omega]','plot');
title('1MHz')

h_fig=plot_nY(p,'{\itp} [Pa]',...
    experiment.PER(:,2),'subtractive method',...
    fem.dielectric_PER(:,2),'FEM model',...
    fem.conductor_PER(:,2),'FEM model',...
    'PER [\Omega]','plot');
title('4MHz')


%% FEM nonuniform
% results
h_fig=plot_nY(1:5,'discrete parts',...
    fem.nonuniform_in4_PER(2:6),'non-uniform',...
    'PER [\Omega]','plot');
hold on
plot(1,fem.nonuniform_in4_PER(2),'-ob')
plot(1,fem.nonuniform_in4_PER(1),'-ok')
axis([1,5,0.2,inf])
legend('non-uniform case','mean n_e, uniform','n_e at r=0, uniform')

h_fig=plot_nY(1:5,'discrete parts',...
    fem.nonuniform_in4_PER(2:6),'non-uniform',...
    fem.nonuniform_in2_PER(2:6),'non-uniform',...
    'PER [\Omega]','plot');
hold on
scatter([1,5],[fem.nonuniform_in4_PER(2),fem.nonuniform_in4_ra_np5(1)],'d','MarkerEdgeColor','r','LineWidth',6)
scatter(1,fem.nonuniform_in2_PER(2),'d','MarkerEdgeColor','b','LineWidth',6)
scatter(1,fem.nonuniform_in4_PER(1),'s','MarkerEdgeColor','r','LineWidth',6)
scatter(1,fem.nonuniform_in2_PER(1),'s','MarkerEdgeColor','b','LineWidth',6)
axis([1,5,0.2,inf])
legend('non-uniform, in4 mesh','non-uniform, in2 mesh','mean n_e, uniform, in4 mesh','mean n_e, uniform, in2 mesh','n_e at r=0, uniform, in4 mesh','n_e at r=0, uniform, in2 mesh')

h_fig=plot_nY(1:5,'discrete parts',...
    fem.nonuniform_in4_Rmetal(2:6),'non-uniform',...
    fem.nonuniform_in2_Rmetal(2:6),'non-uniform',...
    'R_{metal} [\Omega]','plot');
hold on
scatter([1,5],[fem.nonuniform_in4_Rmetal(2),fem.nonuniform_in4_ra_np5(2)],'d','MarkerEdgeColor','r','LineWidth',6)
scatter(1,fem.nonuniform_in2_Rmetal(2),'d','MarkerEdgeColor','b','LineWidth',6)
scatter(1,fem.nonuniform_in4_Rmetal(1),'s','MarkerEdgeColor','r','LineWidth',8)
scatter(1,fem.nonuniform_in2_Rmetal(1),'s','MarkerEdgeColor','b','LineWidth',8)
line([1,5],[mean(input.external.Rmetal_ex(:,1)),mean(input.external.Rmetal_ex(:,1))],'linestyle','--','color','k');
axis([1,5,0,mean(input.external.Rmetal_ex(:,1))])
legend('non-uniform, in4 mesh','non-uniform, in2 mesh','mean n_e, uniform, in4 mesh','mean n_e, uniform, in2 mesh','n_e at r=0, uniform, in4 mesh','n_e at r=0, uniform, in2 mesh')
text(2,0.09,'mean Rmetal-experiment')

h_fig=plot_nY(1:5,'discrete parts',...
    fem.nonuniform_in4_Ls(2:6),'non-uniform',...
    fem.nonuniform_in2_Ls(2:6),'non-uniform',...
    'L_{s} [\muH]','plot');
hold on
scatter([1,5],[fem.nonuniform_in4_Ls(2),fem.nonuniform_in4_ra_np5(3)],'d','MarkerEdgeColor','r','LineWidth',6)
scatter(1,fem.nonuniform_in2_Ls(2),'d','MarkerEdgeColor','b','LineWidth',6)
scatter(1,fem.nonuniform_in4_Ls(1),'s','MarkerEdgeColor','r','LineWidth',8)
scatter(1,fem.nonuniform_in2_Ls(1),'s','MarkerEdgeColor','b','LineWidth',8)
line([1,5],1e6*[input.external.Lcoil_ex,input.external.Lcoil_ex],'linestyle','--','color','k');
line([1,5],1e6*[input.external.Lcoil_th,input.external.Lcoil_th],'linestyle','--','color','k');
text(2,2.95,'Lcoil-experiment')
text(2,2.25,'Lcoil-experience')
axis([1,5,-inf,1e6*input.external.Lcoil_th])
legend('non-uniform, in4 mesh','non-uniform, in2 mesh','mean n_e, uniform, in4 mesh','mean n_e, uniform, in2 mesh','n_e at r=0, uniform, in4 mesh','n_e at r=0, uniform, in2 mesh')

%% sewwp ne
ne=[10000000000000000, 15800000000000000, 25100000000000000, 39800000000000000, 63100000000000000, 100000000000000000, 117000000000000000, 158500000000000000, 243000000000000000, 251200000000000000, 398100000000000000, 631000000000000000, 1000000000000000000];
sigma=[2.2, 3.4, 5.4, 8.5, 13.3, 20.8, 24.3, 32.6, 49.3, 50.9, 79.1, 122.2, 187.5];

h_fig=plot_nY(ne,'n_e [m^{-3}]',...
    sigma,'\sigma',...
    '\sigma [S/m]','semilogx');
xticks('auto')

h_fig=plot_nY(ne,'n_e [m^{-3}]',...
    fem.sweep_ne_PER,'PER',...
    fem.sweep_ne_Rmetal,'Rmetal',...
    'R [\Omega]','semilogx');
xticks('auto')
axis([1e16,1e18,0,2.5])
yyaxis right
semilogx(ne,fem.sweep_ne_Ls,'-oc')
ylabel('Ls [\muH]')
axis([1e16,1e18,0,1e6*input.external.Lcoil_th])
legend('PER','Rmetal','Ls')

h_fig=plot_nY(sigma,'\sigma [S/m]',...
    fem.sweep_ne_PER,'PER',...
    fem.sweep_ne_Rmetal,'Rmetal',...
    'R [\Omega]','plot');
xticks('auto')
axis([1,max(sigma),0,2.5])
yyaxis right
plot(sigma,fem.sweep_ne_Ls,'-oc')
ylabel('Ls [\muH]')
axis([1,max(sigma),0,1e6*input.external.Lcoil_th])
legend('PER','Rmetal','Ls')

h_fig=plot_nY(ne,'n_e [m^{-3}]',...
    fem.sweep_ne_PER,'PER',...
    fem.sweep_ne_Rmetal,'Rmetal',...
    'R [\Omega]','plot');
xticks('auto')
axis([1e16,1e18,0,2.5])
yyaxis right
plot(ne,fem.sweep_ne_Ls,'-oc')
ylabel('Ls [\muH]')
axis([1e16,1e18,0,1e6*input.external.Lcoil_th])
legend('PER','Rmetal','Ls')

%% other
% close all
% load('CHARLIE_raza_sweep210423.mat')
% Rcoil_th=input.external.Rcoil_th(1,:);
% l_wire=150e-3*2; % 单根引线长150mm
% w_RF=input.plasma.w_RF(1,:);
% constants=get_constants();
% C_wire=2*pi*input.geometry.r_wire;
% delta_Cu=sqrt(2./(w_RF*constants.mu0*constants.sigma_Cu));
% S_current_path=delta_Cu.*C_wire;
% Rwire_th=l_wire./(constants.sigma_Cu.*S_current_path);
% Rmetal_th=Rcoil_th+Rwire_th;
% 
% Lcoil_th=1e6*input.external.Lcoil_th;
% % 单匝引线的l极小，电感可忽略

%% 磁化
% 基于解析EM模型做估算
input0=input;
input0.plasma=plasma_raza;

flag.output_plasma_model=true;
[input1, source1]=get_magnetized(flag,input0,source_a);

% 迭代
[input1, source1]=get_magnetized(flag,input1,source1);
[input1, source1]=get_magnetized(flag,input1,source1);
[input1, source1]=get_magnetized(flag,input1,source1);

% output
p1=input1.plasma;
h_fig=plot_nY(p,'{\itp} [Pa]',...
    p1.w_ce(:,1),'wce, 1MHz',...
    p1.w_ce(:,2),'wce, 4MHz',...
    p1.w_ci(:,1),'wci, 1MHz',...
    p1.w_ci(:,2),'wci, 4MHz',...
    p1.nu_eff(:,1),'\nu_{eff}, 1MHz',...
    p1.nu_eff(:,2),'\nu_{eff}, 4MHz',...
    p1.w_RF(:,1),'\omega, 1MHz',...
    p1.w_RF(:,2),'\omega, 4MHz',...
    '\omega [s^{-1}]','loglog');

h_fig=plot_nY(p,'{\itp} [Pa]',...
    experiment.PER(:,1),'subtractive method',...
    fem.dielectric_PER(:,1),'FEM model',...
    source_a.PER(:,1),'analytical EM-za',...
    source1.PER(:,1),'analytical EM-za-magnetized',...
    'PER [\Omega]','plot');
title('1MHz')
h_fig=plot_nY(p,'{\itp} [Pa]',...
    experiment.PER(:,2),'subtractive method',...
    fem.dielectric_PER(:,2),'FEM model',...
    source_a.PER(:,2),'analytical EM-za',...
    source1.PER(:,2),'analytical EM-za-magnetized',...
    'PER [\Omega]','plot');
title('4MHz')


p1.eps_r
p1.sigma

we=-p1.w_RF.*p1.eps_prime;
h_fig=plot_nY(p,'{\itp} [Pa]',...
    p1.sigma(:,1),'\sigma, 1MHz',...
    we(:,1),'-\omega\epsilon, 1MHz',...
    p1.sigma(:,2),'\sigma, 4MHz',...
    we(:,2),'-\omega\epsilon, 4MHz',...
    '','loglog');