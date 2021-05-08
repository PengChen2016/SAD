% analyze the result of CHARLIE_raza_sweep20210415.m

clear
solution_name='for_paper210415';
% pwd == the root of DSA
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
load('experiment_fem_results.mat')
load('CHARLIE_raza_sweep210423.mat')
plasma_raza=input.plasma;
source_raza=source;
plasma_t=input.plasma;
source_t=source;
load('CHARLIE_raza_sweep_analyticalEM210508.mat')
plasma_a=input.plasma;
source_a=source;

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

%1MHz 对比实验与模型（标准）的PER
h_fig=plot_nY(experiment.PER(:,1),'subtractive method',...
    fem.dielectric_PER(:,1),'FEM model',...
    source_t.PER(:,1),'transformer model',...
    source_a.PER(:,1),'analytical EM model',...
    'PER [\Omega]','all');
for i=1:length(h_fig)
    figure(h_fig{i})
    title('1MHz')
%     yticks([0.1,0.2,0.5,1,2,5,10])
%     axis([0.3,10,0.1,10])
end
% yticks([0.1,0.2,0.5,1,2,3]) % for loglog

% 4MHz
h_fig=plot_nY(experiment.PER(:,2),'subtractive method',...
    fem.dielectric_PER(:,2),'FEM model',...
    source_t.PER(:,2),'transformer model',...
    source_a.PER(:,2),'analytical EM model',...
    'PER [\Omega]','all');
for i=1:length(h_fig)
    figure(h_fig{i})
    title('4MHz')
%     yticks([0.1,0.2,0.5,1,2,5,10])
%     axis([0.3,10,0.1,10])
end

h_fig=plot_nY(plasma_raza.ne(:,1),'1MHz',...
    plasma_raza.ne(:,2),'4MHz',...
    'n_e [m^{-3}]','loglog');


%% FEM nonuniform
% results



%% other
close all
load('CHARLIE_raza_sweep210423.mat')
Rcoil_th=input.external.Rcoil_th(1,:);
l_wire=150e-3*2; % 单根引线长150mm
w_RF=input.plasma.w_RF(1,:);
constants=get_constants();
C_wire=2*pi*input.geometry.r_wire;
delta_Cu=sqrt(2./(w_RF*constants.mu0*constants.sigma_Cu));
S_current_path=delta_Cu.*C_wire;
Rwire_th=l_wire./(constants.sigma_Cu.*S_current_path);
Rmetal_th=Rcoil_th+Rwire_th;

Lcoil_th=1e6*input.external.Lcoil_th;
% 单匝引线的l极小，电感可忽略


