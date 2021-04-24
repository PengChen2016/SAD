% analyze the result of CHARLIE_raza_sweep20210415.m
clear

variable_Y_plasma={'ne','nu_m','nu_st',...
    'skin_depth','sigma_dc','eps_r'};
name_Y_plasma={'n_e [m^{-3}]','\nu_m [Hz]','\nu_{st} [Hz]',...
    '\delta_{skin} [m]','\sigma_{dc} [S/m]','\epsilon_r'};
variable_Y_source={'PER','PTE','Xplasma',...
    'transform.Rp'};
name_Y_source={'R_{plasma} [\Omega]','\eta','X_{plasma} [\Omega]'...
    ,'R_p'};

%% 对比ne
% compare CHARLIE_raza_sweep210423 and CHARLIE_r10za_sweep210423
close all
load('CHARLIE_raza_sweep210423.mat')
plasma_raza=input.plasma;
source_raza=source;

load('CHARLIE_r10za_sweep210423.mat')
plasma_r10za=input.plasma;
source_r10za=source;

for i=1:length(variable_Y_plasma)
plot_2Ylog(eval(['plasma_raza.' variable_Y_plasma{i}]), 'raza',...
    eval(['plasma_r10za.' variable_Y_plasma{i}]), 'r10za',...
    name_Y_plasma{i});
end
for i=1:length(variable_Y_source)
plot_2Ylog(eval(['source_raza.' variable_Y_source{i}]), 'raza',...
    eval(['source_r10za.' variable_Y_source{i}]), 'r10za',...
    name_Y_source{i});
end

%% 对比集肤深度公式
% compare CHARLIE_raza_sweep210423 and CHARLIE_raza_vahedi_skin_depth210424
close all
load('CHARLIE_raza_vahedi_skin_depth210424.mat')
plasma_raza_smaller_delta=input.plasma;
source_raza_smaller_delta=source;

for i=1:length(variable_Y_plasma)
plot_2Ylog(eval(['plasma_raza.' variable_Y_plasma{i}]), 'semi-infinite',...
    eval(['plasma_raza_smaller_delta.' variable_Y_plasma{i}]), 'finite-size',...
    name_Y_plasma{i});
end
for i=1:length(variable_Y_source)
plot_2Ylog(eval(['source_raza.' variable_Y_source{i}]), 'semi-infinite',...
    eval(['source_raza_smaller_delta.' variable_Y_source{i}]), 'finite-size',...
    name_Y_source{i});
end