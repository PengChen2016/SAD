% BUG_nonuniform case for paper
% 20210716 created by pengchen2016, Matlab R2017a
% eP-210607-01 双驱串并联融合\

%% initialization
close all
clear
tic
addpath(genpath('./packages'))
now_str=datestr(now,'yyyymmdd_HHMMSS');
%% flag
solution_name='eP-210607-01';
addpath(genpath(['./others/' solution_name '/']))
% flag.using_stored_data=false;
flag.using_stored_data=true;

% program_name='BUG_z_nonuniform210716';
program_name='BUG_z_nonuniform220304';
%%%%%%%% plasma model
flag.type_Xsec='e-H2-Phelps';
flag.input_plasma='BUG_nonuniform_from_Tandem';
% stoc expression
% flag.stoc_model='';
% flag.stoc_model='Vahedi-simplify';
% flag.stoc_model='Cazzador-simplify';
% flag.stoc_model='2018Jainb-simplify';
flag.stoc_model='Cazzador-fit';
flag.medium_approximation='';
% flag.medium_approximation='consider_real(sigma)_only';
% flag.medium_approximation='sigma_dc';
flag.skin_depth='as-medium';
% flag.skin_depth='as-medium-simplified';
% flag.skin_depth='collisionless';
% flag.skin_depth='collision';
% flag.skin_depth='as-medium-simplified-finite-radius';
flag.output_plasma_model=true;
% flag.output_plasma_model=false;
%%%%%%%% electric model
flag.electric_model='';
flag.input_geometry='BUG_base';

%% solving
save_mat_name=['./others/' solution_name '/' program_name '.mat'];
if ~flag.using_stored_data
    % calculate
    log_name=['./others/' solution_name '/' program_name '.log'];
    diary(log_name) % append to the end of the log file
    fprintf('\n-----%s %s-----\n\n',program_name,now_str)

    input=get_input_data( flag );
    input.plasma=plasma_model(flag, input.plasma);
    
    save(save_mat_name)
else
    % load data
    load(save_mat_name)
    disp(now_str)
    warning(['Using data stored in ' save_mat_name])
    if flag.output_plasma_model
        output_plasma_model(flag,input.plasma)
    end
    if flag.output_electric_model
        output_electric_model( flag, source )
    end
end

%% post-processing
input.geometry=get_input_geometry( flag.input_geometry );
z=(0:0.01:1)*input.geometry.l_chamber; % 见input_data

figure
yyaxis left
plot(z,input.plasma.ne)
ylabel('{\itn}_{e} [m^{-3}]');
yyaxis right
plot(z,log10(input.plasma.ne))
ylabel('lg({\itn}_{e}) [m^{-3}]');
yyaxis left
axis([0,input.geometry.l_chamber,0,3.5e17])
xlabel('\itz [m]');
L1=legend('ne','lg(ne)');
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格

figure
yyaxis left
plot(z,input.plasma.sigma)
ylabel('{\it\sigma} [S/m]');
yyaxis right
plot(z,-input.plasma.eps_r)
ylabel('-{\it\epsilon}_{r}');
yyaxis left
axis([0,input.geometry.l_chamber,0,500])
xlabel('\itz [m]');
L1=legend('{\it\sigma}','-{\it\epsilon}_{r}');
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格

% fit
% see .\others\eP-210607-01\fit_Tandem_nonuniform_plasma.sfit

if ~flag.using_stored_data
    fprintf('\n-----END %s-----\n\n',now_str)
    diary off
end