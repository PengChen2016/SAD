% BUG_sweep case for paper
% 20210508 created by pengchen2016, Matlab R2017a

%% initialization
close all
clear
tic
addpath(genpath('./packages'))
now_str=datestr(now,'yyyymmdd_HHMMSS');
%% flag
solution_name='for_paper210415';
addpath(genpath(['./others/' solution_name '/']))
flag.using_stored_data=false;
% flag.using_stored_data=true;

%------------------------------ 1. transformer-base
program_name='paper_BUG_edge_t210517';
%%%%%%%% plasma model
flag.type_Xsec='e-H2-Phelps';
flag.input_plasma='2021Zielke_BUG_sweep';
flag.stoc_model='Cazzador-fit';
flag.medium_approximation='';
flag.skin_depth='as-medium';
flag.output_plasma_model=true;
%%%%%%%% electric model
flag.electric_model='transformer-base';
flag.input_geometry='BUG_base';
flag.Rmetal='measured-Rmetal-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
flag.output_electric_model=true;
%------------------------------ 1. flag end

% %------------------------------ 2. raza, analytical-base
% program_name='paper_BUG_edge_a210517';
% %%%%%%%% plasma model
% flag.type_Xsec='e-H2-Phelps';
% flag.input_plasma='2021Zielke_BUG_sweep';
% flag.stoc_model='Cazzador-fit';
% flag.medium_approximation='';
% flag.skin_depth='as-medium';
% flag.output_plasma_model=true;
% %%%%%%%% electric model
% flag.electric_model='analytical_base';
% flag.input_geometry='BUG_base';
% flag.Rmetal='measured-Rmetal-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';
% flag.output_electric_model=true;
% flag.magnetized='TODO';
% %------------------------------ 2. flag end

% %------------------------------ 3. raza, transformer-2018Jainb
% program_name='paper_BUG_edge_tj210517';
% %%%%%%%% plasma model
% flag.type_Xsec='e-H2-Phelps';
% flag.input_plasma='2021Zielke_BUG_sweep';
% flag.stoc_model='2018Jainb-simplify';
% flag.medium_approximation='';
% flag.skin_depth='as-medium-simplified-finite-radius';
% flag.output_plasma_model=true;
% %%%%%%%% electric model
% flag.electric_model='transformer-2018Jainb';
% flag.input_geometry='BUG_base';
% flag.Rmetal='measured-Rmetal-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';
% flag.output_electric_model=true;
% %------------------------------ 3. flag end


% % program_name='BUG_edge_sweep210508';
% program_name='BUG_analyticalEM210511';
% % program_name='BUG_transformerJ210511';
% 
% %%%%%%%% plasma model
% flag.type_Xsec='e-H2-Phelps';
% 
% flag.input_plasma='2021Zielke_BUG_sweep';
% 
% % stoc expression
% % flag.stoc_model='';
% flag.stoc_model='Vahedi-simplify';
% % flag.stoc_model='Cazzador-simplify';
% % flag.stoc_model='2018Jainb-simplify';
% % flag.stoc_model='Cazzador-fit';
% 
% flag.medium_approximation='';
% % flag.medium_approximation='consider_real(sigma)_only';
% % flag.medium_approximation='sigma_dc';
% 
% % flag.skin_depth='as-medium';
% flag.skin_depth='as-medium-simplified';
% % flag.skin_depth='collisionless';
% % flag.skin_depth='collision';
% % flag.skin_depth='as-medium-simplified-finite-radius';
% 
% flag.output_plasma_model=true;
% % flag.output_plasma_model=false;
% 
% %%%%%%%% electric model
% % flag.electric_model='';
% % flag.electric_model='transformer-base';
% % flag.electric_model='transformer-2011Chabert';
% % flag.electric_model='transformer-2015Bandyopadhyay';
% % flag.electric_model='transformer-2018Jainb';
% flag.electric_model='analytical_base';
% 
% flag.input_geometry='BUG_base';
% 
% flag.Rmetal='measured-Rmetal-woplasma';
% % flag.Rmetal=calculated-Rcoil-woplasma';
% % flag.Lcoil='measured-Lcoil-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';
% 
% flag.output_electric_model=true;
% % flag.output_electric_model=false;

%% solving
save_mat_name=['./others/' solution_name '/' program_name '.mat'];
if ~flag.using_stored_data
    % calculate
    log_name=['./others/' solution_name '/' program_name '.log'];
    diary(log_name) % append to the end of the log file
    fprintf('\n-----%s %s-----\n\n',program_name,now_str)

    input=get_input_data( flag );
    
    input.plasma=plasma_model(flag, input.plasma);
    
    if isfield(flag,'electric_model') && ~isempty(flag.electric_model)
        source=electric_model( flag, input );
    end
    
    % magnetized
    if isfield(flag,'magnetized') && strcmp(flag.magnetized,'TODO')
        input_m=input;
        source_ma=analytical_EM_model( flag, input_m );
        for i=1:4
            [input_m, source_ma]=get_magnetized(flag,input_m,source_ma);
        end
        if isfield(flag,'electric_model') && ~isempty(flag.electric_model)
            source_m=electric_model( flag, input_m );
        end
    end
    
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
% check input
% get_BUG_resultsmat();
load('BUG_results.mat')
d=@(mat,ip) reshape(mat(ip,1,:),1,3);

% ne/Te
figure
plot(d(experiment.Pplasma,1),d(input.plasma.ne,1),'--xy')
hold on
plot(d(experiment.Pplasma,2),d(input.plasma.ne,2),'--xg')
plot(d(experiment.Pplasma,3),d(input.plasma.ne,3),'--xc')
axis([0,50,0,1e18])
xlabel('{\itP}_{plasma} [kW]');
ylabel('n_e [m^{-3}]');
L1=legend('p=0.2Pa','p=0.3Pa','p=0.5Pa');
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格

figure
plot(d(experiment.Pplasma,1),d(input.plasma.Te,1),'--xy')
hold on
plot(d(experiment.Pplasma,2),d(input.plasma.Te,2),'--xg')
plot(d(experiment.Pplasma,3),d(input.plasma.Te,3),'--xc')
axis([0,50,0,16])
xlabel('{\itP}_{plasma} [kW]');
ylabel('T_e [eV]');
L1=legend('p=0.2Pa','p=0.3Pa','p=0.5Pa');
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格

if ~flag.using_stored_data
    fprintf('\n-----END %s-----\n\n',now_str)
    diary off
end

%% disp
% o=@(mat) mat([1,4,7,2,5,8,3,6,9]');
% p=input.plasma;
% o(p.sigma)
% o(p.eps_r)
% o(p.sigma_dc)
% o(p.skin_depth)