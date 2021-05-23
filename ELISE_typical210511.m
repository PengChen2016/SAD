% ELISE_2018Jainb typical case for paper
% 20210511 created by pengchen2016, Matlab R2017a

%% initialization
close all
clear
tic
addpath(genpath('./packages'))
now_str=datestr(now,'yyyymmdd_HHMMSS');
%% flag
solution_name='for_paper210415';
program_name='ELISE_typical210511';

addpath(genpath(['./others/' solution_name '/']))
flag.using_stored_data=false;
% flag.using_stored_data=true;
save_mat_name=['./others/' solution_name '/' program_name '.mat'];
%%%%%%%% plasma model
flag.type_Xsec='e-H2-Phelps';

flag.input_plasma='2018Jainb_ELISE_typical';

% stoc expression
% flag.stoc_model='';
flag.stoc_model='Vahedi-simplify';
% flag.stoc_model='Cazzador-simplify';
% flag.stoc_model='2018Jainb-simplify';
% flag.stoc_model='Cazzador-fit';

flag.medium_approximation='';
% flag.medium_approximation='consider_real(sigma)_only';
% flag.medium_approximation='sigma_dc';

% flag.skin_depth='as-medium';
flag.skin_depth='as-medium';
% flag.skin_depth='collisionless';
% flag.skin_depth='collision';
% flag.skin_depth='as-medium-simplified-finite-radius';

flag.output_plasma_model=true;
% flag.output_plasma_model=false;

%%%%%%%% electric model
% flag.electric_model='';
flag.electric_model='transformer-base';
% flag.electric_model='transformer-2011Chabert';
% flag.electric_model='transformer-2015Bandyopadhyay';
% flag.electric_model='transformer-2018Jainb';
% flag.electric_model='analytical_base';

flag.input_geometry='ELISE_base';

flag.Rmetal='measured-Rmetal-woplasma';
% flag.Rmetal=calculated-Rcoil-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';

flag.output_electric_model=true;
% flag.output_electric_model=false;

%% solving
if ~flag.using_stored_data
    % calculate
    log_name=['./others/' solution_name '/' program_name '.log'];
    diary(log_name) % append to the end of the log file
    fprintf('\n-----%s %s-----\n\n',program_name,now_str)

    input=get_input_data( flag );
    
    % modify geometry
    % 应有flag.skin_depth='as-medium-simplified-finite-radius';
    % 等离子体长度取为线圈长度，密度在线圈长度内平均
    % geometry中的plasma尺寸需改变
    % 区分source中的plasma尺寸（电模型等效尺寸），和plasma中的plasma尺寸（真实尺寸）
    
    input.plasma=plasma_model(flag, input.plasma);
    
    if isfield(flag,'electric_model') && ~isempty(flag.electric_model)
        source=electric_model( flag, input );
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
input.plasma.sigma
input.plasma.eps_r

if ~flag.using_stored_data
    fprintf('\n-----END %s-----\n\n',now_str)
    diary off
end