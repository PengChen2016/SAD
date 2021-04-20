% CHARLIE_raza_sweep case for paper
% 20210415 created by pengchen2016, Matlab R2017a

%% initialization
close all
clear
addpath(genpath('./packages'))
now_str=datestr(now,'yyyymmdd_HHMMSS');
disp(now_str)
tic
%% flag
flag.using_stored_data=false;
% flag.using_stored_data=true;
save_mat_name='./others/CHARLIE_raza_sweep210420.mat';
%%%%%%%% plasma model
flag.input_plasma='CHARLIE_raza_sweep';

% stoc expression
% flag.stoc_model='';
% flag.stoc_model='Vahedi-simplify';
% flag.stoc_model='Cazzador-simplify';
% flag.stoc_model='2018Jainb-simplify';
flag.stoc_model='Cazzador-fit';

flag.medium_approximation='';
% flag.medium_approximation='consider_real(sigma)_only';
% flag.medium_approximation='sigma_dc';

flag.skin_depth='';
% flag.skin_depth='as-medium';
% flag.skin_depth='as-medium-simplified';
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

flag.input_geometry='CHARLIE_base';

flag.Rmetal='measured-Rmetal-woplasma';
% flag.Rmetal=calculated-Rcoil-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';

flag.output_electric_model=true;
% flag.output_electric_model=false;

%% solving
log_name=['.\others\Run' now_str '.log'];
diary(log_name)
if ~flag.using_stored_data
    % calculate
    input=get_input_data( flag );
    input.plasma=plasma_model(flag, input.plasma);
    
    if isfield(flag,'electric_model') && ~isempty(flag.electric_model)
        source=electric_model( flag, input );
    end
    
    save(save_mat_name)
else
    % load data
    load(save_mat_name)
    warning(['Using data stored in ' save_mat_name])
    if flag.output_plasma_model
        output_plasma_model(flag,input.plasma)
    end
    if flag.output_electric_model
        output_electric_model( flag, source )
    end
end

%% post-processing
% check
% 




diary off

%% aid function
% plot Y for X(size specified in this file: p×f)
function plot_1Y()


end