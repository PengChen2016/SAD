% CHARLIE_sweep case for paper
% 20210415 created by pengchen2016, Matlab R2017a

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

% %------------------------------ 1. raza, transformer-base
% program_name='paper_CHARLIE_raza_t210517';
% %%%%%%%% plasma model
% flag.type_Xsec='e-H2-Phelps';
% flag.input_plasma='CHARLIE_raza_sweep';
% flag.stoc_model='Cazzador-fit';
% flag.medium_approximation='';
% flag.skin_depth='as-medium';
% flag.output_plasma_model=true;
% %%%%%%%% electric model
% flag.electric_model='transformer-base';
% flag.input_geometry='CHARLIE_base';
% flag.Rmetal='measured-Rmetal-woplasma';
% flag.Lcoil='measured-Lcoil-woplasma';
% flag.output_electric_model=true;
% %------------------------------ 1. flag end

% %------------------------------ 2. raza, analytical-base
% program_name='paper_CHARLIE_raza_a210517';
% %%%%%%%% plasma model
% flag.type_Xsec='e-H2-Phelps';
% flag.input_plasma='CHARLIE_raza_sweep';
% flag.stoc_model='Cazzador-fit';
% flag.medium_approximation='';
% flag.skin_depth='as-medium';
% flag.output_plasma_model=true;
% %%%%%%%% electric model
% flag.electric_model='analytical_base';
% flag.input_geometry='CHARLIE_base';
% flag.Rmetal='measured-Rmetal-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';
% flag.output_electric_model=true;
% % flag.magnetized='TODO';
% %------------------------------ 2. flag end

% %------------------------------ 3. raza, transformer-2018Jainb
% program_name='paper_CHARLIE_raza_tj210517';
% %%%%%%%% plasma model
% flag.type_Xsec='e-H2-Phelps';
% flag.input_plasma='CHARLIE_raza_sweep';
% flag.stoc_model='2018Jainb-simplify';
% flag.medium_approximation='';
% flag.skin_depth='as-medium-simplified-finite-radius';
% flag.output_plasma_model=true;
% %%%%%%%% electric model
% flag.electric_model='transformer-2018Jainb';
% flag.input_geometry='CHARLIE_base';
% flag.Rmetal='measured-Rmetal-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';
% flag.output_electric_model=true;
% %------------------------------ 3. flag end

% program_name='CHARLIE_raza_sweep_tj210511';
% %%%%%%%% plasma model
% flag.type_Xsec='e-H2-Phelps';
%
% flag.input_plasma='CHARLIE_raza_sweep';
% % flag.input_plasma='2019Raunera_CHARLIE_sweep';
% % flag.input_plasma='CHARLIE_razcoil_sweep';
%
% % stoc expression
% % flag.stoc_model='';
% % flag.stoc_model='Vahedi-simplify';
% % flag.stoc_model='Cazzador-simplify';
% % flag.stoc_model='2018Jainb-simplify';
% flag.stoc_model='Cazzador-fit';
%
% flag.medium_approximation='';
% % flag.medium_approximation='consider_real(sigma)_only';
% % flag.medium_approximation='sigma_dc';
%
% flag.skin_depth='as-medium';
% % flag.skin_depth='as-medium-simplified';
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
% flag.input_geometry='CHARLIE_base';
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
plot_1Y(input.plasma.ne, 'n_e [m^{-3}]');
axis([0.3,10,5e16,5e17])
plot_1Y(input.plasma.Te, 'T_e [eV]');
axis([0.3,10,0,6])
plot_1Y(input.external.Rmetal, 'R_{loss} [\Omega]');
axis([0.3,10,0,0.3])
plot_1Y(input.external.Icoil_rms, 'I_{rms} [A]');
axis([0.3,10,0,60])

% check model
% nu
plot_2Yaxis(input.plasma.nu_m, '\nu_m [Hz]', input.plasma.nu_st, '\nu_{st} [Hz]');
yyaxis left
axis([0.3,10,2e6,2e8])
line([0.3,10],[1e6*2*pi,1e6*2*pi],'linestyle','--')
line([0.3,10],[4e6*2*pi,4e6*2*pi],'linestyle','--')
yyaxis right
axis([0.3,10,2e6,2e8])
plot_1Y(input.plasma.nu_eff, '\nu_{eff} [Hz]');
% skin_depth
plot_2Yaxis(input.plasma.skin_depth, '\delta_{skin} [m]',input.plasma.wavelength, '\lambda [m]');
yyaxis left
line([0.3,8],[input.plasma.r,input.plasma.r],'linestyle','--')
plot_2Yaxis(input.plasma.sigma_dc, '\sigma_{dc} [S/m]', input.plasma.skin_depth, '\delta_{skin} [m]');
% power and impedance
if strfind(flag.electric_model,'analytical')
    plot_2Yaxis(source.Rp, 'R_p [\Omega]',source.PER, 'R_{plasma} [\Omega]');
else
    plot_2Yaxis(source.transformer.Rp, 'R_p [\Omega]',source.PER, 'R_{plasma} [\Omega]');
end
plot_2Yaxis(source.Pplasma, 'P_{plasma} [W]',source.Psys, 'P_{sys} [W]');
plot_2Yaxis(source.Rsys, 'R_{sys} [\Omega]',source.Xsys, 'X_{sys} [\Omega]');

plot_1Ylog(source.PER, 'R_{plasma} [\Omega]');
axis([0.3,10,0.1,10])
plot_1Y(source.PTE, '\eta');
axis([0.3,10,0.5,1])
% meidum
plot_2Yaxis(input.plasma.sigma, '\sigma [S/m]', -input.plasma.w_RF.*input.plasma.eps_prime, '-\omega*\epsilon');
yyaxis left
axis([0.3,10,0,100])
yyaxis right
axis([0.3,10,0,100])
if strfind(flag.electric_model,'analytical')
    plot_1Y(source.Lp, 'L');
else
    plot_1Y(source.transformer.Lp, 'L');
end
line([0.3,10],[source.transformer.Lmp,source.transformer.Lmp],'linestyle','--')
legend('L_p, 1MHz','L_p, 4MHz','L_{mp}');
% other
plot_2Yaxis(input.plasma.wpe2wRF, '\omega_{pe}/\omega_{RF}',input.plasma.wpi2wRF, '\omega_{pi}/\omega_{RF}');

if ~flag.using_stored_data
    fprintf('\n-----END %s-----\n\n',now_str)
    diary off
end