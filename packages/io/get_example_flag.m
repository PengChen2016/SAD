function [ flag ] = get_example_flag(code)
% example for flag
%% plasma model
switch code
    case 0
        flag.input_plasma='2018Jainb_ELISE_typical';
    case 2
        flag.input_plasma='2020Chen_NIS_sweep1';
    case 3
        flag.input_plasma='2021Chen_NIS_sweep_all';
    case {'2018Jainb','2018jainb','J','j'}
        flag.input_plasma='2018Jainb_ELISE_sweep_f';
        flag.type_Xsec='e-H2-Phelps';
        flag.stoc_model='2018Jainb-simplify';
        flag.medium_approximation='';
        flag.skin_depth='as-medium-simplified-finite-radius';
        flag.output_plasma_model=false;
        flag.electric_model='transformer-2018Jainb';
        flag.input_geometry='ELISE_base';
        flag.Rmetal='measured-Rmetal-woplasma';
        flag.Lcoil='calculated-Lcoil-woplasma';
        flag.output_electric_model=false;
        return
    case {'2011Chabert','2011chabert','c','C','a'}
        % TODO: 2011Chabert可能没有使用nu_eip！
        flag.type_Xsec='e-Ar-Biagi';
        flag.input_plasma='2011Chabert';
        flag.stoc_model='';
        flag.medium_approximation='';
        flag.skin_depth='as-medium-simplified'; 
        flag.output_plasma_model=false;
        flag.electric_model='analytical_base';
        flag.input_geometry='2011Chabert';
        flag.Rmetal='calculated-Rcoil-woplasma';
        flag.Lcoil='calculated-Lcoil-woplasma';
        flag.output_electric_model=false;
        return
end
flag.type_Xsec='e-H2-Phelps';
flag.stoc_model='Cazzador-fit';
flag.medium_approximation='';
flag.output_plasma_model=false;

%% electric model
flag.electric_model='transformer-base';
flag.input_geometry='ELISE_base';
flag.Rmetal='measured-Rmetal-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
flag.output_electric_model=false;
end