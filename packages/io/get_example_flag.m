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