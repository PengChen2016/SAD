function [ flag ] = get_flag( )
% example for flag
flag.input_plasma='2018Jainb_ELISE_typical';
% flag.input_plasma='2020Chen_NIS_sweep1';
flag.input_geometry='ELISE_base';
flag.Rmetal='measured-Rmetal-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
flag.electric_model='transformer_base';
flag.stoc_model='Cazzador-fit';
flag.medium_approximation='';
end