function [ plasma ] = ICP_heating_model( flag, plasma)
% ICP heating model
% 主要基于2014Cazzador、1995Vahedia，借鉴2018Jain
if strcmp(flag.input_plasma,'given_directly')
    fprintf('[WARN] Use ICP heating model: Set ν directly \n');
    given_nu_m=plasma.nu_m;
    plasma=Ohmic_heating_model(plasma, 'Phelps-m');
    plasma.nu_m=given_nu_m;
    if isfield(flag,'stoc_model') && ~isempty(flag.stoc_model)
        given_nu_st=plasma.nu_st;
        plasma=stochastic_heating_model(flag.stoc_model, plasma);
        plasma.nu_st=given_nu_st;
    end
else
    plasma=Ohmic_heating_model(plasma, 'Phelps-m');
    if isfield(flag,'stoc_model') && isempty(flag.stoc_model)
        fprintf('[INFO] Use ICP heating model: Ohmic heating. \n');
    else
        fprintf('[INFO] Use ICP heating model: Ohmic+stochastic heating \n');
        plasma=stochastic_heating_model(flag.stoc_model, plasma);
    end
end

% 考虑两种加热机制后
if ~isfield(flag,'stoc_model') || isempty(flag.stoc_model)
    plasma.nu_eff=plasma.nu_m;
else
    plasma.nu_eff=plasma.nu_st+plasma.nu_m;%有效碰撞频率
end

% 分析等效碰撞频率中的主导加热机制
plasma.nu_m2nu_eff=plasma.nu_m./plasma.nu_eff;
% 若上述比值大于1，则欧姆加热占主

% 分析带电粒子是否响应电磁场
plasma.wpi2wRF=plasma.wpi./plasma.w_RF;
% 若驱动频率较低，上述比值远大于1，则粒子可以响应电磁场
end