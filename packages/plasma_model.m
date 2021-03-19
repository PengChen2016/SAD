function [ plasma ] = plasma_model( flag, plasma)
% ICP的等效电磁媒质模型 equivalent_medium_model_of_plasma
constants=get_constants();

%% ICP heating model
% 主要基于2014Cazzador、1995Vahedia，借鉴2018Jain
if strcmp(flag.input_plasma,'given_directly')
    fprintf('[WARN] Use ICP heating model: Set ν directly \n');
    given_nu_m=plasma.nu_m;
    plasma=Ohmic_heating_model(plasma);
    plasma.nu_m=given_nu_m;
    if ~isempty(flag.stoc_model)
        given_nu_st=plasma.nu_st;
        plasma=stochastic_heating_model(flag.stoc_model, plasma);
        plasma.nu_st=given_nu_st;
    end
else
    plasma=Ohmic_heating_model(plasma);
    if isempty(flag.stoc_model)
        fprintf('[INFO] Use ICP heating model: Ohmic heating. \n');
    else
        fprintf('[INFO] Use ICP heating model: Ohmic+stochasti heating \n');
        plasma=stochastic_heating_model(flag.stoc_model, plasma);
    end
end

% 考虑两种加热机制后
if isempty(flag.stoc_model)
    plasma.nu_eff=plasma.nu_m;
else
    plasma.nu_eff=plasma.nu_st+plasma.nu_m;%有效碰撞频率
end

% 分析等效碰撞频率中的主导加热机制
plasma.nu_m2nu_eff=plasma.nu_m./plasma.nu_eff;
% 若上述比值大于1，则欧姆加热占主

% 分析带电粒子是否响应电磁场
plasma.wpe2wRF=plasma.wpe./plasma.w_RF;
plasma.wpi2wRF=plasma.wpi./plasma.w_RF;
% 若驱动频率较低，上述比值远大于1，则粒子可以响应电磁场

% 分析是否可以忽略heating model中位移电流，使用复电导率
plasma.ratio_displacement_current=(plasma.w_RF.*...
    sqrt(plasma.w_RF.^2+plasma.nu_eff.^2))./plasma.wpe.^2;
% 若驱动频率较低，上述比值远小于1，则可忽略位移电流

%% 等离子体等效电磁媒质模型
% 等离子体 电磁模型
plasma.mu_c_r=1; %等离子体复相对磁导率，目前取为真空 
plasma.eps_c_r=1-plasma.wpe.^2./plasma.w_RF./...
    (plasma.w_RF-1i*plasma.nu_eff);  % 等离子体复相对介电常数，即常见的eps_p

%分析是否频率非常低到可以使用直流电导率-忽略了复电导率虚部
plasma.v2wRF=plasma.nu_eff./plasma.w_RF;
plasma.sigma_dc=constants.eps0*plasma.wpe.^2./plasma.nu_eff; %直流电导率

switch flag.medium_approximation
    case ''
        % 有损介质，do nothing
    case 'consider_real(sigma)_only'
        % 复电导率忽略虚部，即Re(eps_c_r)=1,Im(eps_c_r)不变
        disp('[WARN] consider_real(sigma)_only before solving EM field.')
        plasma.eps_c_r=1+1i*imag(plasma.eps_c_r);
    case 'sigma_dc'
        % 复电导率近似为直流电导率，即Re(eps_c_r)=1,Im(eps_c_r)由sigma_dc计算
        disp('[WARN] Use sigma_dc approximation.')
        plasma.eps_c_r=1+plasma.sigma_dc./(1i*plasma.w_RF*constants.eps0);
    otherwise
        error('No such type.')
end

plasma.mu_c=plasma.mu_c_r*constants.mu0;
plasma.eps_c=constants.eps0*plasma.eps_c_r; %等离子体复介电常数
plasma.sigma_c=1i*plasma.w_RF.*plasma.eps_c; %复电导率，考虑eps0
plasma.sigma_p=1i*plasma.w_RF.*(plasma.eps_c...
    -constants.eps0); % 等离子体复电导率，即常见的sigma_p

% 电磁 有损介质模型ε'-jε''
plasma.eps_prime=real(plasma.eps_c); %即实数eps
plasma.eps_double_prime=-imag(plasma.eps_c);
plasma.tan_delta=plasma.eps_double_prime./plasma.eps_prime; %损耗正切tanδ=ε''/ε'

% 电磁 三参数模型
plasma.mu_r=real(plasma.mu_c_r);
plasma.eps_r=real(plasma.eps_c_r);
plasma.sigma=real(plasma.sigma_c); % 目前未考虑复磁导率

%% 电磁波分析
if strcmp(flag.input_plasma,'given_directly') && ~strcmp(flag.electric_model,'analytical_base')
    warning('使用给定的δ_skin数据')
    assert(isfield(plasma,'skin_depth'))
    given_skin_depth=plasma.skin_depth;
    plasma = wave_analysis( plasma );
    plasma.skin_depth=given_skin_depth;
else
    plasma = wave_analysis( plasma );
end

idx=find(plasma.r<3*plasma.skin_depth);
if ~isempty(idx)
    if strfind(flag.electric_model,'transformer')...
            || isempty(flag.stoc_model)
        disp('[WARN] δ>≈R的元素的索引')
        disp(idx')
        warning('随机加热模型、变压器模型长直圆柱涡流问题电阻计算 有bug')
%         pause
    end
end

idx=find(plasma.wavelength<3*plasma.r);
if ~isempty(idx)
    if strfind(flag.electric_model,'transformer')
        disp('[WARN] λ＜≈R的元素的索引')
        disp(idx')
        warning('集中参数电路模型有bug')
%         pause
    end
end

%% 输出与可视化
if flag.output_plasma_model
    output_plasma_model(flag, plasma)
end

end