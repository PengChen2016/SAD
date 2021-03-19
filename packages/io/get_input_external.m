function [ input ] = get_input_external( flag, geometry, w_RF )
% 返回预置 ICP源外电路数据 结构体

% input.Rcoil_th % 线圈电阻计算值
% input.Lcoil_th % 线圈自感计算值
% input.Rmetal_ex % 金属部件电阻测量值
% input.Rcoil_ex % 线圈电阻测量值
% input.Lcoil_ex % 线圈电感测量值
% input.Icoil_rms % 射频电流[A]
% input.Rmetal % 实际使用的线圈电阻
% input.Lcoil % 实际使用的线圈自感
% input.Qcoil % 品质因素

%% 线圈阻抗计算值-几何决定
constants=get_constants();% 全局常数 结构体
% TODO:考虑FEM中引线长度，计算一个理论电阻
% TODO: 验算以下公式是否在适用范围内
% 计算电阻：考虑导线表面集肤效应，未考虑邻近效应
l_wire=geometry.N_coil*2*pi*geometry.r_coil;
C_wire=2*pi*geometry.r_wire;
delta_Cu=sqrt(2./(w_RF*constants.mu0*constants.sigma_Cu));
S_current_path=delta_Cu.*C_wire;
input.Rcoil_th=l_wire./(constants.sigma_Cu.*S_current_path);

% 线圈电感
switch flag.electric_model
    case 'transformer-2018Jainb'
        input.Lcoil_th=geometry.N_coil^2*0.002*pi*(2*geometry.r_coil*100)* ...
            (log(4*2*geometry.r_coil/geometry.l_coil)-0.5)*1e-6;
    otherwise
        L_infinite_long=constants.mu0*pi*geometry.r_coil^2*geometry.N_coil^2/geometry.l_coil;
        input.Lcoil_th=L_infinite_long*get_Nagaoka(2*geometry.r_coil/geometry.l_coil);
        % 20210304 pengchen reviewed these formulas: ok,
        % same as previous version by pengzhao/chenzuo
end

%% 导入Zloss实验值
input.Icoil_rms=nan; % 计算功率绝对值时需要电流实验值
% 实验测量线圈阻抗（一般按照相减方法来测量）
switch geometry.flag
    case 'ELISE_base'
        % 2018Jain
        input.Rcoil_ex=0.5;
        input.Rmetal_ex=nan;
        input.Lcoil_ex=7.5e-6;
        % 2018Jainb在变压器模型中加入了一个表征其他金属损耗的电阻
    case 'HUST_large_driver_base'
        % TODO: 待校正
        input.Rcoil_ex=nan;
        input.Rmetal_ex=1.45;
        input.Lcoil_ex=16.04e-6;
    case 'BATMAN_base'
        warning('no data')
    case 'HUST_small_driver_base'
        input.Rcoil_ex=nan;
        input.Rmetal_ex=0.8;
        input.Lcoil_ex=18.7e-6;
    case 'HUST_small_driver_ZL'
        input.Rcoil_ex=70.22e-3;
        input.Rmetal_ex=75.22e-3;
        input.Lcoil_ex=3.67e-6;
        % ZL测得带线圈、FS、背板小激励器电感
        % input.Ls_ex=3.21e-6;
    case 'CHARLIE_base'
        % 背板等损耗可忽略
        input.Rcoil_ex=[[0.0970835, 0.0955728, 0.0889425, 0.0892782, 0.0940621, 0.08953]',...
            [0.245888, 0.209715, 0.180256, 0.160029, 0.183613, 0.213659]'];
        input.Rmetal_ex=input.Rcoil_ex;
        input.Lcoil_ex=2.2e-6;
        input.Icoil_rms=[[51.9436, 35.4548, 28.1894, 28.3519, 32.6933, 44.2288]',...
            [17.5453, 14.0711, 12.2062, 10.0207, 9.26352, 9.66274]'];
    case 'small_source1_LZS'
        % 线圈以外无其他金属
        input.Rcoil_ex=1;
        input.Rmetal_ex= input.Rcoil_ex;
        input.Lcoil_ex=2.2e-6; %原文未给出，不关心
    case 'NIO1_base'
        warning('no data')
end

%% 用户选择使用
fprintf('[INFO] Choose %s as Rloss with plasma\n',flag.Rmetal)
switch flag.Rmetal
    case 'measured-Rmetal-woplasma'
        if isnan(input.Rmetal_ex)
            disp('No measured-Rmetal-woplasma. ')
            if isnan(input.Rcoil_ex)
                error('no data')
            else
                disp('Use measured-Rcoil-woplasma instead.')
                input.Rmetal=input.Rcoil_ex;
            end
        else
            input.Rmetal=input.Rmetal_ex;
        end
    case 'calculated-Rcoil-woplasma'
        input.Rmetal=input.Rcoil_th;
end

fprintf('[INFO] Use %s as Lcoil with plasma\n',flag.Lcoil)
switch flag.Lcoil
    case 'measured-Lcoil-woplasma'
        input.Lcoil=input.Lcoil_ex;
    case 'calculated-Lcoil-woplasma'
        input.Lcoil=input.Lcoil_th;
end

input.Qcoil=w_RF.*input.Lcoil./input.Rmetal;

end