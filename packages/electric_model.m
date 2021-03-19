function [ source ] = electric_model( flag, input )
% 螺线管线圈ICP源的电模型
constants=get_constants();

%% 几何校正
% l_equ=l_chamber; %等离子体、线圈均与腔室同长
% %         N_equ=N_coil;
% N_equ=N_coil*l_equ/l_coil; %变换线圈长度
% r_p=r_chamber; %等离子体与腔室同半径

    % 增大线圈阻抗的校正方法
    %             % 理论计算线圈阻抗
    %             fprintf('几何校正前线圈阻抗：')
    %             fprintf('%s = %.2e , ','Rcoil',Rcoil);
    %             fprintf('%s = %.2e \n','Lcoil',Lcoil);
    %             Rcoil_modified=N_equ*2*pi*r_coil*sqrt(constants.mu0*w_RF*constants.sigma_Cu/2)/2/pi/r_wire_of_coil/constants.sigma_Cu;%已考虑集肤效应，未考虑邻近效应
    %             Lcoil_modified=constants.mu0*check_Nagaoka(2*r_coil/l_equ)*pi*r_coil^2*N_equ^2/l_equ;  %一次螺线管线圈电感理论计算表达式
    %             fprintf('几何校正后，理论计算线圈阻抗：')
    %             fprintf('%s = %.2e , ','Rcoil',Rcoil_modified);
    %             fprintf('%s = %.2e \n','Lcoil',Lcoil_modified);
    %             if Rcoil<Rcoil_modified
    %                 Rcoil=Rcoil_modified*ones(size(Rcoil));
    %                 fprintf('使用几何校正后')
    %             else
    %                 fprintf('使用几何校正前')
    %             end
    %             fprintf('%s = %.2e \n','Rcoil',Rcoil);
    %             if abs(Lcoil)<abs(Lcoil_modified)
    %                 Lcoil=Lcoil_modified;
    %                 fprintf('使用几何校正后')
    %             else
    %                 fprintf('使用几何校正前')
    %             end
    %             fprintf('%s = %.2e \n','Lcoil',Lcoil);
    %             Qcoil=w_RF*Lcoil/Rcoil;

%% 电模型
fprintf('[INFO] Use electric model: %s\n',flag.electric_model);
if isnan(input.external.Icoil_rms)
    input.external.Icoil_rms=1;
    disp('[INFO] Use Icoil_rms=1A to calculate power.')
end
if strfind(flag.electric_model,'transformer')
    source=transformer_model( flag, input );
else
    switch flag.electric_model
        case 'analytical_base'
            source=analytical_EM_model( flag, input );
        case 'multi-filament'
            error('To be realized.')
        otherwise
            warning('Unexpected electric model. Please stop and check.')
            pause
    end
end

%% 功率与阻抗计算

source.PTE=source.PER./source.Rsys;   %射频功率传输效率 PTE
source.PCF=source.Rsys./source.Xsys;  %激励器射频功率耦合因数

%% 输出与可视化
if flag.output_electric_model
    output_electric_model( flag, source )
end

end