function [ plasma ] = equivalent_EM_medium_model( flag, plasma)
% equivalent EM medium model of plasma
% Tested in test_plasma_model.m.
constants=get_constants();
if ~isfield(plasma,'nu_eff') || isempty(plasma.nu_eff)
    error('No nu_eff value. Please run ICP heating model first.')
end
%% description from the electrodynamic perspective
% assumed: not magnetized
plasma.mu_c_r=1; %等离子体复相对磁导率，目前取为真空
plasma.eps_c_r=1-plasma.wpe.^2./plasma.w_RF./...
    (plasma.w_RF-1i*plasma.nu_eff);  % 等离子体复相对介电常数，即常见的eps_p

if isfield(flag,'medium_approximation')
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
end

plasma.mu_c=plasma.mu_c_r*constants.mu0;
plasma.eps_c=constants.eps0*plasma.eps_c_r; %等离子体复介电常数
plasma.sigma_c=1i*plasma.w_RF.*plasma.eps_c; %复电导率，考虑eps0
plasma.sigma_p=1i*plasma.w_RF.*(plasma.eps_c...
    -constants.eps0); % 等离子体复电导率，即常见的sigma_p

% 分析是否可以忽略位移电流，使用复电导率
plasma.ratio_displacement_current=(plasma.w_RF.*...
    sqrt(plasma.w_RF.^2+plasma.nu_eff.^2))./plasma.wpe.^2;
% 若驱动频率较低，上述比值远小于1，则可忽略位移电流

%分析是否频率非常低到可以使用直流电导率-忽略了复电导率虚部
plasma.v2wRF=plasma.nu_eff./plasma.w_RF;
plasma.sigma_dc=constants.eps0*plasma.wpe.^2./plasma.nu_eff; %直流电导率

% 分析带电粒子是否响应电磁场
plasma.wpe2wRF=plasma.wpe./plasma.w_RF;
plasma.wpi2wRF=plasma.wpi./plasma.w_RF;
% 若驱动频率较低，上述比值远大于1，则粒子可以响应电磁场

%% equivalent EM medium
% 目前未考虑复磁导率
% lossy medium with ε'-jε''
plasma.eps_prime=real(plasma.eps_c); %即实数eps
plasma.eps_double_prime=-imag(plasma.eps_c);
plasma.tan_delta=plasma.eps_double_prime./plasma.eps_prime; %损耗正切tanδ=ε''/ε'

% medium with 3 parameters, σ, ε, μ
plasma.mu_r=real(plasma.mu_c_r);
plasma.eps_r=real(plasma.eps_c_r);
plasma.sigma=real(plasma.sigma_c);

%% wave analysis
if ~isfield(flag,'skin_depth')
    flag.skin_depth='';
end

if strcmp(flag.input_plasma,'given_directly') && ...
        isfield(plasma,'skin_depth')
    disp('[WARN] Use given δ_skin data.');
    given_skin_depth=plasma.skin_depth;
    plasma = wave_analysis( plasma, flag.skin_depth );
    plasma.skin_depth=given_skin_depth;
else
    fprintf('[INFO] Use skin depth formula: %s \n',flag.skin_depth );
    plasma = wave_analysis( plasma, flag.skin_depth );
end

end