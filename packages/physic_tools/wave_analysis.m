function [ plasma ] = wave_analysis( plasma )
% 电磁波分析
constants=get_constants();
%% 基础
plasma.k_wave=plasma.w_RF.*sqrt(plasma.mu_c.*plasma.eps_c); %wave number
plasma.gamma_wave=1i*plasma.k_wave; %传播常数
plasma.alpha_wave=-imag(plasma.k_wave); %衰减常数
plasma.beta_wave=real(plasma.k_wave); %相位常数
        
%% 等效集肤深度
plasma.skin_depth=1./plasma.alpha_wave; %一般介质中经典集肤深度
% plasma.skin_depth=get_plasma_skin_depth('as-medium-simplified',...
%     plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);

if ~isempty(find(plasma.r<3*plasma.skin_depth, 1))
    disp('[WARN] 存在δ>≈R ，电磁波穿透等离子体，电场基本均匀,集肤深度概念不适用');
end

%% 波长
plasma.wavelength=2*pi./plasma.beta_wave; %一般介质中波长

if ~isempty(find(plasma.wavelength<3*plasma.r, 1))
    disp('[WARN] 存在λ＜≈R ，需要考虑相位变化');
end

end