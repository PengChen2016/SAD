function [ skin_depth ] = get_skin_depth( type, f, sigma, eps_r, mu_r, r_geo )
% 计算集肤深度 calculate skin depth
% TODO 电磁集肤深度计算
sigma,f,mu_r

constants=get_constants();

%% 不同的计算公式
switch type
    case 'conductor'
        skin_depth=sqrt(1/(pi*f*mu_r*mu0*sigma));
    case 'medium'
        % 垂直入射半无限大等离子体
        
        alpha_wave=-imag(k_p_waplasma.ve); %衰减常数
        %             skin_depth_eff(X1i,X2i)=1/alpha_wave; %一般介质中经典集肤深度
        

end



%% 公式的适用性
if r_plasma<3*skin_depth_eff(X1i,X2i)
    warning('δ>≈R ，电磁波穿透等离子体，电场基本均匀,集肤深度概念与随机加热模型不适用')
    pause
end

end