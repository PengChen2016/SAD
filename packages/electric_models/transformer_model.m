function [ source ] = transformer_model( flag, input )
% transformer model
% Ref: 1992Piejak, 2011Chabert, 2015Bandyopadhyay, 2018Jain
constants=get_constants();

%% in data struct disassembly
geometry=input.geometry;
external=input.external;
plasma=input.plasma;
w_RF=plasma.w_RF;

%% circuit element
switch flag.electric_model
    case 'transformer-base'
        % HUST NIS Group
        % 二次侧等离子体电阻
        % 根据长直导体（sigma_p）圆柱周向涡流模型，即变压器模型常用的Rp推导方法，Rp应带2
        % sigma_dc是用sigma_p推导后，为了简化表达式而使用。
        Rp=2*pi*geometry.r_plasma_eff./...
            (plasma.skin_depth.*plasma.sigma_dc*geometry.l_plasma);
        % 20210511 改用2而非1
        % 2011Chabert中区分了高气压极限和低气压极限（ν与ω比值）下的不同Rp简化表达式，
        % 在低气压极限下时使用无碰撞趋肤深度可以得到一个不带2的简化表达式，
        % 但在高气压极限下时使用碰撞趋肤深度会得到一个带2的简化表达式。
        % 随ne增大，可能Rp减小，而PER增大。此处使用1/2会对PER带来显著影响
        
        % 二次侧等离子体中电子惯性形成的电感。
        Lp=Rp./plasma.nu_eff;
        
        % 二次侧等离子体感应电感
        % 无限长电感公式
        Lmp_infinite_long=constants.mu0*pi*geometry.r_plasma_eff^2/geometry.l_plasma;
        % 长冈系数做有限长度校正
        Lmp=get_Nagaoka(2*geometry.r_plasma_eff/geometry.l_plasma)*Lmp_infinite_long;
        
        % 线圈与等离子体互感
        % 耦合系数，1992Piejak
        kp=(geometry.r_plasma_eff/geometry.r_coil)^2;
        % 耦合系数对理想变互感公式做校正
        M=kp*sqrt(external.Lcoil*(Lmp+Lp));
    case 'transformer-2011Chabert'
        % according to the analytical EM model
        error('Wait')
        
        
        
    case 'transformer-2018Jainb'
        % 长直圆柱周向涡流模型-等离子体导体视角
        Rp=2*pi*geometry.r_plasma_eff./...
            (plasma.skin_depth.*plasma.sigma_dc*geometry.l_plasma);
        Lp=Rp./plasma.nu_eff;
        
        % 考虑几何效应的自感 % formula 75 in 2018Jainb
        Lmp=0.002*pi*(2*geometry.r_plasma_eff*100)*...
            (log(4*2*geometry.r_plasma_eff/geometry.l_plasma)-0.5)*1e-6;
        % 该表达式不适用于D/L<0.4的长线圈
        if Lmp<0
            warning('Lmp<0. Use Lmp=Nagaoka*Lmp_infinite_long.')
            % 无限长电感公式
            Lmp_infinite_long=constants.mu0*pi*geometry.r_plasma_eff^2/geometry.l_plasma;
            % 长冈系数做有限长度校正
            Lmp=get_Nagaoka(2*geometry.r_plasma_eff/geometry.l_plasma)*Lmp_infinite_long;
        end
        % 考虑有限长度的线圈互感 % formula 76 in 2018Jainb
        M=0.0095*geometry.N_coil*1e-6*...
            (2*geometry.r_plasma_eff*100)^2/...
            sqrt((2*geometry.r_coil*100)^2+(geometry.l_coil*100)^2);
        kp=M./sqrt(external.Lcoil*(Lmp+Lp));
    case 'transformer-2015Bandyopadhyay'
        % 2015Bandyopadhyay % ignore the difference of skin_depth
        % Re(σ_p) consider_real(sigma)_only, wrong
        Rp=2*pi*geometry.r_plasma_eff/...
            (plasma.skin_depth.*plasma.sigma*geometry.l_plasma);
        Lp=Rp./plasma.nu_eff;
        
        % 二次侧等离子体感应电感，未考虑有限长度
        % 耦合系数，1992Piejak
        kp=(geometry.r_plasma_eff/geometry.r_coil)^2;
        % 由一次侧线圈与耦合系数计算 1992Piejak
        Lmp=external.Lcoil*kp/geometry.N_coil^2;
        M=kp*sqrt(external.Lcoil*(Lmp+Lp));
    otherwise
        error('No such type.')
end

%% circuit
% 等离子体换算到一次侧, 2011Chabert
source.PER=w_RF.^2.*M.^2.*Rp./...
    (Rp.^2+w_RF.^2.*(Lmp+Lp).^2);  % 一次侧等离子体等效电阻
source.Lplasma=-w_RF.^2.*M.^2.*(Lmp+Lp)./...
    (Rp.^2+w_RF.^2.*(Lmp+Lp).^2); % 一次侧等离子体等效电感
% 注意：与解析电磁模型中Lplasma定义不同
idx=find(source.Lplasma>0,1);
if ~isempty(idx)
    disp('Lplasma of the following idx >0:')
    disp(idx')
    error('Lplasma of ICP source should < 0')
end
source.Iplasma_rms=external.Icoil_rms.*...
    (-1i*w_RF.*M./(1i*w_RF.*(Lmp+Lp)+Rp)); % phase of Icoil is assumed as 0.
source.Rmetal=external.Rmetal;
source.Xplasma=source.Lplasma.*w_RF;

%系统等效阻抗
source.Rsys=source.Rmetal+source.PER; % 系统电阻
source.Lsys=external.Lcoil+source.Lplasma;   %系统电感
source.Xsys=source.Lsys.*w_RF;   %系统电抗

%% 功率
source.Pplasma=source.PER.*external.Icoil_rms.^2; % 等离子体吸收功率
source.Psys=source.Rsys.*external.Icoil_rms.^2;

%% out data struct assembly
source.transformer.Rp=Rp;
source.transformer.Lp=Lp;
source.transformer.Lmp=Lmp;
source.transformer.kp=kp;
source.transformer.M=M;

end