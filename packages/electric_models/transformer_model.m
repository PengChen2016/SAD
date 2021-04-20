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
        % 长直圆柱周向涡流模型-等离子体电介质视角，2011Chabert
        Rp=pi*geometry.r_plasma_eff./...
            (plasma.skin_depth.*plasma.sigma_dc*geometry.l_plasma);
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
        % 长直圆柱周向涡流模型-等离子体电介质视角
        error('Wait')
        
        
        
    case 'transformer-2018Jainb'
        % 长直圆柱周向涡流模型-等离子体导体视角
        Rp=2*pi*geometry.r_plasma_eff/...
            (plasma.skin_depth.*plasma.sigma_dc*geometry.l_plasma);
        Lp=Rp./plasma.nu_eff;
        
        % 考虑几何效应的自感
        Lmp=0.002*pi*(2*geometry.r_plasma_eff*100)*...
            (log(4*2*geometry.r_plasma_eff/geometry.l_plasma)-0.5)*1e-6;
        %该表达式结果可能为负值，则bug; 即使取绝对值，也超出适用范围
        if Lmp<0
            warning('Lmp<0. Use Lmp=-Lmp.')
            Lmp=-Lmp;
        end
        % 考虑有限长度的线圈互感
        M=0.0095*geometery.N_coil*1e-6*...
            (2*geometry.r_plasma_eff*100)^2/...
            sqrt((2*geometery.r_coil*100)^2+(geometry.l_coil*100)^2);
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
idx=find(source.Lplasma>0,1);
if ~isempty(idx)
    disp('Lplasma of the following idx >0:')
    disp(idx')
    error('Lplasma of ICP source should < 0')
end
source.Rmetal=external.Rmetal;
source.Xplasma=source.Lplasma.*w_RF;

%系统等效阻抗
source.w_RF=w_RF;
source.Rsys=external.Rmetal+source.PER; % 系统电阻
source.Lsys=external.Lcoil+source.Lplasma;   %系统电感
source.Xsys=source.Lsys.*w_RF;   %系统电抗

%% 功率
source.P_abs=source.PER.*external.Icoil_rms.^2; % 等离子体吸收功率

%% out data struct assembly
source.transform.Rp=Rp;
source.transform.Lp=Lp;
source.transform.Lmp=Lmp;
source.transform.kp=kp;
source.transform.M=M;

end