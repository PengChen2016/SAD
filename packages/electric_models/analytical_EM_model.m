function [ source ] = analytical_EM_model( flag, input )
% transformer model
% Ref: 2011Chabert
constants=get_constants();

%% in data operation
geometry=input.geometry;
external=input.external;
plasma=input.plasma;
w_RF=plasma.w_RF;

% one part (l_part=l_coil=l_plasma) of the infinite long model
if geometry.l_plasma~=geometry.l_coil
    fprintf('[WARN] l_coil≠l_plasma. Use the shorter one: ')
    if geometry.l_plasma>geometry.l_coil
        l_part=geometry.l_coil;
        fprintf('l_part=l_coil. \n')
    else
        l_part=geometry.l_plasma;
        fprintf('l_part=l_plasma. \n')
    end
else
    l_part=geometry.l_coil;
end
% Then, 下文只能出现l_part，不能出现l_plasma或l_coil

%% EMF
Hz0=geometry.N_coil*external.Icoil_rms/l_part; % H_rms at r_coil
% check the condition for uniform H
eps_tube_r=9; % eps_r of dielectric tube. 2 for glass, 9 for ceramic.
k_wave_tube=w_RF*sqrt(constants.mu0*constants.eps0*eps_tube_r);
if ~isempty(find(k_wave_tube*geometry.r_plasma_eff>0.2,1))
    error('介质管中磁场不是常数，该模型不适用')
end

% r<=r_plasma
J0_kr0=besselj(0,plasma.k_wave*geometry.r_plasma_eff);
Hz_plasma=@(r) Hz0.*besselj(0,plasma.k_wave*r)./J0_kr0;
Etheta_plasma=@(r) -1i*plasma.k_wave.*Hz0./(w_RF.*plasma.eps_c).*...
    besselj(1,plasma.k_wave*r)./J0_kr0;
S_plasma=@(r) -Hz_plasma(r).*Etheta_plasma(r);
source.PQplasma=S_plasma(geometry.r_plasma_eff)*2*pi*geometry.r_plasma_eff*l_part;
source.Pplasma=real(source.PQplasma);
source.Xplasma=imag(source.PQplasma)./external.Icoil_rms.^2;
source.Lplasma=source.Xplasma./w_RF;
% 注意：与变压器模型中Lplasma定义不同

% impedance
source.PER=source.Pplasma./external.Icoil_rms.^2; % PER=Rind
source.Rmetal=external.Rmetal;
source.Rsys=source.Rmetal+source.PER; 
source.Psys=source.Rsys.*external.Icoil_rms.^2;

% r<=r_coil
% EMF is confined in this region, so this region represents the whole system
Hz_no_plasma=@(r) Hz0; % H=Hz0 where no plasma
Etheta0=Etheta_plasma(geometry.r_plasma_eff)...
    *geometry.r_plasma_eff/geometry.r_coil...
    -1i*w_RF*constants.mu0.*Hz0...
    .*(geometry.r_coil^2-geometry.r_plasma_eff^2)/(2*geometry.r_coil);...
    % E_rms at r_coil
S0=-Etheta0.*Hz0;
source.PQ0=S0*2*pi*geometry.r_coil*l_part;
source.Xsys=imag(source.PQ0)./external.Icoil_rms.^2; %Xind
source.Lsys=source.Xsys./w_RF; % Lind
source.Lm=source.Lsys-source.PER./plasma.nu_eff;

% transformer
Iplasma=l_part*Hz0.*(1./J0_kr0-1); % Attention: I and l are different characters.
source.Iplasma_rms=abs(Iplasma);
source.Rp=source.Pplasma./source.Iplasma_rms.^2; 
source.Lp=source.Rp./plasma.nu_eff;
% TODO: transformer model 2011chabert
source.transformer.Lmp=constants.mu0*pi*geometry.r_plasma_eff^2/l_part;
% The low-pressure regime
skin_depth1=get_plasma_skin_depth('collisionless',...
    plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);
source.transformer.Rp1=pi*geometry.r_plasma_eff./(plasma.sigma_dc.*skin_depth1*l_part); 
% The high-pressure regime
skin_depth2=get_plasma_skin_depth('collision',...
    plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);
source.transformer.Rp2=2*pi*geometry.r_plasma_eff./(plasma.sigma_dc.*skin_depth2*l_part); 

%% out data struct assembly
source.emf.Hz_no_plasma=Hz_no_plasma;
source.emf.Hz_plasma=Hz_plasma;
source.emf.Etheta_plasma=Etheta_plasma;
source.emf.S_plasma=S_plasma;

% 逻辑表达式实现匿名分段函数
source.emf.Bzm_r=@(r) constants.mu0*abs((0<=r & r<=geometry.r_plasma_eff)*Hz_plasma(r)...
    +(geometry.r_plasma_eff<r & r<=geometry.r_coil)*Hz_no_plasma(r));

end