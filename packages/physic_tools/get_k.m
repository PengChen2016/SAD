function [ k ] = get_k( Te, Xsec )
% 无电磁场时，Maxwellian分布电子撞击重粒子反应的速率系数
constants=get_constants();%
E=Xsec(:,1);
a=Xsec(:,2);
k=zeros(size(Te));
for i=1:length(Te(:))
    %截面对Maxwellian分布积分
    k(i)=trapz(E,a.*E.*exp(-E/Te(i)))/Te(i)^1.5*sqrt(8*constants.e/pi/constants.me); 
end
end

% 参见e:\GitRepos\Nix_M\packages\plasma_physic_tools\get_k_e_heavy.m