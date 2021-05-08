function [ skin_depth ] = get_plasma_skin_depth( type, f, nu_c, wpe, r_plasma )
% calculate skin depth of (cylinder) plasma from frequency
% 
% type: type of different formula
% f: RF frequency
% nu_c: collision frequency
% wpe: electron plasma frequency
% r_plasma: radius of cylinder plasma
% TODO: 在考虑几何效应时，分情况选用公式
constants=get_constants();
w_RF=2*pi*f;

%% 计算
switch type
    case 'as-medium'
        % 电磁集肤深度
        % 2014Cazzador表达式
        temp=nu_c.^2+w_RF.^2;
        A=constants.mu0*constants.eps0*w_RF.^2.*(1-wpe.^2./temp);
        B=constants.mu0*constants.eps0*nu_c.*w_RF.*wpe.^2./temp;
        skin_depth=(2./B).*sqrt((A+sqrt(A.^2+B.^2))/2);
    case 'as-medium-simplified'
        % wpe>>w,v时，不考虑有限半径的电磁集肤深度
        % 1995Vahedi/2018Jain等使用
        temp=1+nu_c.^2./w_RF.^2;
        skin_depth=(constants.c./wpe).*sqrt(2*temp./(1+sqrt(temp)));
        flag_suitable=isempty(find(wpe<10*w_RF, 1)) ...
            && isempty(find(wpe<10*nu_c, 1))...
            && isempty(find(skin_depth>r_plasma, 1));
    case 'collisionless'
        % collisionless skin depth δ, applicable to ω_RF>>ν_c
        skin_depth=constants.c./wpe;
    case 'collision'
        % collisionless skin depth δ, applicable to ω_RF<<ν_c
        skin_depth=(constants.c./wpe).*(2*nu_c./w_RF).^0.5;
    case  'as-medium-simplified-finite-radius'
        % wpe>>w,v时，考虑有限半径的电磁集肤深度。适用范围存疑
        % 1995Vaheda推导
        % 但在r_plasma<≈skin depth时，skin_depth不再能代表波穿透深度
        temp=1+nu_c.^2./w_RF.^2;
        a_geo=(constants.c./wpe).^2.*((3.83/r_plasma)^2-(w_RF/constants.c).^2);
        b_geo=1+temp.*a_geo;
        skin_depth=(constants.c./wpe).*sqrt(2*temp./b_geo./...
            (1+sqrt(1+nu_c.^2./w_RF.^2./b_geo.^2)) );
    otherwise
        error('No such type.')
end

end