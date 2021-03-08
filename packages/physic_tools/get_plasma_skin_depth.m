function [ skin_depth ] = get_plasma_skin_depth( type, f, vc, wpe, r_plasma )
% 由频率计算等离子体集肤深度 calculate plasma skin depth from frequency
constants=get_constants();
% TODO：重构为矢量化；适用性检查

%% 处理输入
% 若输入不同size，则做扩展
num=length(wpe(:));
if num>1
    size_in=size(wpe);
    assert(isequal(size(vc),size_in))
    if length(f(:))<num
        assert(1==length(f(:)))
        f=f*ones(size_in);
    end
end

w_RF=2*pi*f;

%% 计算
switch type
    case 'as-medium'
        % 电磁集肤深度
        % 2014Cazzador表达式
        temp=vc.^2+w_RF.^2;
        A=constants.mu0*constants.eps0*w_RF.^2.*(1-wpe.^2./temp);
        B=constants.mu0*constants.eps0*vc.*w_RF.*wpe.^2./temp;
        skin_depth=(2./B).*sqrt((A+sqrt(A.^2+B.^2))/2);
    case 'as-medium-simplified'
        % wpe>>w,v时，不考虑有限半径的电磁集肤深度
        % 1995Vahedi/2018Jain等使用
        temp=1+vc.^2./w_RF.^2;
        skin_depth=(constants.c./wpe).*sqrt(2*temp./(1+sqrt(temp)));
    case 'collisionless'
        % wpe>>w,v时，无碰撞集肤深度
        skin_depth=constants.c./wpe;
    case 'collision'
        % TODO
    case  'as-medium-simplified-finite-radius'
        % wpe>>w,v时，考虑有限半径的电磁集肤深度。适用范围存疑
        % 1995Vaheda推导
        % 但在r_plasma<≈skin depth时，skin_depth不再能代表波穿透深度
        temp=1+vc.^2./w_RF.^2;
        a_geo=(constants.c./wpe).^2.*((3.83/r_plasma)^2-(w_RF/constants.c).^2);
        b_geo=1+temp.*a_geo;
        skin_depth=(constants.c./wpe).*sqrt(2*temp./b_geo./...
            (1+sqrt(1+vc.^2./w_RF.^2./b_geo.^2)) );
    otherwise
        error('No such type.')
end


end