% -*- coding: utf-8 -*-
% ----------------------------------------------
%{
 * brief #Abstract
 * Created 13:10:30 2022/11/07
 * author
   Zengshan Li, HUST: origin
        新的可用的整体模型代码.rar 2020/09/03
        global_model_LZS_20200801\Untitled2.m ver2019
   Peng Chen, HUST: modified 整合进SAD
 * note #Detail
    粒子种类：H H2 H+ H2+ H3+ H-
 * #TODO
%}
% ----------------------------------------------
function [Xt]=plasma_evolution(flag,input)

if isfield(flag,'global_model') && ~isempty(flag.global_model)
    constants=get_constants();
    
    % 时间积分
    X0=[input.X.nHi input.X.nH2i input.X.nH3i input.X.nHNi...
        1.5*input.X.Te*constants.e*input.X.ne...
        input.X.nH input.X.nHH]; % 初值
    
    [Xt.t, X] = ode15s(@(t,X) one_timestep(t,X,flag,input), input.X.tspan, X0);
    % 求解 X'=fun(t,X) (上式中为匿名函数) 在tspan上的积分
    % ode15s 是基于 1 到 5 阶数值微分公式 (NDF) 的可变步长、可变阶次 (VSVO) 解算器。
    % ode15s自动调整积分的时间步长
    
    % 由于 myOde2 实际不含时，求得微分量乘以步长（可能其他处理），
    % 即作为下次myOde2的初始值
    
    Xt.nHi=X(:,1);
    Xt.nH2i=X(:,2);
    Xt.nH3i=X(:,3);
    Xt.nHNi=X(:,4);
    Xt.ne=Xt.nHi+Xt.nH2i+Xt.nH3i-Xt.nHNi;
    Xt.Te=X(:,5)./1.5./constants.e./Xt.ne;
    Xt.nH=X(:,6);
    Xt.nHH=X(:,7);
end

end

function [Res]=one_timestep(t,X,flag,input)

if isfield(flag,'electric_model') && ~isempty(flag.electric_model)
    %%%%%%% 电模型
    % 整体模型提供新等离子体参数
    input.plasma.ne=X.ne;
    input.plasma.Te=X.Te;
    input.plasma.ng=input.plasma.nHH;
    input.plasma=get_characterstic_parameters( flag, input.plasma );
    
    input.plasma=plasma_model(flag, input.plasma);
    % TODO：检查更新后的input，是否等价于完整的新input
    source=electric_model( flag, input );
    % 得到的source.Pplasma是根据所输入Icoil而非Ps计算的，因此只用source.PTE
else
    source.PTE=1;
end

%%%%%% 整体模型
input.plasma.Pabs=input.plasma.Pin*source.PTE;

%     flag.global_model='H2_LZS2020';
[Res]=myOde2(t,X,input.plasma.Pabs,input.plasma.p,...
    input.plasma.Tg,input.plasma.Q0,input.plasma.r,input.plasma.l);
end
