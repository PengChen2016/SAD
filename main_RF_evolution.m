% -*- coding: utf-8 -*-
% ----------------------------------------------
%{
 * brief #Abstract
 * Created 13:10:30 2022/11/07
 * author
   Zengshan Li, HUST: origin 新的可用的整体模型代码.rar 2020/09/03
   Peng Chen, HUST: modified 整合进SAD
 * note #Detail
 * #TODO
%}
% ----------------------------------------------

%% initialization
% close all
clear
addpath(genpath('./packages'))
now_str=datestr(now,'yyyymmdd_HHMMSS');

%% flag
solution_name='test_gloabl_model';
addpath(genpath('./packages'))
addpath(genpath(['./others/' solution_name '/']))
flag.using_stored_data=false;
% flag.using_stored_data=true;

%------------------------------ 1. test_new1
program_name='test_new1';

%%%%%%%% plasma model
flag.type_Xsec='e-H2-Phelps';

% stoc expression
% flag.stoc_model='';
% flag.stoc_model='Vahedi-simplify';
% flag.stoc_model='Cazzador-simplify';
% flag.stoc_model='2018Jainb-simplify';
flag.stoc_model='Cazzador-fit';

flag.medium_approximation='';
% flag.medium_approximation='consider_real(sigma)_only';
% flag.medium_approximation='sigma_dc';

flag.skin_depth='as-medium';
% flag.skin_depth='as-medium-simplified';
% flag.skin_depth='collisionless';
% flag.skin_depth='collision';
% flag.skin_depth='as-medium-simplified-finite-radius';

% flag.output_plasma_model=true;
flag.output_plasma_model=false;

%%%%%%%% electric model
flag.electric_model='';
% flag.electric_model='transformer-base';
% flag.electric_model='analytical_base';

flag.input_geometry='BUG_base';

flag.Rmetal='measured-Rmetal-woplasma';
% flag.Rmetal='calculated-Rcoil-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';

flag.output_electric_model=true;
% flag.output_electric_model=false;
%------------------------------ 1. flag end

save_mat_name=['./others/' solution_name '/' program_name '.mat'];
if ~flag.using_stored_data
    % calculate
    log_name=['./others/' solution_name '/' program_name '.log'];
    diary(log_name) % append to the end of the log file
    fprintf('\n-----%s %s-----\n\n',program_name,now_str)
    
    %% Input
    % flag.global_model='';
    flag.global_model='H2_LZS2020';
    %%%%%%%%%%%% D1
    flag1=flag;
    flag1.input_plasma='BUG_0.3Pa_1MHz_55kW';
    flag1.electric_model='';
    input1=get_input_data( flag1 );
    % modify input   
    % 气体
    input1.plasma.Tg=600;                        %温度
    input1.plasma.p=0.25;                        %压强
    input1.plasma.Q0=35; % sccm
    % 几何
    input1.plasma.r=0.1;%小源参数0.06
    input1.plasma.l=0.3;%小源参数0.3
    % 电参数
    input1.plasma.Pin=2000;
    
    input1.X.tspan=[0 1]; % 求解时间区间
    % 等离子体参数 初值
    input1.X.Te=2;
    input1.X.nHi=2E15;
    input1.X.nH2i=4E15;
    input1.X.nH3i=4E15;
    input1.X.nHNi=1E14;
    input1.X.ne=input1.X.nHi+input1.X.nH2i+input1.X.nH3i-input1.X.nHNi;
    input1.X.nH=3E17;
    input1.X.nHH=30E19;

    % %%%%%%%%%%%% D2
    % flag2=flag;
    % flag2.input_plasma='zero';
    % plasma2.size=1;
    % plasma2.f=1e6;
    % plasma2.ne=1.5e14;
    % plasma2.Te=2;
    % plasma2.p=0.3;
    % plasma2.Tg=630;
    % plasma2.Pin=[];
      

    %% solving
    tic
    [Xt]=plasma_evolution(flag1,input1);
    toc
    
    save(save_mat_name)
else
    % load data
    load(save_mat_name)
    disp(now_str)
    warning(['Using data stored in ' save_mat_name])
    % % %     if flag.output_plasma_model
    % % %         output_plasma_model(flag,input.plasma)
    % % %     end
    % % %     if flag.output_electric_model
    % % %         output_electric_model( flag, source )
    % % %     end
end

%% post-processing

figure
subplot(1,2,1);
plot(Xt.t,Xt.Te)
ylabel('{\itT}_{e} [eV]')
xlabel('{\itt} [s]');
grid on
L1=legend('{\itT}_{e}');
set(L1,'location','south');
set(L1,'box','off')
set(L1,'AutoUpdate','off')

subplot(1,2,2);
semilogy(Xt.t,Xt.ne)
hold on
semilogy(Xt.t,Xt.nHH,'-.')
semilogy(Xt.t,Xt.nH,'--')
ylabel('{\itn} [m^{-3}]')
xlabel('{\itt} [s]');
grid on
L1=legend('{\itn}_{e}','{\itn}_{H_2}','{\itn}_{H}');
set(L1,'location','south');
set(L1,'box','off')
set(L1,'AutoUpdate','off')

figure
subplot(1,2,1);
semilogx(Xt.t,Xt.Te)
ylabel('{\itT}_{e} [eV]')
xlabel('{\itt} [s]');
grid on
L1=legend('{\itT}_{e}');
set(L1,'location','south');
set(L1,'box','off')
set(L1,'AutoUpdate','off')
axis([-inf,inf,-inf,inf])
subplot(1,2,2);
loglog(Xt.t,Xt.ne)
hold on
loglog(Xt.t,Xt.nHH,'-.')
loglog(Xt.t,Xt.nH,'--')
ylabel('{\itn} [m^{-3}]')
xlabel('{\itt} [s]');
grid on
L1=legend('{\itn}_{e}','{\itn}_{H_2}','{\itn}_{H}');
set(L1,'location','south');
set(L1,'box','off')
set(L1,'AutoUpdate','off')
axis([-inf,inf,-inf,inf])

%% End
if ~flag.using_stored_data
    fprintf('\n-----END %s-----\n\n',now_str)
    diary off
end
