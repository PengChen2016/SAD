%% Abstract
% RF-ICP源功率沉积模型，负源激励器功率耦合（激励器等效阻抗与RF传输效率）分析代码

close all
clear
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
% now_str=datestr(now,'yyyy-mm-dd');
disp(now_str)

%物理常量
e=1.6e-19;            %电子电量
me=9.1e-31;           %电子质量
mi=1.67e-27;          %氢离子质量
eps0=8.85e-12;          %真空介电常数
c=3e8;                %真空光速[m/s]
kB=1.381e-23;         %玻尔兹曼常数 J/K
mu0=4*pi*1e-7;         %真空磁导率
sigma_Cu=5.8e7;       %铜的电导率

%% 输入
%%%% 控制位
flag_parametric_analysis=false;
flag_parametric_analysis=true; % 参数化分析则为真

flag_parameters_coupling=false;
flag_using_stored_data=false;
if flag_parametric_analysis
    %             flag_parameters_coupling=true; % 不同参数一一对应分析则为真
    %         flag_using_stored_data=true; % 使用存储数据则为真
end

% 实验数据、输入条件来源
% flag_experiment_data='ELISE_base'; %Ref from Jain by CP
% flag_experiment_data='BATMAN_base'; %no data
% flag_experiment_data='HUST_small_driver_base'; %by ZP
% flag_experiment_data='CHARLIE_base'; %Ref from Rauner by CP
% flag_experiment_data='small_source1_LZS'; %Ref from 李增山-整体模型耦合解析电模型-2020.03.30\ by CP
% flag_experiment_data='NIO1_base'; %no data
flag_experiment_data='HUST_large_driver_base'; %by CP, LJW 201029

flag_Rmetal='experiment-measured-woplasma';
% flag_Rmetal='theory-coil-woplasma';
flag_Lcoil='experiment-measured';
% flag_Lcoil='theory-coil';

% stoc表达式
% flag_vst_expression='Vahedi-simplify';
% flag_vst_expression='Cazzador-simplify';
flag_vst_expression='Cazzador-fit';
flag_electric_model='transformer_base';
% flag_electric_model='analytical_base';
% flag_electric_model='transformer_2011Chabert';

flag_good_conductor_approximation=false;
% flag_good_conductor_approximation=true;

flag_output_plasma_model=true;
flag_output_plasma_model=false;
flag_output_electric_model=true;
% flag_output_electric_model=false;
flag_output_for_paper=true;
flag_output_for_paper=false;
if flag_output_for_paper
    flag_experiment_data='ELISE_base';
    flag_parametric_analysis=true; 
%     flag_parametric_analysis=false;
    
    flag_parameters_coupling=false;
    flag_using_stored_data=false;
    
    flag_vst_expression='Cazzador-fit';
            flag_vst_expression='Vahedi-simplify'; %用于对比st模型
    flag_electric_model='transformer_base';
    flag_good_conductor_approximation=false;
    
%     flag_experiment_data='CHARLIE_base';
%     flag_parametric_analysis=false; 
end

%%%% 等离子体与ICP源参数
fprintf('Used experiment data: %s \n',flag_experiment_data);
if flag_parameters_coupling && ~strcmp(flag_experiment_data,'CHARLIE_base')
    error('no data')
end
% input of  equivalent_medium_model_of_plasma
%%%%%%%%%%%%%%%%%%%%%%%%%%%等离子体参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=1e6;               %驱动频率，单位Hz
if ~flag_parametric_analysis
    %用于单点计算
    switch flag_experiment_data
        case 'ELISE_base'
            %2018Jain的ELISE case使用5e18，15eV，1200-1500K
            ne=5e18;                 %电子密度[m^-3]
            %             Te=15;                     %电子温度[eV] 2018Jain
            Te=9;                     %电子温度[eV] 2018Jainb
            p=0.3;                      %气压[Pa]
            Tg=1200;                     %气体温度[K] %2018Jainb
            %             Tg=630;                     %气体温度[K] %2018Fantz
        case 'BATMAN_base'
            %2010Mcneely中为中心1.5e18，>10eV
            ne=1.5e18;                 %电子密度[m^-3]
            Te=10;                     %电子温度[eV]
            %2006Fantz报道激励器内1000K，2018Fantz报道整体630K
            p=0.3;                      %气压[Pa]
            Tg=630;                     %气体温度[K]
        case 'HUST_large_driver_base'
             ne=5e17;                 %电子密度[m^-3]
            Te=9;                     %电子温度[eV]
            %2006Fantz报道激励器内1000K，2018Fantz报道整体630K
            p=0.3;                      %气压[Pa]
            Tg=630;                     %气体温度[K]
        case 'HUST_small_driver_base'
            warning('no data')
        case 'CHARLIE_base'
            %             %2019Rauner中CHARLIE为4.2eV，平均密度1.8e17，r=40mm处密度4.8e16
            %             ne=4.8e16;                 %电子密度[m^-3]
            %             Te=4.2;                     %电子温度[eV]
            
            %200531 2019Raunera
            %             %axial center radial center
            %             ne=1.7e17;                 %电子密度[m^-3]
            %             Te=5.5;                     %电子温度[eV]
            %axial average radial center
            ne=7e16;                 %电子密度[m^-3]
            Te=5.5;                     %电子温度[eV]
            %             %axial center radial edge
            %             ne=9.2e16;                 %电子密度[m^-3]
            %             Te=5.5;                     %电子温度[eV]
            
            %2019Rauner中的一组参数
            %         p=1;                      %气压[Pa]
            %         Tg=560;                     %气体温度[K]
            p=0.3;                      %气压[Pa]
            Tg=128.7*log10(p)+568.5;                     %气体温度[K]
            
            if flag_output_for_paper
                                %论文使用 10Pa，520W，4MHz
                p=10;
                Te=2.1;
                f=4e6;
                ne=4.1e17;
                I_rms=9.66;
                Rmetal_ex=0.21;
                Tg=128.7*log10(p)+568.5;
            end
            
        case 'small_source1_LZS'
            ne=1.6e17;
            Te=6.5;
            p=0.6;
            Tg=600;
        case 'NIO1_base'
            warning('no data')
            f=2.1e6;
            p=1;
            Tg=400;
        otherwise
            error('Unexpected experiment_result.')
    end
else
    % 注意后文进行两处相应修改
    switch flag_experiment_data
        case 'BATMAN_base'
            % 聚变负源典型参数范围
            ne=[1e17,5e17,1e18];                 %电子密度[m^-3]
            %以下仅为画图好看。实际上这里是温度不是能量，聚变负源有100eV的电子但没有100eV的温度
            Te=[0.1:0.2:0.9,1:1:100];                 %电子温度[eV]
            
            p=0.3;                      %气压[Pa]
            Tg=630;                     %气体温度[K]
        case 'ELISE_base'
            if flag_output_for_paper
                dne=16:0.1:19;             %指数
                ne=10.^dne;                 %电子密度[m^-3]
                Te=5:5:25;                 %电子温度[eV]
                
                %                 dne=16:1:19;             %指数
                %                 ne=10.^dne;                 %电子密度[m^-3]
                %                 Te=5:0.1:25;                 %电子温度[eV]
                
                p=0.3;                      %气压[Pa]
                Tg=630;                     %气体温度[K]
            else
                dne=16:0.1:19;             %指数
                ne=10.^dne;                 %电子密度[m^-3]
                %             ne=[5e17,5e18,1e19];
                
                %             Te=10;
                %             Te=9; % 2018Jainb使用
                Te=1:4:21;                 %电子温度[eV]
                
                p=0.3;                      %气压[Pa]
                %             p=0.3:0.3:10; % 2019Raunera使用0.3~10Pa
                %                          p=0.1:0.1:1; % 2018Jainb使用0.1~1Pa
                
                %Tg=1200;                     %气体温度[K]
                Tg=630;                     %气体温度[K]
            end
        case 'CHARLIE_base'
            if flag_parameters_coupling
                %参数耦合情况:一一对应关系
                % 注意后文进行两处相应修改
                %                 p=[0.3;0.5;1;3;5;10];
                %                 Te=[5.50 ;4.79 ;4.06 ;3.00 ;2.62 ;2.10];
                %                 f=1e6;
                %                 ne=[7.0E+16;1.3E+17;1.7E+17;2.4E+17;2.4E+17;2.2E+17];
                %                 I_rms=[51.94 ;35.45 ;28.19 ;28.35 ;32.69 ;44.23];
                %                 Rloss=[0.0970835;0.0955728;0.0889425;0.0892782;0.0940621;0.08953];
                %                 f=4e6;
                %                 ne=[1.2E+17;1.6E+17;2.3E+17;2.7E+17;3.3E+17;3.7E+17];
                %                 I_rms=[17.55 ;14.07 ;12.21 ;10.02 ;9.26 ;9.66];
                %                 Rloss=[0.245888;0.209715;0.180256;0.160029;0.183613;0.213659];
                %论文使用 10Pa，520W，4MHz
                p=10*ones(8,1);
                Te=2.1*ones(8,1);
                f=4e6;
                ne=[8.26E+16;1.82E+17;2.22E+17;2.59E+17;2.86E+17;6.31E+17;7.68E+17;4.1e17];
                I_rms=9.66*ones(8,1);
                Rmetal_ex=0.21*ones(8,1);
            else
                % 参数解耦情况
                ne=[5e16,1e17,5e17,1e18,5e18,1e19,5e19];                 %电子密度[m^-3]
                
                Te=[10,15,20,25,30];                 %电子温度[eV]
                %Te=10;                 %电子温度[eV]
                
                p=0.3;                      %气压[Pa]
            end
            Tg=128.7*log10(p)+568.5;                     %气体温度[K]
        case 'NIO1_base'
            % 2014Cazzador使用参数范围
            ne=1.8e16;                 %电子密度[m^-3]
            Te=[0.1:0.2:0.9,1:3:100];                 %电子温度[eV]
        case 'HUST_large_driver_base'
                dne=14:0.1:18;             %指数
                ne=10.^dne;                 %电子密度[m^-3]
                Te=6:2:10;                 %电子温度[eV]                
                p=0.3;                      %气压[Pa]
                Tg=630;                     %气体温度[K]
        otherwise
            error('Unexpected experiment_result.')
    end
end

%字符数组不能够存储不同长度的字符串
X_var={'ne';'Te';'p'};
name_X_var={'\itn\rm_e';'\itT\rm_e';'\itp'};
unit_X_var={'m^{-3}';'eV';'Pa'};
if flag_parameters_coupling
    %参数耦合情况:一一对应关系
    %或仅单个画图参数？？
    %     X1=p;
    %     idx_X1=3;
    %     num_X1=length(p);
    %     num_X2=1;
    X1=ne;
    idx_X1=1;
    num_X1=length(ne);
    num_X2=1;
else
    %参数解耦情况
    %后处理与主次自变量解耦合代码
    %若前文已给定输入，在这里的两处修改，即可实现主次变量的改变
    
    %%%%%%%%%%%%%%%手动修改1%%%%%%%%%%%%%%%
    % ne-Te
    X1=ne; %绘图时第一自变量，即横轴
    idx_X1=1;
    X2=Te; %绘图时第二自变量，即legend
    idx_X2=2;
    %     X1=Te;
    %     idx_X1=2;
    %     X2=ne;
    %     idx_X2=1;
    X3=p; %不参与绘图的其他变量，实际上是单点值
    idx_X3=3;
    
    %     % p-ne
    % %     X1=p; %绘图时第一自变量，即横轴
    % %     idx_X1=3;
    % %     X2=ne; %绘图时第二自变量，即legend
    % %     idx_X2=1;
    %     X1=ne;
    %     idx_X1=1;
    %     X2=p;
    %     idx_X2=3;
    %     X3=Te; %不参与绘图的其他变量，实际上是单点值
    %     idx_X3=2;
    %%%%%%%%%%%%%%%手动修改1%%%%%%%%%%%%%%%
    
    %紧跟着的后续运算
    %         fprintf('Used multi X: %s \n',[X_var{idx_X1} ' & ' X_var{idx_X2}]);%弃用
    num_X1=length(X1);
    num_X2=length(X2);
    no_mid_X2=ceil(num_X2/2);
    
    if length(X3)>1
        error([X_var{idx_X3} ' must be single value'])
    end
    
    
    clear ne Te
    for X1i=1:num_X1
        for X2i=1:num_X2
            %%%%%%%%%%%%%%%手动修改2%%%%%%%%%%%%%%%
            % ne-Te
            ne(X1i,X2i)=X1(X1i);
            Te(X1i,X2i)=X2(X2i);
            %                                     Te(X1i,X2i)=X1(X1i);
            %                                     ne(X1i,X2i)=X2(X2i);
            p(X1i,X2i)=X3;
            
            %             % p-ne
            % %             p(X1i,X2i)=X1(X1i);
            % %             ne(X1i,X2i)=X2(X2i);
            %                     ne(X1i,X2i)=X1(X1i);
            %                     p(X1i,X2i)=X2(X2i);
            %             Te(X1i,X2i)=X3;
            %%%%%%%%%%%%%%%手动修改2%%%%%%%%%%%%%%%
            legend_X2{X2i}=[name_X_var{idx_X2} '\rm=' num2str(X2(X2i),'%.1e') '\rm' unit_X_var{idx_X2}];
        end
    end
end

w_RF=2*pi*f; %驱动角频率，单位Hz
ng=p./(kB*Tg);  %由理想气体状态方程得到中性气体分子密度[m^-3]，用于碰撞频率计算

% input of electric model
%%%%%%%%%%%%%%%%%%%%%%%%%%%几何参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch flag_experiment_data
    case 'ELISE_base'
        r_chamber=0.14;             %放电腔内径[m]
        r_plasma=r_chamber; %简化：以放电腔半径为等离子体半径
        r_coil=0.158;             %线圈绕线半径[m]
        r_wire_of_coil=0.004;             %线圈导线半径[m] %200625改正，之前为0.003
        N_coil=6;                  %线圈匝数
        l_chamber=0.14;               %线圈、等离子体长度[m]
        l_plasma=l_chamber;
        l_coil=0.08;
        % 一次线圈为短螺线管，二次线圈为长-短螺线管
    case 'HUST_large_driver_base'
        r_chamber=0.142;             %放电腔内径[m]
        r_plasma=r_chamber; %以放电腔半径为等离子体半径
        r_coil=0.155;             %线圈绕线半径[m]
        r_wire_of_coil=0.003;             %线圈导线半径[m] 
        N_coil=9;                  %线圈匝数
        l_chamber=0.147;               %放电腔长度[m] 
        l_plasma=l_chamber; %以放电腔长度为等离子体长度
        l_coil=0.08935;   %线圈长度     
    case 'BATMAN_base'
        warning('no data')
        %TODO: 待修改。目前是copy ELISE数据
        r_chamber=0.14;             %放电腔内径[m]
        r_plasma=r_chamber; %简化：以放电腔半径为等离子体半径
        r_coil=0.158;             %线圈绕线半径[m]
        r_wire_of_coil=0.003;             %线圈导线半径[m]
        N_coil=6;                  %线圈匝数
        l_chamber=0.14;               %线圈、等离子体长度[m]
        l_plasma=l_chamber;
        l_coil=0.08;
        % 一次线圈为短螺线管，二次线圈为长-短螺线管
    case 'HUST_small_driver_base'
        r_chamber=0.051;             %放电腔内径[m]
        r_plasma=r_chamber; %简化：以放电腔半径为等离子体半径
        r_coil=0.063;             %线圈绕线半径[m]
        r_wire_of_coil=0.003;             %线圈导线半径[m]
        N_coil=6;                  %线圈匝数
        l_chamber=0.06;               %线圈、等离子体长度[m]
        l_plasma=l_chamber;
        l_coil=l_chamber;
    case 'CHARLIE_base'
        r_chamber=0.05;             %放电腔内径[m]
        r_plasma=r_chamber-0.005; %简化：以放电腔半径为等离子体半径
        r_coil=0.055;             %线圈绕线半径[m]
        r_wire_of_coil=0.003;             %线圈导线半径[m]
        N_coil=5;                  %线圈匝数
        l_chamber=0.4;               %线圈、等离子体长度[m]
        l_plasma=l_chamber;
        %         l_plasma=0.2; %几何校正
        % 当Lmp≥Rp/veff，可见Rp对PER影响小，则lplasma对PER影响小
        l_coil=0.1;
    case 'small_source1_LZS'
        r_chamber=0.06;             %放电腔内径[m]
        r_plasma=r_chamber; %简化：以放电腔半径为等离子体半径
        r_coil=r_chamber;             %线圈绕线半径[m]
        r_wire_of_coil=0.003;             %线圈导线半径[m]
        N_coil=5;                  %线圈匝数
        l_chamber=0.16;               %线圈、等离子体长度[m]
        l_plasma=l_chamber;
        l_coil=l_chamber;
    case 'NIO1_base'
        warning('no data')
        %TODO: 待修改。目前是copy ELISE数据
        r_chamber=0.14;             %放电腔内径[m]
        r_plasma=r_chamber; %简化：以放电腔半径为等离子体半径
        r_coil=0.158;             %线圈绕线半径[m]
        r_wire_of_coil=0.003;             %线圈导线半径[m]
        N_coil=6;                  %线圈匝数
        l_chamber=0.14;               %线圈、等离子体长度[m]
        l_plasma=l_chamber;
        l_coil=0.08;
        % 一次线圈为短螺线管，二次线圈为长-短螺线管
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%电参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 理论计算线圈阻抗
disp('理论计算线圈阻抗')
% TODO:考虑FEM中引线长度，计算一个理论电阻
% TODO: 验算以下公式是否在适用范围内
Rcoil_th=N_coil*2*pi*r_coil*sqrt(mu0*w_RF*sigma_Cu/2)/2/pi/r_wire_of_coil/sigma_Cu;%已考虑集肤效应，未考虑邻近效应
Lcoil_th=mu0*check_Nagaoka(2*r_coil/l_coil)*pi*r_coil^2*N_coil^2/l_coil;  %一次螺线管线圈电感理论计算表达式
fprintf('%s = %.2e MHz\n','f',f/1e6);
fprintf('%s = %.2e \n','Rcoil_th',Rcoil_th);
fprintf('%s = %.2e \n','Lcoil_th',Lcoil_th);
% 实验测量线圈阻抗（一般按照相减方法来测量）
Icoil_rms=1;  %射频电流[A] 计算功率绝对值时使用――因此需要根据实验测得值输入
switch flag_experiment_data
    case 'ELISE_base'
        % 2018Jain
        Rmetal_ex=0.5;
        Lcoil_ex=7.5e-6;
        % 2018Jainb在变压器模型中加入了一个表征其他金属损耗的电阻
    case 'HUST_large_driver_base'
        Rmetal_ex=1.45;
        Lcoil_ex=16.04e-6;
    case 'BATMAN_base'
        warning('no data')
        %TODO: 待修改。目前是copy ELISE数据
        Rmetal_ex=0.5;
        Lcoil_ex=7.5e-6;
    case 'HUST_small_driver_base'
        Rmetal_ex=0.8;
        Lcoil_ex=18.7e-6;
    case 'CHARLIE_base'
        %         %1MHz,0.3Pa
        %         Rcoil=0.1;
        %         I_coil=51.9;  %射频电流[A]
        
%         %2019Raunera
%         %4MHz,0.3Pa
%         Rmetal_ex=0.25;
%         Icoil_rms=17.5;  %射频电流[A]
        
        Lcoil_ex=2.2e-6;
        if flag_parameters_coupling
            Icoil_rms=I_rms;
        end
    case 'small_source1_LZS'
        Rmetal_ex=1;
        Lcoil_ex=2.2e-6; %原文未给出，不关心
    case 'NIO1_base'
        warning('no data')
        %TODO: 待修改。目前是copy ELISE数据
        Rmetal_ex=0.5;
        Lcoil_ex=7.5e-6;
end

switch flag_Rmetal
    case 'experiment-measured-woplasma'
        Rmetal=Rmetal_ex;
        disp('基于Subtractive method，使用实验测得Rs wo plasma作为Rmetal')
        disp('有可能部分实验数据是Rcoil而非Rs')
    case 'theory-coil-woplasma'
        Rmetal=Rcoil_th;
        disp('使用理论计算Rcoil wo plasma作为Rmetal')
end

switch flag_Lcoil
    case 'experiment-measured'
        Lcoil=Lcoil_ex;
        disp('使用实验测得的Lcoil')
case 'theory-coil'
Lcoil=Lcoil_th;
        disp('使用理论计算的Lcoil')
end

Qcoil=w_RF.*Lcoil./Rmetal;

if ~flag_parameters_coupling
    X_temp1=Rmetal;
    X_temp2=Icoil_rms;
    for X1i=1:num_X1
        for X2i=1:num_X2
            Rmetal(X1i,X2i)=X_temp1;
            Icoil_rms(X1i,X2i)=X_temp2;
        end
    end
end

if ~flag_using_stored_data
    %% ICP的等效电磁媒质模型 equivalent_medium_model_of_plasma
    %%%% 如果改成函数，虽然简洁不需要多处修改了，但传递数据似乎比较麻烦？不能改成函数，顶多多脚本文件
    % 主要基于2014Cazzador、1995Vahedia，借鉴2018Jain、
    disp('## ICP的等效电磁媒质模型')
    % aid function
    % classical skin depth
    % wpe>>w,v时，不考虑有限半径的经典集肤深度
    delta_simplified_fun=@(vc,wpe)(c/wpe)*sqrt(...
        2*(1+vc^2/w_RF^2)/(1+sqrt(1+vc^2/w_RF^2)));
    % wpe>>w,v时，考虑有限半径的经典集肤深度。适用范围存疑
    b_geo=@(vc,wpe)1+(1+vc^2/w_RF^2)*... %几何效应因子 by 1995Vahedi
        (c/wpe)^2*((3.83/r_chamber)^2-(w_RF/c)^2); %a_geo
    delta_geo_simplified_fun=@(vc,wpe)(c/wpe)*sqrt(...
        2*(1+vc^2/w_RF^2)/b_geo(vc,wpe)/...
        (1+sqrt(1+vc^2/w_RF^2/b_geo(vc,wpe)^2)));
    % 经典集肤深度
    A_delta=@(vc,wpe)mu0*eps0*w_RF^2*(1-wpe^2/(vc^2+w_RF^2));
    B_delta=@(vc,wpe)mu0*eps0*vc*w_RF*wpe^2/(vc^2+w_RF^2);
    delta_fun=@(vc,wpe)(2/B_delta(vc,wpe))*sqrt(...
        (A_delta(vc,wpe)+sqrt(A_delta(vc,wpe)^2+B_delta(vc,wpe)^2))/2);
    
%     % 经典集肤深度随vc变化
%     d_vc=4:0.1:8;
%     num_x=length(d_vc);
%     vc=10.^d_vc;
%     for i=1:num_x
%         delta1(i)=delta_simplified_fun(vc(i),1e9);
%         delta2(i)=delta_geo_simplified_fun(vc(i),1e9);
%         delta3(i)=delta_fun(vc(i),1e9);
%     end
%     handle_fig=figure;
%     loglog(vc,delta1,'-');
%     hold on
%     loglog(vc,delta2,'--s');
%     loglog(vc,delta3,'-.o');
%     legend('simplified','geo simplified','full')
%     % 可见经典集肤深度基本与vc正相关。
%     % delta_geo_simplified_fun有问题
    
    % stoc model
    fprintf('vst use %s\n',flag_vst_expression);
    alpha_st_fun=@(delta_st,Ve)4*delta_st^2*w_RF^2/pi/(Ve^2);
    % vst Cazzador fit
    K0=me*me/(mu0*e^3); %与频率无关 %Te取eV，不需要再乘以e
    f_Cazzador=2.1e6; %驱动频率，单位Hz
    w_RF_Cazzador=2*pi*f_Cazzador;
    vst_Cazzador_fit_fun=@(x)10^(-244.1+48.1*log10(x)-3.467*log10(x)^2+...
        0.1113*log10(x)^3-0.001336*log10(x)^4); % 2.1MHz下的vst-x关系，x=neTe
    va_Cazzador_fit_fun=@(X)vst_Cazzador_fit_fun(K0*w_RF_Cazzador^2*X)/w_RF_Cazzador; %与频率无关的va-X关系
    vst_new_f_fit_fun=@(x)va_Cazzador_fit_fun(x/K0/w_RF^2)*w_RF;
    vst_new_f_fit_fun1=@(x)vst_Cazzador_fit_fun(x*w_RF_Cazzador^2/w_RF^2)*w_RF/w_RF_Cazzador; %测试确定，该变换表达式正确
    % 测试确定，通过变换方法得到的结果合理
    
    for X1i=1:num_X1
        for X2i=1:num_X2
            % 考虑在这里传递数据，即将数组数据取出作为单个变量，并在循环结束后再赋给数组；以方便解耦合
            wpe(X1i,X2i)=sqrt(ne(X1i,X2i)*e^2/me/eps0);     %电子等离子体频率
            wpi(X1i,X2i)=sqrt(ne(X1i,X2i)*e^2/mi/eps0);     %离子等体频率
            Ve(X1i,X2i)=sqrt(8*Te(X1i,X2i)*e/(pi*me));    %电子平均速率，计算自由程用
            %             wce(X1i,X2i)=e;     %电子等离子体频率
            
            %以下求解弹性碰撞频率vm
            kenp(X1i,X2i)=k_enp(Te(X1i,X2i));
            venp(X1i,X2i)=ng(X1i,X2i)*kenp(X1i,X2i); %电子与中性分子动量转移碰撞
            lambda_coll_enp=Ve(X1i,X2i)/venp(X1i,X2i); %平均自由程
            keniz(X1i,X2i)=k_eniz(Te(X1i,X2i));
            veniz(X1i,X2i)=ng(X1i,X2i)*keniz(X1i,X2i); %电子与中性分子电离碰撞
            lambda_coll_eniz=Ve(X1i,X2i)/veniz(X1i,X2i); %平均自由程
            A(X1i,X2i)=4*pi*mi*(eps0*e*Te(X1i,X2i)^(3/2))/(e^3)/(me+mi)/sqrt(ne(X1i,X2i)); %库仑算子中Λ
            veip(X1i,X2i)=ne(X1i,X2i)*4*sqrt(2*pi)/3*((e^2/4/pi/eps0)^2)*log(A(X1i,X2i))/sqrt((e*Te(X1i,X2i))^3*me); %库仑碰撞
            lambda_coll_eip=Ve(X1i,X2i)/veip(X1i,X2i); %平均自由程
            vm(X1i,X2i)=venp(X1i,X2i)+veniz(X1i,X2i)+veip(X1i,X2i);
            delta_m(X1i,X2i)=delta_simplified_fun(vm(X1i,X2i),wpe(X1i,X2i));
            
            %以下求解随机加热频率vst
            switch flag_vst_expression
                case 'Cazzador-fit'
                    X_temp=ne(X1i,X2i)*Te(X1i,X2i)*w_RF_Cazzador^2/w_RF^2;
                    if X_temp<1e14||X_temp>1e21
                        warning(['X=' num2str(X_temp) ' out of the range of Cazzador-fit'])
                        pause
                    end
                    vst(X1i,X2i)=vst_new_f_fit_fun(ne(X1i,X2i)*Te(X1i,X2i));
                    %                     test_temp=vst_new_f_fit_fun1(ne(X1i,X2i)*Te(X1i,X2i));
                    % 测试驱动频率对Cazzador fit关系的影响
                    %         vst4_bug=vst_fit_Cazzador(ne(X1i,X2i)*Te(X1i,X2i)); %有较大差别
                    delta_st(X1i,X2i)=delta_simplified_fun(vst(X1i,X2i),wpe(X1i,X2i));
                    alpha_st(X1i,X2i)=alpha_st_fun(delta_st(X1i,X2i),Ve(X1i,X2i));
                case 'Vahedi-simplify' %case 'Cazzador-fit'时，欲测试stoc则注释掉，顺序执行两种vst计算
                    % 基于1995Vahedia model，略有修正
                    delta_st_Vahedi1=(c^2*Ve(X1i,X2i)/(wpe(X1i,X2i)^2)/pi/w_RF)^(1/3); %趋肤深度
                    delta_st_Vahedi2=c/wpe(X1i,X2i);
                    alpha_st_Vahedi1=4*delta_st_Vahedi1^2*w_RF^2/pi/(Ve(X1i,X2i)^2); %表征电子穿过趋肤层耗时与RF周期比值的参数
                    alpha_st_Vahedi2=4*delta_st_Vahedi2^2*w_RF^2/pi/(Ve(X1i,X2i)^2);
                    if alpha_st_Vahedi1<0.03
                        vst(X1i,X2i)=Ve(X1i,X2i)/(2*pi*delta_st_Vahedi1); %190922PC，根据1995Vahedia
                        % 测试stoc model
                        alpha_st(X1i,X2i)=alpha_st_Vahedi1;
                        delta_st(X1i,X2i)=delta_st_Vahedi1;
                        %             vst2(X1i,X2i)=(4/pi)^0.2*(w_RF^0.4)*(Ve(X1i,X2i)/delta_st_Vahedi1)^0.6; %待询问2019Zuo
                        %             vst3(X1i,X2i)=-pi*Ve(X1i,X2i)/4/delta_st_Vahedi1/(log(alpha_st_Vahedi1)+1.58); %待询问2019Zuo
                        %             fprintf('%s = %.2e \n','vst',vst(X1i,X2i));
                        %             fprintf('%s = %.2e \n','vst2',vst2(X1i,X2i));
                        %             fprintf('%s = %.2e \n','vst3',vst3(X1i,X2i));
                    else
                        if alpha_st_Vahedi2<0.03
                            warning('α分段存在问题')
                            %                 pause %测试时使用
                            % 200422 ne=1e16,Te=15/20两次,出现bug
                            % 仍然使用α>0.03时的δst算出来的α
                            vst(X1i,X2i)=Ve(X1i,X2i)/delta_st_Vahedi2/4;
                        else
                            if alpha_st_Vahedi2<10
                                vst(X1i,X2i)=Ve(X1i,X2i)/delta_st_Vahedi2/4;
                            else
                                vst(X1i,X2i)=pi/4/(w_RF^2)*(Ve(X1i,X2i)/delta_st_Vahedi2)^3;
                            end
                        end
                        alpha_st(X1i,X2i)=alpha_st_Vahedi2;
                        delta_st(X1i,X2i)=delta_st_Vahedi2;
                    end
                case 'Cazzador-simplify'
                    % 基于2018Jain model，略有修正。copy LZS代码
                    %         试探逻辑与'Vahedi-simplify'类似，待修改
                    %
                    %         delta_st_Vahedi1=(c^2*Ve(X1i,X2i)/(wpe(X1i,X2i)^2)/pi/w_RF)^(1/3); %趋肤深度
                    %         delta_st_Vahedi2=c/wpe(X1i,X2i);
                    %         alpha_st_Vahedi1=4*delta_st_Vahedi1^2*w_RF^2/pi/(Ve(X1i,X2i)^2); %表征电子穿过趋肤层耗时与RF周期比值的参数
                    %         alpha_st_Vahedi2=4*delta_st_Vahedi2^2*w_RF^2/pi/(Ve(X1i,X2i)^2);
                    %         if alpha_st_Vahedi1<0.03
                    %             vst(X1i,X2i)=Ve(X1i,X2i)/(2*pi*delta_st_Vahedi1); %190922PC，根据1995Vahedia
                    %             delta_st(X1i,X2i)=delta_st_Vahedi1;
                    %             %             vst2(X1i,X2i)=(4/pi)^0.2*(w_RF^0.4)*(Ve(X1i,X2i)/delta_st_Vahedi1)^0.6; %待询问2019Zuo
                    %             %             vst3(X1i,X2i)=-pi*Ve(X1i,X2i)/4/delta_st_Vahedi1/(log(alpha_st_Vahedi1)+1.58); %待询问2019Zuo
                    %             %             fprintf('%s = %.2e \n','vst',vst(X1i,X2i));
                    %             %             fprintf('%s = %.2e \n','vst2',vst2(X1i,X2i));
                    %             %             fprintf('%s = %.2e \n','vst3',vst3(X1i,X2i));
                    %         else
                    %             if alpha_st_Vahedi2<0.03
                    %                 warning('α分段存在问题')
                    %                 %                 pause %测试时使用
                    %                 % 200422 ne=1e16,Te=15/20两次,出现bug
                    %                 % 仍然使用α>0.03时的δst算出来的α
                    %                 vst(X1i,X2i)=Ve(X1i,X2i)/delta_st_Vahedi2/4;
                    %             else
                    %                 if alpha_st_Vahedi2<10
                    %                     vst(X1i,X2i)=Ve(X1i,X2i)/delta_st_Vahedi2/4;
                    %                 else
                    %                     vst(X1i,X2i)=pi/4/(w_RF^2)*(Ve(X1i,X2i)/delta_st_Vahedi2)^3;
                    %                 end
                    %                 delta_st(X1i,X2i)=delta_st_Vahedi2;
                    %             end
                    %         end
            end
            %case {'Vahedi-simplify','Cazzador-simplify'}时测试stoc模型用
            switch flag_vst_expression
                case {'Vahedi-simplify','Cazzador-simplify'}
                    vst_fit(X1i,X2i)=vst_new_f_fit_fun(ne(X1i,X2i)*Te(X1i,X2i));
                    delta_st_fit(X1i,X2i)=delta_simplified_fun(vst_fit(X1i,X2i),wpe(X1i,X2i));
                    alpha_st_fit(X1i,X2i)=alpha_st_fun(delta_st_fit(X1i,X2i),Ve(X1i,X2i));
            end
            %             %使用外部的vm、vst数据
            %             if flag_single_independent_variable
            %                 fprintf('%s = %.2e , ','ourcode，vm',vm);
            %                 fprintf('%s = %.2e \n','vst',vst);
            %                 switch flag_experiment_data
            %                     case 'ELISE_base'
            %                         % 2018Jain中间数据作为变压器输入，用于变压器模型benchmark
            %                         disp('使用2018Jain的vm和vst')
            %                         vm=5.6e6;
            %                         vst=4.3e7;
            %                     case 'small_source1_LZS'
            %                         % 李增山-整体模型耦合解析电模型-2020.03.30\中间数据作为输入，用于解析电模型benchmark
            %                         disp('使用李增山-整体模型耦合解析电模型-2020.03.30的vm和vst')
            %                         vm=1.1238e7;
            %                         vst=3.6537e7;
            %                 end
            %             end
            
            %考虑两种加热机制后
            veff(X1i,X2i)=vst(X1i,X2i)+vm(X1i,X2i);%有效碰撞频率
            
            % 等离子体等效电磁媒质参数
            epsp_r(X1i,X2i)=1-wpe(X1i,X2i)^2/w_RF/(w_RF-1i*veff(X1i,X2i)); %等离子体复相对介电常数,即epsp
            % 复介电常数epsc=eps0*epsp_r
            sigmap(X1i,X2i)=eps0*wpe(X1i,X2i)^2/(1i*w_RF+veff(X1i,X2i)); % 复电导率
            mu_p=mu0; %等离子体磁导率取为真空磁导率
            sigma_medium(X1i,X2i)=eps0*veff(X1i,X2i)*wpe(X1i,X2i)^2/(w_RF^2+veff(X1i,X2i)^2); %与Re(sigmap)一致
            eps_medium(X1i,X2i)=eps0*(1-wpe(X1i,X2i)^2/(w_RF^2+veff(X1i,X2i)^2)); %与Re(eps0*epsp_r)一致
            if flag_good_conductor_approximation
                % 复电导率忽略虚部，即epsp_real=1,epsp_imag不变
                disp('等离子体等效电磁媒质模型中使用良导体近似')
                epsp_r(X1i,X2i)=1+1i*imag(epsp_r(X1i,X2i));
                sigmap(X1i,X2i)=real(sigmap(X1i,X2i));
                eps_medium(X1i,X2i)=1;
            end
            
            epsp_r_real(X1i,X2i)=real(epsp_r(X1i,X2i)); %即eps_r
            epsp_r_imag(X1i,X2i)=imag(epsp_r(X1i,X2i));
            tandelta(X1i,X2i)=-epsp_r_imag(X1i,X2i)/epsp_r_real(X1i,X2i); %损耗正切，ε'-jε''，tanδ=ε''/ε'
            %         %检验表达式用-断点调试
            %         epsp_r_real(X1i,X2i)=1-wpe(X1i,X2i)^2/(veff(X1i,X2i)^2+w_RF^2);         %复介电常数实部
            %         epsp_r_imag(X1i,X2i)=-veff(X1i,X2i)*wpe(X1i,X2i)^2/(veff(X1i,X2i)^2+w_RF^2)/w_RF;         %复介电常数虚部
            %         sigma_from_epsp=-w_RF*eps0*epsp_r_imag(X1i,X2i); %与sigmap_real一致
            sigmap_real(X1i,X2i)=real(sigmap(X1i,X2i)); %即sigma
            sigmap_imag(X1i,X2i)=imag(sigmap(X1i,X2i));
            %         %检验表达式用-断点调试
            %         sigmap_real(X1i,X2i)=eps0*wpe(X1i,X2i)^2*veff(X1i,X2i)/(veff(X1i,X2i)^2+w_RF^2);
            %         sigmap_imag(X1i,X2i)=eps0*wpe(X1i,X2i)^2*(-w_RF)/(veff(X1i,X2i)^2+w_RF^2);
            %         eps_r_from_sigmap(X1i,X2i)=sigmap_imag(X1i,X2i)/w_RF/eps0+1; %与epsp_r_real一致
            %         tandelta_from_sigmap(X1i,X2i)=sigmap_real(X1i,X2i)/(w_RF*eps0*eps_r_from_sigmap(X1i,X2i)); %与tandelta一致
            
            % 电磁波分析
            k_p_wave(X1i,X2i)=w_RF*sqrt(mu_p*eps0*epsp_r(X1i,X2i)); %wave number
            %             gamma_wave=1i*k_p_wave(X1i,X2i); %传播常数
            
            % 等效集肤深度
            skin_depth_eff(X1i,X2i)=delta_simplified_fun(veff(X1i,X2i),wpe(X1i,X2i));
            if ~flag_parametric_analysis&&flag_output_plasma_model
                if r_plasma<3*skin_depth_eff(X1i,X2i)
                    warning('δ>≈R ，电磁波穿透等离子体，电场基本均匀,集肤深度概念与随机加热模型不适用')
                    pause
                end
            end
            %             alpha_wave=-imag(k_p_wave(X1i,X2i)); %衰减常数
            %             skin_depth_eff(X1i,X2i)=1/alpha_wave; %一般介质中经典集肤深度
            % 对比不同趋肤深度表达式
            if ~flag_parametric_analysis&&flag_output_plasma_model
                % 对比不同适用情况的集肤深度表达式
                fprintf('%s = %.2e\n','EM经典δ',delta_fun(veff(X1i,X2i),wpe(X1i,X2i)));
                disp('若wpe >> w,vc(碰撞频率)，集肤深度表达式可简化') %vm/vst<veff<<wpe
                fprintf('%s = %.2e , ','wpe',wpe(X1i,X2i));
                fprintf('%s = %.2e \n','veff',veff(X1i,X2i));
                fprintf('%s = %.2e \n','简化经典δ',delta_simplified_fun(veff(X1i,X2i),wpe(X1i,X2i)));
                %                 fprintf('%s = %.2e \n','平面线圈ICP源考虑几何效应的简化δ',delta_geo_simplified_fun(veff(X1i,X2i),wpe(X1i,X2i)));
                
                % 检验表达式：相同参数下，delta_fun结果应与后文skin_depth_general_medium一致
                %             % 使用良导体参数，验证一般介质集肤深度公式与良导体集肤深度公式，在良导体参数下一致
                %                             sigmap_real=sigma_Cu;
                %                             f=50;
                %                             w_RF=2*pi*f;
                %                             epsp_r_real=1;
                %                             epsp_r=epsp_r_real-1i*sigma_Cu/w_RF/eps0;
                %                             tandelta=sigma_Cu/w_RF/eps0;
                
                skin_depth_good_conductor=skin_depth(sigmap_real(X1i,X2i),f,1); %良导体集肤深度
                k_p_wave(X1i,X2i)=w_RF*sqrt(mu_p*eps0*epsp_r(X1i,X2i)); %wave number
                alpha_wave=-imag(k_p_wave(X1i,X2i)); %对于良导体，衰减常数≈相位常数
                beta_wave=real(k_p_wave(X1i,X2i)); %对于良导体，衰减常数≈相位常数
                %                 wavelength_wave=2*pi/beta_wave; %导体趋肤深度等于波长的2π倍
                %             % 检验表达式-断点调试
                %             % RF-ICP一般复介电常数实部为负值，因此表达式与常见表达式略有不同
                %             alpha_wave=w_RF*sqrt(-mu_p*eps0*epsp_r_real*(sqrt(1+tandelta^2)+1)/2); %衰减常数
                skin_depth_general_medium=1/alpha_wave;
                fprintf('%s = %.2e , ','skin depth: 良导体表达式',skin_depth_good_conductor);
                fprintf('%s = %.2e \n','一般介质中表达式',skin_depth_general_medium);
                % 使用良导体参数，结果一致为skin_depth(5.8e7,50,1)=0.0093
                % 使用ELISE等体参数，集肤深度良导体公式结果略大于一般介质公式结果，可忽略差别
            end
            
            %             %使用外部的delta_eff数据
            %             if ~flag_parametric_analysis&&flag_output_plasma_model
            %                 fprintf('%s = %.2e , ','ourcode，veff',veff);
            %                 fprintf('%s = %.2e \n','skin_depth_eff',skin_depth_eff);
            %                 switch flag_experiment_data
            %                     case 'ELISE_base'
            %                         % 2018Jain中间数据作为变压器输入，用于变压器模型benchmark
            %                         disp('使用2018Jain的δeff')
            %                         skin_depth_eff=8.3e-3;
            %                 end
            %             end
            
            %波长
            beta_wave=real(k_p_wave(X1i,X2i)); %相位常数
            wavelength_wave(X1i,X2i)=2*pi/beta_wave; %一般介质中波长
            %             % 对比不同波长表达式
            %             if ~flag_parametric_analysis&&flag_output_plasma_model
            %                 %             % 使用无损介质参数，以验证一般介质中波长公式适用于无损介质
            %                 %             f=1e6;
            %                 %             w_RF=2*pi*f;
            %                 %             epsp_r_real=1;
            %                 %             epsp_r=epsp_r_real;
            %                 %             tandelta=0;
            %
            %                 wavelength_lossless_medium=lambda(epsp_r_real(X1i,X2i),f,1); %无损介质中波长
            %                 k_p_wave(X1i,X2i)=w_RF*sqrt(mu_p*eps0*epsp_r(X1i,X2i)); %无损时 波数=相位常数
            %                 alpha_wave=-imag(k_p_wave(X1i,X2i)); %无损时衰减常数=0
            %                 beta_wave=real(k_p_wave(X1i,X2i)); %无损时 波数=相位常数
            % %                 %检验表达式用-断点调试
            % %                 % RF-ICP一般复介电常数实部为负值，因此表达式与常见表达式略有不同
            % %                 beta_wave=w_RF*sqrt(-mu_p*eps0*epsp_r_real*(sqrt(1+tandelta^2)-1)/2); %传播常数
            %                 wavelength_general_medium=2*pi/beta_wave; %见前文
            %                 fprintf('%s = %.2e , ','wavelength: 无损介质中',wavelength_lossless_medium);
            %                 fprintf('%s = %e \n','一般介质中',wavelength_general_medium);
            %                 % 使用无损介质参数，结果一致为lambda(1,1e6,1)=300
            %                 % 使用ELISE等离子体参数，波长一般介质公式结果约为无损介质公式结果的一半，考虑导电性后分布参数效应更明显
            %             end
            %             %测试时应在此设置断点
            
            if ~flag_parametric_analysis&&flag_output_electric_model
                if  wavelength_wave(X1i,X2i)<3*r_plasma
                    warning('λ＜≈R ，需要考虑相位变化,变压器模型不适用')
                    pause
                end
            end
            
            %分析等效碰撞频率中的主导加热机制
            vm_per_veff(X1i,X2i)=vm(X1i,X2i)/veff(X1i,X2i);
            
            % 分析带电粒子是否响应电磁场
            wpe_per_w(X1i,X2i)=wpe(X1i,X2i)/w_RF;
            wpi_per_w(X1i,X2i)=wpi(X1i,X2i)/w_RF;
            
            %分析是否可以忽略heating model中位移电流，使用复电导率
            veff_displacement_ratio(X1i,X2i)=wpe(X1i,X2i)^2/(w_RF*sqrt(w_RF^2+veff(X1i,X2i)^2));
            
            %分析是否频率非常低到可以使用直流电导率-忽略了复电导率虚部
            veff_per_w(X1i,X2i)=veff(X1i,X2i)/w_RF;
            sigmaeff_dc(X1i,X2i)=eps0*wpe(X1i,X2i)^2/veff(X1i,X2i);         %直流电导率
        end
    end
end

file_name='stored_data_test200529_2.mat';
% 已存储则注释掉下列语句
if flag_parametric_analysis && flag_using_stored_data
    warning('Store data first.')
    error('end')
end
% 在flag_using_stored_data=false时运行以下三行一次,然后将flag_using_stored_data改回true
% if flag_parametric_analysis
%     flag_using_stored_data=true;
%     save(file_name)
%     error('end')
% end

if flag_parametric_analysis && flag_using_stored_data
    load(file_name)
    warning('using stored data')
end
%%%%%%%%%%%%%%%%%%%%%%%%%% 后处理 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_output_plasma_model
    if ~flag_parametric_analysis
        %用于单点计算
        disp('等离子体参数')
        fprintf('%s = %e , ','ne',ne);
        fprintf('%s = %.2e \n','Te',Te);
        fprintf('%s = %.2e ,','f',f);
        fprintf('%s = %.2e , ','p',p);
        fprintf('%s = %.2e , ','Tg',Tg)
        fprintf('%s = %e \n','ng',ng);
        disp('特征频率')
        fprintf('%s = %.2e , ','vm',vm);
        fprintf('%s = %.2e , ','vst',vst);
        fprintf('%s = %.2e \n','veff',veff);
        fprintf('%s = %.2e , ','ω_RF',w_RF);
        fprintf('%s = %.2e , ','ωpe',wpe);
        fprintf('%s = %.2e \n','ωpi',wpi);
        disp('模型适用条件') %参考利伯曼教材
        disp('if wp/w_RF > 1, 带电粒子响应RF电场')
        fprintf('%s = %.2e , ','wpe/w_RF',wpe_per_w);
        fprintf('%s = %.2e \n','wpi/w_RF',wpi_per_w);
        if ~flag_parametric_analysis && wpi_per_w>1
            disp('1995Vahedia的ICP heating model要求wpi<w_RF,但当前等离子体参数集并不满足')
        end
        disp('if vm_per_veff << 1, 加热机制以随机加热为主.')
        fprintf('%s = %.2e \n','vm_per_veff',vm_per_veff);
        disp('')
        disp('复介电常数')
        disp(['epsp_r = ' num2str(epsp_r,'%.2e')])
        fprintf('%s = %.2e , ','eps_r',epsp_r_real);
        fprintf('%s = %.2e \n','tandelta',tandelta);
        
        disp('复电导率')
        disp(['sigmap = ' num2str(sigmap,'%.2e')])
        disp('if v/w_RF >> 1, v+jω≈v,可以使用sigma_dc表达式')
        fprintf('%s = %.2e \n','veff/w_RF',veff_per_w);
        fprintf('%s = %.2e \n','sigma_dc',sigmaeff_dc);
        disp('if wpe^2/(w_RF*sqrt(w_RF^2+v^2)) >> 1, then sigmap without jw*eps0 can be used.')
        fprintf('%s = %.2e \n','for v=veff,',veff_displacement_ratio);
        disp('')
        disp('特征长度')
        fprintf('%s = %.2e \n','debye length',sqrt(eps0*Te*e/(ne*e*e)));
        fprintf('%s = %.2e , ','effective_skin_depth',skin_depth_eff);
        fprintf('%s = %.2e \n','wavelength_from_epsD',wavelength_wave);
        fprintf('%s = %.2e \n','r_plasma',r_plasma);
        % TODO：待输出δstoc
        disp('特征时间')
        fprintf('%s = %.2e , ','趋肤层渡越时间τ',skin_depth_eff/Ve);
        fprintf('%s = %.2e \n','RF周期T',2*pi/w_RF);
        
        disp('碰撞分析')
        fprintf('%s = %.2e , ','venp',venp);
        fprintf('%s = %.2e , ','veniz',veniz);
        fprintf('%s = %.2e \n','veip',veip);
        fprintf('%s = %.2e , ','λenp',lambda_coll_enp);
        fprintf('%s = %.2e , ','λeniz',lambda_coll_eniz);
        fprintf('%s = %.2e \n','λeip',lambda_coll_eip);
    else
        %         % 输出到未打开的excel表
        %         filename = 'testdata.xlsx';
        % sheet = 2;
        % xlRange = 'E1';
        % xlswrite(filename,A,sheet,xlRange)
        
        if flag_output_for_paper
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output for paper
            plot_line_width=3;
            gca_line_width=1;
            marker_size=8;
            font_size=15;
            marker_indices=1:5:length(X1);
            left_color = [0 0 0];
            right_color =left_color;
            %             right_color = [0.8500    0.3250    0.0980];
            group_color_list={'r','b','g','c','m','k','y'};
            group_style_list={'-','--','-.',':'};
            
            flag_multi_Te=true;
%             flag_multi_Te=false;
            
            
            % 等离子体等效电磁参数
            name_Y='等离子体等效电磁参数';
            Y1=sigma_medium;
            name_Y1='\it{\bf\sigma}\rm';
            Y2=-eps_medium/eps0;
            name_Y2='-\it{\bf\epsilon}\rm_r';
            handle_fig=figure;
            set(handle_fig,'defaultAxesColorOrder',[left_color; right_color]);
            % marker会遮掩线型，而且也不好看
            flag_logy=true;
            %             flag_logy=false;
            if flag_logy
                % 左右轴取值：为了便于联系两个曲线族与相应坐标轴
                no_group=1;
                yyaxis right;
                loglog(X1,Y1(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                ylabel([name_Y1 ' [S/m]']);
                yyaxis left;
                loglog(X1,Y2(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                ylabel(name_Y2);
                hold on
                for no_group=2:num_X2
                    yyaxis right;
                    loglog(X1,Y1(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    yyaxis left;
                    loglog(X1,Y2(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                end
            else
            end
            xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
            set(gca,'FontSize',font_size)
            set(gca, 'LineWidth',gca_line_width)
            %             title([name_Y ' \rmat \rm' now_str]);
            grid on%显示网格
            yyaxis left;
            text(0.4*X1(1),0.5e5,'(b)','FontSize',font_size)
            
            % legend按默认顺序：先左边全部，再右边全部
            % 自定义legend
            legend_style_group = zeros(1,num_X2);
            legend_text_group=cell(1,num_X2);
            for no_group=1:num_X2
                legend_style_group(no_group) = plot(NaN,NaN,['s' group_color_list{no_group}],'MarkerSize',marker_size,'MarkerEdgeColor',group_color_list{no_group},'MarkerFaceColor',group_color_list{no_group});
                legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
            end
            L1=legend(legend_style_group, legend_text_group{:});
            set(L1,'FontSize',font_size);
            set(L1,'location','southeast');
            set(L1,'box','off')
            set(L1,'AutoUpdate','off')
     
            legend_style_group = zeros(1,2);
            legend_text_group={name_Y1,name_Y2};
            for no_group=1:2
                legend_style_group(no_group) = plot(NaN,NaN,[group_style_list{no_group} 'k'],'LineWidth',plot_line_width);
            end
            axes2 = axes('position',get(gca,'position'),'visible','off');
            L2=legend(axes2,legend_style_group, legend_text_group{:});
            set(L2,'FontSize',font_size);
            set(L2,'location','south');
            set(L2,'box','off')
            
            fprintf('%s = %.2e \n',[ num2str(X1(21)) ',15eV,σ= '],sigma_medium(21,3));
            fprintf('%s = %.2e \n',[ num2str(X1(11)) ',15eV,σ= '],sigma_medium(11,3));
            
                        
            name_Y='等离子体等效电磁参数比值';
            Y1=sigma_medium./(-w_RF*eps_medium);
            name_Y1='\it{\bf\sigma}/{\bf\omega}|{\bf\epsilon}|\rm';
            handle_fig=figure;
            flag_logy=true;
            %             flag_logy=false;
            if flag_logy
                for no_group=1:num_X2
                    loglog(X1,Y1(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    hold on
                end
                ylabel(name_Y1);
            else
            end
            xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
            set(gca,'FontSize',font_size)
            set(gca, 'LineWidth',gca_line_width)
            axis([X1(1),X1(end),1,2e1])
            line([X1(1),X1(end)],[1,1],'linestyle',':','linewidth',1*plot_line_width,'color','k');
            line([X1(1),X1(end)],[10,10],'linestyle',':','linewidth',1*plot_line_width,'color','k');
            grid on%显示网格
            text(0.4*X1(1),0.6e0,'(a)','FontSize',font_size)
            %             title([name_Y1 ' \rmat \rm' now_str]);
            legend_text_group=cell(1,num_X2);
            for no_group=1:num_X2
                legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
            end
            L1=legend(legend_text_group{:});
            set(L1,'FontSize',font_size);
            set(L1,'location','northwest');
            set(L1,'box','off')
            set(L1,'AutoUpdate','off')
            
            no_X2=no_mid_X2;
            if ~flag_parameters_coupling
                fprintf(['论文绘图用：' X_var{idx_X2} '=' num2str(X2(no_X2)) ' ' unit_X_var{idx_X2} '\n'])
            end
            
            % 碰撞频率
            %             name_Y='碰撞频率';
            %             handle_fig=figure;
            %             loglog(X1,venp(:,no_X2),'-r','LineWidth',plot_line_width);
            %             %                 axis equal
            %             %                 axis([X1(1),X1(end),1e4,1e8]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            %             hold on
            %             loglog(X1,veniz(:,no_X2),'--r','LineWidth',plot_line_width);
            %             loglog(X1,veip(:,no_X2),'-.r','LineWidth',plot_line_width);
            %             line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',0.5*plot_line_width,'color','k');
            %             loglog(X1,vm(:,no_X2),'-ob','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','b',...
            %                 'LineWidth',plot_line_width);
            %             loglog(X1,vst(:,no_X2),'-sb','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','b',...
            %                 'LineWidth',plot_line_width);
            %             L1=legend('\it{\bf\nu}\rm_{en}^{p}','\it{\bf\nu}\rm_{en}^{iz}','\it{\bf\nu}\rm_{ei}^{p}','{\it\bf\omega}','\it{\bf\nu}\rm_{m}','\it{\bf\nu}\rm_{st}');
            %
            %             set(L1,'FontSize',font_size);
            %             set(L1,'location','best');
            %             set(L1,'Orientation','horizontal');
            %             xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
            %             ylabel('\it{\bf\nu} \rm[Hz]');
            %             set(gca,'FontSize',font_size)
            %             set(gca, 'LineWidth',gca_line_width)
            %
            %             %         set(L1,'box','off')
            %             %             title(['碰撞频率 \rmat \rm' now_str]);
            %             grid on%显示网格
            
            name_Y='等效碰撞频率';
            handle_fig=figure;
            if flag_multi_Te
                for no_group=1:num_X2
                    loglog(X1,vm(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    hold on
                    loglog(X1,vst(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                end
                line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',1*plot_line_width,'color','k');
                
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('\it{\bf\nu} \rm[Hz]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                %             title(['碰撞频率 \rmat \rm' now_str]);
                grid on%显示网格
                axis([X1(1),X1(end),2e6,2e8])
                text(0.4*X1(1),1e6,'(a)','FontSize',font_size)
                
                % 自定义legend
                legend_style_group = zeros(1,num_X2);
                legend_text_group=cell(1,num_X2);
                for no_group=1:num_X2
                    legend_style_group(no_group) = plot(NaN,NaN,['s' group_color_list{no_group}],'MarkerSize',marker_size,'MarkerEdgeColor',group_color_list{no_group},'MarkerFaceColor',group_color_list{no_group});
                    legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
                end
                L1=legend(legend_style_group, legend_text_group{:});
                set(L1,'FontSize',font_size);
                set(L1,'location','northwest');
                set(L1,'box','off')
                %                 set(L1,'Orientation','horizontal');
                set(L1,'AutoUpdate','off')
                text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
                
                legend_style_group = zeros(1,2);
                legend_text_group={'\it{\bf\nu}\rm_{m}','\it{\bf\nu}\rm_{st}'};
                for no_group=1:2
                    legend_style_group(no_group) = plot(NaN,NaN,[group_style_list{no_group} 'k'],'LineWidth',plot_line_width);
                end
                axes2 = axes('position',get(gca,'position'),'visible','off');
                L2=legend(axes2,legend_style_group, legend_text_group{:});
                set(L2,'FontSize',font_size);
                set(L2,'location','north');
%                 set(L2,'Orientation','horizontal');
                set(L2,'box','off')
            else
            end
            
            flag_multi_Te=false;
            
            name_Y='碰撞频率';
            handle_fig=figure;
            if flag_multi_Te
                for no_group=1:num_X2
                    loglog(X1,venp(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    hold on
                    loglog(X1,veniz(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                    loglog(X1,veip(:,no_group),[group_style_list{3} group_color_list{no_group}],'LineWidth',plot_line_width);
                end
                line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',1*plot_line_width,'color','k');
                
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('\it{\bf\nu} \rm[Hz]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                %             title(['碰撞频率 \rmat \rm' now_str]);
                grid on%显示网格
                axis([X1(1),X1(end),1e4,2e8])
                
                % 自定义legend
                legend_style_group = zeros(1,num_X2);
                legend_text_group=cell(1,num_X2);
                for no_group=1:num_X2
                    legend_style_group(no_group) = plot(NaN,NaN,['s' group_color_list{no_group}],'MarkerSize',marker_size,'MarkerEdgeColor',group_color_list{no_group},'MarkerFaceColor',group_color_list{no_group});
                    legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
                end
                L1=legend(legend_style_group, legend_text_group{:});
                set(L1,'FontSize',font_size);
                set(L1,'location','north');
                set(L1,'box','off')
                set(L1,'Orientation','horizontal');
                set(L1,'AutoUpdate','off')
                text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
                
                legend_style_group = zeros(1,3);
                legend_text_group={'\it{\bf\nu}\rm_{en}^{(p)}','\it{\bf\nu}\rm_{en}^{(iz)}','\it{\bf\nu}\rm_{ei}^{(p)}'};
                for no_group=1:3
                    legend_style_group(no_group) = plot(NaN,NaN,[group_style_list{no_group} 'k'],'LineWidth',plot_line_width);
                end
                axes2 = axes('position',get(gca,'position'),'visible','off');
                L2=legend(axes2,legend_style_group, legend_text_group{:});
                set(L2,'FontSize',font_size);
                set(L2,'location','north');
                set(L2,'Orientation','horizontal');
                set(L2,'box','off')
            else
                loglog(X1,venp(:,no_X2),'-r','LineWidth',plot_line_width);
                %                 axis equal
                axis([X1(1),X1(end),1e4,5e7]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
                hold on
                loglog(X1,veniz(:,no_X2),'--r','LineWidth',plot_line_width);
                loglog(X1,veip(:,no_X2),'-.b','LineWidth',plot_line_width);
                line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',plot_line_width,'color','k');
                L1=legend('\it{\bf\nu}\rm_{en}^{(p)}','\it{\bf\nu}\rm_{en}^{(iz)}','\it{\bf\nu}\rm_{ei}^{(p)}');
                text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
                set(L1,'FontSize',font_size);
                set(L1,'location','southeast');
                set(L1,'Orientation','horizontal');
                set(L1,'box','off');
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('\it{\bf\nu} \rm[Hz]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                
                %         set(L1,'box','off')
                %             title(['碰撞频率 \rmat \rm' now_str]);
                grid on%显示网格
                text(0.4*X1(1),0.2e4,'(b)','FontSize',font_size)
            end            
            %计算线性性
            fprintf('%s = %.2e \n','15eV,k= ',(log(veip(end,no_X2))-log(veip(1,no_X2)))/(log(X1(end))-log(X1(1))));
            fprintf('%s = %.2e \n','15eV,k= ',(log(veip(11,no_X2))-log(veip(1,no_X2)))/(log(X1(11))-log(X1(1))));
            fprintf('%s = %.2e \n','15eV,k= ',(log(veip(21,no_X2))-log(veip(11,no_X2)))/(log(X1(21))-log(X1(11))));
            fprintf('%s = %.2e \n','15eV,k= ',(log(veip(31,no_X2))-log(veip(21,no_X2)))/(log(X1(31))-log(X1(21))));
            
            % 碰撞截面随能量变化
%             [E_enp,a_enp]=textread('ELASTIC.txt','%n %n','headerlines',10);
%             [E_eniz,a_eniz]=textread('IONIZATION.txt','%n %n','headerlines',10);
%             handle_fig=figure;
%             loglog(E_enp,a_enp,'-r','LineWidth',plot_line_width)
%             hold on
%             loglog(E_eniz,a_eniz,'-.b','LineWidth',plot_line_width)
%             % line([X1(1),X1(end)],[1,1],'linestyle',':','linewidth',1*plot_line_width,'color','k');
%             
%             ylabel('Cross section [m^2]');
%             xlabel('Energy [eV]')
%             set(gca,'FontSize',font_size)
%             set(gca, 'LineWidth',gca_line_width)
%             %             title([name_Y ' \rmat \rm' now_str]);
%             grid on%显示网格
%             % text(0.4*X1(1),0.5e5,'(b)','FontSize',font_size)
%             
%             L1=legend('{\it\bf\sigma}_{en}^{(p)}','{\it\bf\sigma}_{en}^{(iz)}');
%             set(L1,'FontSize',font_size);
%             set(L1,'location','northwest');
%             set(L1,'box','off')
%             set(L1,'AutoUpdate','off')
            
            % 特征尺寸
            handle_fig=figure;
            if flag_multi_Te
                for no_group=1:num_X2
                    loglog(X1,skin_depth_eff(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    hold on
                    loglog(X1,delta_st(:,no_group),[group_style_list{3} group_color_list{no_group}],'LineWidth',plot_line_width);
                    loglog(X1,wavelength_wave(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                end
                line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle',':','linewidth',0.5*plot_line_width,'color','k');
                
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('Characteristic length [m]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                %             title(['特征尺寸 \rmat \rm' now_str]);
                grid on%显示网格
                axis([X1(1),X1(end),5e-3,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
                
                % 自定义legend
                legend_style_group = zeros(1,num_X2);
                legend_text_group=cell(1,num_X2);
                for no_group=1:num_X2
                    legend_style_group(no_group) = plot(NaN,NaN,['s' group_color_list{no_group}],'MarkerSize',marker_size,'MarkerEdgeColor',group_color_list{no_group},'MarkerFaceColor',group_color_list{no_group});
                    legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
                end
                L1=legend(legend_style_group, legend_text_group{:});
                set(L1,'FontSize',font_size);
                set(L1,'location','north');
                set(L1,'box','off')
                set(L1,'Orientation','horizontal');
                set(L1,'AutoUpdate','off')
                text(1.1*X1(1),1.5*r_plasma,'\itr\rm_{plasma}','FontSize',font_size)
                
                legend_style_group = zeros(1,3);
                legend_text_group={'\it{\bf\delta}_{\rmeff}','\it{\bf\lambda}','\it{\bf\delta}_{\rmst}'};
                for no_group=1:3
                    legend_style_group(no_group) = plot(NaN,NaN,[group_style_list{no_group} 'k'],'LineWidth',plot_line_width);
                end
                axes2 = axes('position',get(gca,'position'),'visible','off');
                L2=legend(axes2,legend_style_group, legend_text_group{:});
                set(L2,'FontSize',font_size);
                set(L2,'location','north');
                set(L2,'Orientation','horizontal');
                set(L2,'box','off')
            else
                loglog(X1,skin_depth_eff(:,no_X2),'-r','LineWidth',plot_line_width);
                hold on
                loglog(X1,delta_st(:,no_X2),'--r','LineWidth',plot_line_width);
                loglog(X1,delta_m(:,no_X2),'-.r','LineWidth',plot_line_width);
                loglog(X1,wavelength_wave(:,no_X2),'-ob','LineWidth',plot_line_width);
                line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle',':','linewidth',plot_line_width,'color','k');
                %                     line([X1(1),5e17],[0.05,0.05],'linestyle','-.','linewidth',1.5,'color','k');
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('Characteristic length [m]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                %                 L1=legend('\it{\bf\delta}_{\rmeff}','\it{\bf\delta}_{\rmst}','\it{\bf\lambda}','\itr\rm_{plasma}');
%                 L1=legend('\it{\bf\delta}_{\rmeff}','\it{\bf\delta}_{\rmst}','\it{\bf\lambda}');
                L1=legend('\it{\bf\delta}_{\rmeff}','\it{\bf\delta}_{\rmst}','\it{\bf\delta}_{\rmm}','\it{\bf\lambda}');
                text(1.1*X1(1),1.5*r_plasma,'\itr\rm_{plasma}','FontSize',font_size)
                set(L1,'FontSize',font_size);
                set(L1,'location','northeast');
                %                 set(L1,'Orientation','horizontal');
                set(L1,'box','off')
                %             title(['特征尺寸 \rmat \rm' now_str]);
                grid on%显示网格
                axis([X1(1),X1(end),5e-3,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            end
            
            % 特征频率
            name_Y='特征频率';
            handle_fig=figure;
            loglog(X1,wpi(:,no_X2),'-r','LineWidth',plot_line_width);
            %             axis equal
            %             axis([X1(1),X1(end),1e6,1e12]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            hold on
            loglog(X1,wpe(:,no_X2),'--b','LineWidth',plot_line_width);
            xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
            ylabel('\it{\bf\omega} \rm[Hz]');
            set(gca,'FontSize',font_size)
            set(gca, 'LineWidth',gca_line_width)
            line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',plot_line_width,'color','k');
            %             L1=legend('\it{\bf\omega}_{\rmpi}','\it{\bf\omega}_{\rmpe}','{\it\bf\omega}');
            L1=legend('\it{\bf\omega}_{\rmpi}','\it{\bf\omega}_{\rmpe}');
            text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
            set(L1,'FontSize',font_size);
            set(L1,'location','northwest');
            set(L1,'box','off')
            %             title(['特征频率 \rmat \rm' now_str]);
            grid on%显示网格
            
            % vst简化表达式与Cazzador fit表达式对比
            switch flag_vst_expression
                case {'Vahedi-simplify','Cazzador-simplify'}
                    name_Y='对比st表达式';
                    handle_fig=figure;
                    set(handle_fig,'defaultAxesColorOrder',[left_color; right_color]);
                    % marker会遮掩线型，而且也不好看
                    flag_logy=true;
                    %             flag_logy=false;
                    if flag_logy
                        yyaxis left;
                        h1=loglog(X1,vst_fit(:,no_X2),'-dr','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','r',...
                            'LineWidth',plot_line_width);
                        yyaxis right;
                        h2=loglog(X1,delta_st_fit(:,no_X2),'-b','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','r',...
                            'LineWidth',plot_line_width);
                        hold on
                        yyaxis left;
                        h3=loglog(X1,vst(:,no_X2),'-sr','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','r',...
                            'LineWidth',plot_line_width);
                        line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',plot_line_width,'color','r');
                        axis([X1(1),X1(end),4e6,1e8])
                        ylabel('{\it\bf\nu}_{st} [Hz]');
                        yyaxis right;
                        h4=loglog(X1,delta_st(:,no_X2),'-ob','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','b',...
                            'LineWidth',plot_line_width);
                        line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle',':','linewidth',plot_line_width,'color','b');
                        axis([X1(1),X1(end),1e-3,3e-1])
                        ylabel('{\it\bf\delta}_{st} [m]');
                    else
                    end
                    xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                    set(gca,'FontSize',font_size)
                    set(gca, 'LineWidth',gca_line_width)
                    grid on%显示网格
                    L1=legend([h1,h3],'{\it\bf\nu}_{st}, numerical','{\it\bf\nu}_{st}, piecewise');
                    set(L1,'FontSize',font_size);
                    set(L1,'location','southwest');
                    set(L1,'Orientation','horizontal');
                    set(L1,'box','off')
                    
                    yyaxis left
                    text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
                    yyaxis right
                     text(1.05*X1(end),1.2*r_plasma,'\itr\rm_{plasma}','FontSize',font_size)
                    
                    axes2 = axes('position',get(gca,'position'),'visible','off');
                    L2=legend(axes2,[h2,h4], '{\it\bf\delta}_{st}, numerical','{\it\bf\delta}_{st}, piecewise');
                    set(L2,'FontSize',font_size);
                    set(L2,'location','northeast');
                    set(L2,'box','off')
                    set(L2,'Orientation','horizontal');
                   
                    
                    no_X1=11;
                    (vst(no_X1)-vst_fit(no_X1))/(vst_fit(no_X1))
                    
                    
            end
            
            return
            
            % 部分图需要手动调整图形大小以避免legend覆盖曲线，因此不自动保存
            %             save_path='d:\School\DoctorProgram\eP-项目笔记文件夹\eP-190821-01激励器FEM模型\200221期刊论文\';
            %             saveas(gcf,[save_path name_Y '.svg'],'svg')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output for paper
        else
            if flag_parameters_coupling
                %         % benchmark with 2019Raunera fig7a
                %         vm_ref=[5.12E+06,7.71E+06,1.52E+07,4.45E+07,7.27E+07,1.37E+08];
                %         vst_ref=[7.77E+06,9.19E+06,9.58E+06,9.58E+06,9.17E+06,8.27E+06];
                %                 figure
                %         loglog(X1,vm,'-r','LineWidth',1.5);
                %         hold on
                %         loglog(X1,vst,'-b','LineWidth',1.5);
                %         hold on
                %         loglog(X1,vm_ref,'--dr','LineWidth',1.5);
                %         hold on
                %         loglog(X1,vst_ref,'--db','LineWidth',1.5);
                %         xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
                %         ylabel('\it\nu\rm[Hz]');
                %         L1=legend('\it\nu\rm_{m}','\it\nu\rm_{st}','\it\nu\rm_{m}-Rauner','\it\nu\rm_{st}-Rauner');
                %         set(L1,'FontSize',10);
                %         set(L1,'location','southeast');
                %         %         set(L1,'box','off')
                %         title(['碰撞频率 \rmat \rm' now_str]);
                %         %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
                %         grid on%显示网格
                %         pause
                
                % 等离子体等效电磁参数
                figure
                yyaxis left
                semilogx(X1,sigma_medium,'-r','LineWidth',1.5);
                ylabel('\it\sigma_{\rmeff}');
                yyaxis right
                semilogx(X1,-eps_medium/eps0,'-b','LineWidth',1.5);
                ylabel('\it-\epsilon_{\rmr-eff}');
                xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
                title(['等离子体等效电磁参数 \rmat \rm' now_str]);
                L1=legend('\it\sigma_{\rmeff}','\it-\epsilon_{\rmr-eff}');
                set(L1,'FontSize',10);
                %                     set(L1,'location','southeast');
                %         set(L1,'box','off')
                %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
                grid on%显示网格
                %             fig_background_transparent(gcf,gca)  %临时使用
                
                figure
                semilogx(X1,sigma_medium./(-w_RF*eps_medium),'-r','LineWidth',1.5);
                ylabel('\it\sigma_{\rmeff}/\it\omega|\epsilon_{\rmeff}|');
                xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
                title(['\it\sigma_{\rmeff}/\it\omega|\epsilon_{\rmeff}|' ' \rmat \rm' now_str]);
                line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
                line([X1(1),X1(end)],[10,10],'linestyle','-.','linewidth',1.5,'color','k');
                %                     L1=legend('\it\sigma_{\rmeff}','\it-\epsilon_{\rmr-eff}');
                %                     set(L1,'FontSize',10);
                %                     set(L1,'location','southeast');
                %         set(L1,'box','off')
                %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
                grid on%显示网格
                pause
            else
                disp('模型参数')
                fprintf('%s = %.2e ,','f',f);
                fprintf('%s = %.2e , ',X_var{idx_X3},X3);
                fprintf('%s = %.2e , ','Tg',Tg)
                fprintf('%s = %.2e~%.2e , ',X_var{idx_X1},X1(1),X1(end))
                fprintf('%s = %.2e~%.2e\n',X_var{idx_X2},X2(1),X2(end))
                
                % TODO：除Te外，全部为SI单位，待添加进图
                %                     name_Y='等离子体等效电磁参数';
                %                     Y1=sigma_medium;
                %                     name_Y1='\it\sigma_{\rmeff}';
                %                     Y2=-w_RF*eps_medium;
                %                     name_Y2='\it-\omega\epsilon_{\rmeff}';
                %                     handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                %                         X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                %                     % 200529 试过logy的图，不怎么好看，只是看上去近了一点
                %                     L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                %                         [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
                %                     set(L1,'location','northwest');
                
                name_Y='等离子体等效电磁参数';
                Y1=sigma_medium;
                name_Y1='\it\sigma_{\rmeff}';
                Y2=-eps_medium/eps0;
                name_Y2='\it-\epsilon_{\rmr-eff}';
                handle_fig=plot_parametric_2Ylogaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
                set(L1,'location','northwest');
                %         ne_center=5e18;
                %         ne_skin=ne_center*0.55;
                %         ne_axi_ave=ne_center*0.42;
                %         line([,X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
                
                Y1=sigma_medium./(-w_RF*eps_medium);
                name_Y1='\it\sigma_{\rmeff}/\it\omega|\epsilon_{\rmeff}|';
                handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                title([name_Y1 ' \rmat \rm' now_str]);
                line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
                line([X1(1),X1(end)],[10,10],'linestyle','-.','linewidth',1.5,'color','k');
                L1=legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
                set(L1,'location','northwest');
                
                pause
                
                %                 name_Y='有损介质参数';
                %                 Y1=-epsp_r_real;
                %                 name_Y1='\it-\epsilon_{\rmr}';
                %                 Y2=-tandelta;
                %                 name_Y2='\rm-tan\it\delta';
                %                 handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                %                     X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                %                 L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                %             [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
                %                 set(L1,'location','northwest');
                %                 pause
                
                %         name_Y='复电导率\it\sigma_{\rmp}';
                %         Y1=sigmap_real;
                %         name_Y1='\rmRe(\it\sigma_{\rmp}\rm)[S/m]';
                %         Y2=-sigmap_imag;
                %         name_Y2='-\rmIm(\it\sigma_{\rmp}\rm) = - \it\omega\epsilon[S/m]';
                %         handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                %             X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                %         name_Y2='-\rmIm(\it\sigma_{\rmp}\rm)[S/m]';
                %         handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                %             X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                %         pause
                
                Y1=skin_depth_eff;
                name_Y1='\it\delta_{\rmeff}\rm[m]';
                handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                title([name_Y1 ' \rmat \rm' now_str ', if \it\delta_{\rmeff}>\itr_{\rmchamber},电磁波穿透等离子体']);
                line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle','-.','linewidth',1.5,'color','k');
                legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end},'\itr_{\rmchamber}');
                
                Y1=wavelength_wave;
                name_Y1='\it\lambda_{\rmwave}\rm[m]';
                handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                title([name_Y1 ' \rmat \rm' now_str ', if \it\lambda_{\rmwave}<\itr_{\rmchamber},各处不同相位']);
                line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle','-.','linewidth',1.5,'color','k');
                legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end},'\itr_{\rmchamber}');
                pause
                
                name_Y='带电粒子对电磁场的响应';
                Y1=wpi_per_w;
                name_Y1='\it\omega_{\rmpi}/\it\omega_{\rmRF}';
                Y2=wpe_per_w;
                name_Y2='\it\omega_{\rmpe}/\it\omega_{\rmRF}';
                handle_fig=plot_parametric_2Ylog(name_Y,Y1,name_Y1,Y2,name_Y2,...
                    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
                L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                    '\it\omega_{\rmpi}/\it\omega_{\rmRF}\rm=1');
                set(L1,'location','northwest');
                pause
                
                name_Y='碰撞频率';
                Y1=vm;
                name_Y1='\nu_{\rmm}\rm[Hz]'; %可能是rad/s
                Y2=vst;
                name_Y2='\nu_{\rmst}\rm[Hz]';
                handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                line([X1(1),X1(end)],[w_RF,w_RF],'linestyle','-.','linewidth',1.5,'color','k');
                L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                    '\it\omega_{\rmRF}');
                set(L1,'location','northwest');
                set(L1,'box','off')
                pause
                
                Y1=veff_per_w;
                name_Y1='\it\nu_{\rmeff}/\it\omega_{\rmRF}';
                handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                title([name_Y1 ' \rmat \rm' now_str ', if>>1, \sigma_{\rmDC}可用']);
                line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
                line([X1(1),X1(end)],[10,10],'linestyle','-.','linewidth',1.5,'color','k');
                legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
                pause
                
                % 欧姆加热不同种类碰撞频率-for 聚变负源典型参数
                figure
                loglog(X1,venp(:,1),'-r','LineWidth',1.5);
                hold on
                loglog(X1,veniz(:,1),'--r','LineWidth',1.5);
                hold on
                loglog(X1,veip(:,1),'-.r','LineWidth',1.5);
                hold on
                loglog(X1,veip(:,no_mid_X2),'-.b','LineWidth',1.5);
                hold on
                loglog(X1,veip(:,end),'-.g','LineWidth',1.5);
                xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
                ylabel('碰撞频率');
                L1=legend('venp','veniz',...
                    ['veip at ' legend_X2{1}],['veip at ' legend_X2{no_mid_X2}],['veip at ' legend_X2{end}]);
                set(L1,'FontSize',10);
                title(['不同种类碰撞频率 \rmat \rm' now_str]);
                axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
                grid on%显示网格
                
                %                         % 欧姆加热不同种类碰撞频率-for 2014Cazzador
                %                         figure
                %                         loglog(X1,venp(:,1),'-r','LineWidth',1.5);
                %                         hold on
                %                         loglog(X1,veniz(:,1),'--r','LineWidth',1.5);
                %                         hold on
                %                         loglog(X1,veip(:,1),'-.r','LineWidth',1.5);
                %                         xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
                %                         ylabel('碰撞频率');
                %                         L1=legend('venp','veniz',['veip at ' legend_X2{1}]);
                %                         set(L1,'FontSize',10);
                %                         title(['不同种类碰撞频率 \rmat \rm' now_str]);
                %                         axis([X1(1),X1(end),1e6,1e8]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
                %                         grid on%显示网格
                
                % 测试stoc模型
                % stoc模型参数
                name_Y='stoc模型参数';
                Y1=alpha_st;
                name_Y1='\rm\alpha_{\rmst}';
                Y2=delta_st;
                name_Y2='\rm\delta_{\rmst}';
                handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                yyaxis left;
                line([X1(1),X1(end)],[0.03,0.03],'linestyle','-.','linewidth',1.5,'color','k');
                L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                    '\it\alpha\rm=0.03',...
                    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
                %         set(L1,'location','northwest');
                
                % vst简化表达式与Cazzador fit表达式对比
                switch flag_vst_expression
                    case {'Vahedi-simplify','Cazzador-simplify'}
                        name_Y='\it\nu_{\rmst}';
                        Y1=vst_fit;
                        name_Y1='Cazzador fit';
                        Y2=vst;
                        name_Y2=flag_vst_expression;
                        handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                            X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                        line([X1(1),X1(end)],[w_RF,w_RF],'linestyle','-.','linewidth',1.5,'color','k');
                        L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                            [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                            '\it\omega_{\rmRF}');
                        set(L1,'location','northwest');
                        
                        name_Y='\it\delta_{\rmst}';
                        Y1=delta_st_fit;
                        name_Y1='Cazzador-fit';
                        Y2=delta_st;
                        name_Y2=flag_vst_expression;
                        handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                            X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                        %                     line([X1(1),X1(end)],[w_RF,w_RF],'linestyle','-.','linewidth',1.5,'color','k');
                        %                     L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                        %                         [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                        %                         '\it\omega_{\rmRF}');
                        %                     set(L1,'location','northwest');
                end
                pause
            end
        end
    end
end

% error('end')

%% 螺线管线圈ICP源的电模型
% 改成函数，变量传递是可以接受的。但作为一个实验室代码，可能这样更便于阅读――做好版本管理就可以了
switch flag_electric_model
    case 'transformer_base'
        % 变压器模型
        disp('## ICP源电模型：变压器模型')
        if l_plasma<r_plasma
            warn('plasma为短螺线管')
        end
        % fprintf('%s = %.2e \n','IC的长冈系数',check_Nagaoka(2*r_plasma/l_plasma));
        Lmp=mu0*pi*r_plasma^2*check_Nagaoka(2*r_plasma/l_plasma)/l_plasma;  %等离子体感应电感-螺线管电感经验公式。相较于ZP，考虑了长冈系数
        % Lmp1=mu0*pi*r_plasma^2*(1-8*r_plasma/(3*pi*l_plasma))/l_plasma; % 某论文
        % (1-8*r_plasma/(3*pi*l_plasma))与长冈系数结果有显著差别。长冈系数计算结果更接近2018Jain
        for X1i=1:num_X1%不同的ne
            for X2i=1:num_X2%不同的Te
                Rp(X1i,X2i)=pi*r_plasma/(skin_depth_eff(X1i,X2i)*sigmap_real(X1i,X2i)*l_plasma); %二次侧等离子体电阻，2011Chabert
                Lp(X1i,X2i)=Rp(X1i,X2i)/veff(X1i,X2i); %等离子体中电子惯性形成的电感。
                if flag_good_conductor_approximation
                    % 这似乎并不是忽略复电导率虚部
                    disp('变压器模型中使用良导体近似')
                    Rp(X1i,X2i)=2*pi*r_plasma/(skin_depth_eff(X1i,X2i)*sigmap_real(X1i,X2i)*l_plasma); %二次侧等离子体电阻，2011Chabert
                    Lp(X1i,X2i)=0;
                end
                
                M=sqrt(Lcoil*(Lmp+Lp(X1i,X2i)))*r_plasma^2/r_coil^2; %线圈和等离子体之间的互感
                %         M=(N-0)*Lmp;          %by ZP,但该表达式未考虑线圈半径
                %         fprintf('%s = %.2e , %s = %.2e , %s = %.2e\n','1995Vehadi,Rp',Rp,'Lmp',Lmp,'M',M);
                
                %         %2018Jainb中
                %         if flag_single_independent_variable
                %             disp('使用2018Jainb的Rp、Lmp、M表达式')
                %             Rp(X1i,X2i)=2*pi*r_plasma/(skin_depth_eff(X1i,X2i)*sigmaeff_dc(X1i,X2i)*l_plasma);
                %             Lmp=0.002*pi*(2*r_plasma*100)*(log(4*2*r_plasma/l_plasma)-0.5)*1e-6; %该表达式结果可能为负值，则bug;即使取绝对值，也超出适用范围
                %             if Lmp<0
                %                 warning('Lmp<0. We will use Lmp=-Lmp.')
                %                 pause
                %                 Lmp=-Lmp;
                %             end
                %             M=0.0095*N*1e-6*(2*r_plasma*100)^2/sqrt((2*r_coil*100)^2+(l_coil*100)^2);
                %             fprintf('%s = %.2e , %s = %.2e , %s = %.2e\n','2018Jainb,Rp',Rp,'Lmp',Lmp,'M',M);
                %         end
                
                PER(X1i,X2i)=Rp(X1i,X2i)*w_RF^2*M^2/(Rp(X1i,X2i)^2+w_RF^2*(Lmp+Lp(X1i,X2i))^2);  %换算到一次侧的等离子体等效电阻
                Rsys(X1i,X2i)=PER(X1i,X2i)+Rmetal(X1i,X2i); %系统等效阻抗的电阻分量
                Lplasma(X1i,X2i)=-M^2*w_RF^2*(Lmp+Lp(X1i,X2i))/(Rp(X1i,X2i)^2+w_RF^2*(Lmp+Lp(X1i,X2i))^2); %换算到一次侧的等离子体等效电感
                if ~flag_parametric_analysis && Lplasma(X1i,X2i)>0
                    warning('Lplasma应<0, please stop and check')
                    pause
                end
                Xplasma(X1i,X2i)=w_RF*Lplasma(X1i,X2i);
                Lsys(X1i,X2i)=Lcoil+Lplasma(X1i,X2i);   %系统等效阻抗的电感分量
                Xsys(X1i,X2i)=w_RF*Lsys(X1i,X2i); %系统等效阻抗的电抗分量
                P_abs(X1i,X2i)=PER(X1i,X2i)*Icoil_rms(X1i,X2i)^2;            %吸收功率 即Pcp_real
            end
        end
    case 'analytical_base'
        % 解析模型 2011Chabert
        disp('## ICP源电模型：轴对称平行平面解析模型')
        % 几何校正
        l_equ=l_chamber; %等离子体、线圈均与腔室同长
        %         N_equ=N_coil;
        N_equ=N_coil*l_equ/l_coil; %变换线圈长度
        r_p=r_chamber; %等离子体与腔室同半径
        
        if flag_parametric_analysis
            
        else
            % 增大线圈阻抗的校正方法
            %             % 理论计算线圈阻抗
            %             fprintf('几何校正前线圈阻抗：')
            %             fprintf('%s = %.2e , ','Rcoil',Rcoil);
            %             fprintf('%s = %.2e \n','Lcoil',Lcoil);
            %             Rcoil_modified=N_equ*2*pi*r_coil*sqrt(mu0*w_RF*sigma_Cu/2)/2/pi/r_wire_of_coil/sigma_Cu;%已考虑集肤效应，未考虑邻近效应
            %             Lcoil_modified=mu0*check_Nagaoka(2*r_coil/l_equ)*pi*r_coil^2*N_equ^2/l_equ;  %一次螺线管线圈电感理论计算表达式
            %             fprintf('几何校正后，理论计算线圈阻抗：')
            %             fprintf('%s = %.2e , ','Rcoil',Rcoil_modified);
            %             fprintf('%s = %.2e \n','Lcoil',Lcoil_modified);
            %             if Rcoil<Rcoil_modified
            %                 Rcoil=Rcoil_modified*ones(size(Rcoil));
            %                 fprintf('使用几何校正后')
            %             else
            %                 fprintf('使用几何校正前')
            %             end
            %             fprintf('%s = %.2e \n','Rcoil',Rcoil);
            %             if abs(Lcoil)<abs(Lcoil_modified)
            %                 Lcoil=Lcoil_modified;
            %                 fprintf('使用几何校正后')
            %             else
            %                 fprintf('使用几何校正前')
            %             end
            %             fprintf('%s = %.2e \n','Lcoil',Lcoil);
            %             Qcoil=w_RF*Lcoil/Rcoil;
        end
        
        % 介质窗相对介电常数，玻璃约2，氧化铝陶瓷约9
        epst_r=9;
        k_t_wave=w_RF*sqrt(mu0*eps0*epst_r);
        if k_t_wave*r_p>1
            warn('介质管中磁场不是常数，该模型不适用')
        end
        % 重构：将以上参数代码提到最前面输入部分去
        Hz0=N_equ*Icoil_rms(X1i,X2i)/l_equ;
        for X1i=1:num_X1%不同的ne
            for X2i=1:num_X2%不同的Te
                Hz_p=@(r)Hz0*besselj(0,k_p_wave(X1i,X2i)*r)/besselj(0,k_p_wave(X1i,X2i)*r_p);
                % 与2018Zhao中带FS的FEM模型的20A时磁场结果相比，小了
                Etheta_p=@(r)-1i*k_p_wave(X1i,X2i)*Hz0*besselj(1,k_p_wave(X1i,X2i)*r)/...
                    (w_RF*eps0*epsp_r_real(X1i,X2i)+besselj(0,k_p_wave(X1i,X2i)*r_p));
                % 与2018Zhao中带FS的FEM模型的20A时电场结果相比，大了
                Hz_t=@(r)Hz0; %介质管中磁场均匀
                
                %检验表达式用-断点调试
                % 复坡印廷定理计算复功率
                Etheta_t_at_rc=Etheta_p(r_p)*r_p/r_coil-1i*w_RF*mu0*Hz0*(r_coil^2-r_p^2)/(2*r_coil); %rc处电场
                Pcp1=-pi*r_coil*l_equ*Etheta_t_at_rc*Hz_t(r_coil);
                % ERROR 该表达式与简化表达式结果不一致
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %简化表达式计算复功率
                Pcp_part1=k_p_wave(X1i,X2i)*r_p*besselj(1,k_p_wave(X1i,X2i)*r_p)/...
                    (w_RF*eps0*epsp_r(X1i,X2i)*besselj(0,k_p_wave(X1i,X2i)*r_p));
                Pcp_part2=Pcp_part1+w_RF*mu0*(r_coil^2-r_p^2)/2;
                Pcp=1i*pi*N_equ^2*Icoil_rms(X1i,X2i)^2*Pcp_part2/l_equ;
                
                Pcp_real=real(Pcp); %即Pabs
                Pabs(X1i,X2i)=Pcp_real;
                %                 %检验表达式用-断点调试
                %                 Pcp_real1_part1=1i*k_p_wave(X1i,X2i)*r_p*besselj(1,k_p_wave(X1i,X2i)*r_p)/(w_RF*eps0*epsp_r(X1i,X2i)*besselj(0,k_p_wave(X1i,X2i)*r_p));
                %                 Pcp_real1=pi*N_equ^2*I_coil^2*real(Pcp_real1_part1)/l_equ; %与Pcp_real一致，说明上下表达式自洽
                PER(X1i,X2i)=2*Pcp_real/Icoil_rms(X1i,X2i)^2; %换算到一次侧的等离子体等效电阻
                I_p=l_equ*Hz0*(1/besselj(0,k_p_wave(X1i,X2i)*r_p)-1); %等离子体中电流
                Rp(X1i,X2i)=2*Pcp_real/abs(I_p)^2; %二次侧等离子体电阻
                %                 %检验表达式用-断点调试
                % 与LZS实现结果一致
                %                 Rp1_part1=1i*k_p_wave(X1i,X2i)*r_p*besselj(1,k_p_wave(X1i,X2i)*r_p)/(w_RF*eps0*epsp_r(X1i,X2i)*besselj(0,k_p_wave(X1i,X2i)*r_p));
                % 可能少了一个N^2
                %                 Rp1=2*pi*real(Rp1_part1)/l_equ/(1/besselj(0,k_p_wave(X1i,X2i)*r_p)-1)^2; %与Rp一致，说明上下表达式自洽
                Rsys(X1i,X2i)=PER(X1i,X2i)+Rmetal(X1i,X2i); %系统等效阻抗的电阻分量
                Xsys(X1i,X2i)=2*imag(Pcp)/Icoil_rms(X1i,X2i)^2; %系统等效阻抗的电抗分量
                Lsys(X1i,X2i)=Xsys(X1i,X2i)/w_RF; %系统等效阻抗的电感分量
                Lplasma(X1i,X2i)=Lsys(X1i,X2i)-Lcoil;  % 一次侧的等离子体电感（定义的）
                if ~flag_parametric_analysis && Lplasma(X1i,X2i)>0
                    warning('Lplasma应<0, please stop and check')
                    pause
                end
                Xplasma(X1i,X2i)=w_RF*Lplasma(X1i,X2i); %一次侧的等离子体导致电抗
            end
        end
    case 'transformer_2011Chabert'
        warning('施工中. Please stop and check.')
    otherwise
        warning('Unexpected electric model. Please stop and check.')
        pause
end
PTE=PER./Rsys;   %射频功率传输效率 PTE
PCF=Rsys./Xsys;  %激励器射频功率耦合因数

%%%%%%%%%%%%%%%%%%%%%%%%%%% 后处理 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_output_electric_model
    if ~flag_parametric_analysis
        disp('等离子体等效阻抗')
        fprintf('%s = %.2e Ω, ','PER',PER);
        fprintf('%s = %.2e , ','Lplasma',Lplasma);
        fprintf('%s = %.2e \n','Xplasma',Xplasma);
        disp('系统等效阻抗Zsys')
        fprintf('%s = %.2e , ','Rsys',Rsys);
        fprintf('%s = %.2e , ','Lsys',Lsys);
        fprintf('%s = %.2e \n','Xsys',Xsys);
        fprintf('%s = %.2e \n','Rmetal',Rmetal);
        fprintf('%s = %.2e \n','射频功率传输效率 PTE',PTE);
        fprintf('%s = %.2e \n','射频功率耦合因数 PCF',PCF);
    else
        if flag_parameters_coupling
            % 等离子体等效阻抗
            figure
            semilogx(X1,PER,'-r','LineWidth',1.5);
            hold on
            semilogx(X1,-Xplasma,'-b','LineWidth',1.5);
            ylabel('\itZ_{\rmeff}');
            xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            title(['等离子体等效阻抗 \rmat \rm' now_str]);
            L1=legend('\itPER','\it-X_{\rmplasma}');
            set(L1,'FontSize',10);
            %                     set(L1,'location','southeast');
            %         set(L1,'box','off')
            %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            grid on%显示网格
            %             fig_background_transparent(gcf,gca)  %临时使用
            
            % 系统等效阻抗
            figure
            semilogx(X1,Rsys,'-r','LineWidth',1.5);
            hold on
            semilogx(X1,Xsys,'-b','LineWidth',1.5);
            ylabel('\itZ_{\rmeff}');
            xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            title(['系统等效阻抗 \rmat \rm' now_str]);
            L1=legend('\itR_{\rmsys}','\itX_{\rmsys}');
            set(L1,'FontSize',10);
            %                     set(L1,'location','southeast');
            %         set(L1,'box','off')
            %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            grid on%显示网格
            %             fig_background_transparent(gcf,gca)  %临时使用
            
            % PTE与PCF
            figure
            yyaxis left
            semilogx(X1,PTE,'-r','LineWidth',1.5);
            ylabel('\itPTE');
            yyaxis right
            semilogx(X1,PCF,'-b','LineWidth',1.5);
            ylabel('\itPCF');
            xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            title(['PTE&PCF \rmat \rm' now_str]);
            L1=legend('\itPTE','\itPCF');
            set(L1,'FontSize',10);
            %                     set(L1,'location','southeast');
            %         set(L1,'box','off')
            %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            grid on%显示网格
            %             fig_background_transparent(gcf,gca)  %临时使用
            
        else
            name_Y='等离子体等效阻抗\itZ_{\rmplasma}';
            Y1=PER;
            name_Y1='\itPER\rm, namely \itR_{\rmplasma}\rm[\Omega]';
            Y2=Xplasma;
            name_Y2='\itX_{\rmplasma}\rm[\Omega]';
            handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            line([X1(1),X1(end)],[0,0],'linestyle','-.','linewidth',1.5,'color','k');
            L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                'Z=0');
            set(L1,'location','southwest');
            
            Y2=-Lplasma*1e6;
            name_Y2='-\itL_{\rmplasma}\rm[\muH]';
            handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            pause
            
            name_Y='系统等效阻抗\itZ_{\rmsys}即实验可测的端口阻抗';
            Y1=Rsys;
            name_Y1='\itR_{\rmsys}\rm[\Omega]';
            Y2=Xsys;
            name_Y2='\itX_{\rmsys}\rm[\Omega]';
            handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            Y2=Lsys*1e6;
            name_Y2='\itL_{\rmsys}\rm[\muH]';
            handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            pause
            

            
            
            name_Y='系统功率传输效率PTE和耦合因数PCF';
            Y1=PTE;
            name_Y1='\itPTE';
            Y2=PCF;
            name_Y2='\itPCF';
            handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
            set(L1,'location','northwest');
            pause
            
            %     figure(2)
            %     [c,h]=contour(ne,Te,PTE,'LevelList',[0.7:0.05:1]);
            %     clabel(c,h,'FontSize',11);
            %     title('Power tranfer efficinecy, at \itT_{\rmgas}\rm=400 K，\itp_{gas}\rm=0.3 Pa');
            %     xlabel('\itn_{\rme} \rm[m^{-3}]');
            %     ylabel('\itT_{\rme} \rm[eV]');
            %     set(gca,'XScale','log');
            %     grid on
            %
            %     figure(3)
            %     plot(Te,PER(11,:),'-g','LineWidth',1.5);
            %     hold on
            %     plot(Te,PER(21,:),'-b','LineWidth',1.5);
            %     hold on
            %     plot(Te,PER(31,:),'-k','LineWidth',1.5);
            %     hold on
            %     grid on
            %     L1=legend('\itn_{e}\rm=1e16 m^{-3}','\itn_{e}\rm=1e17 m^{-3}','\itn_{e}\rm=1e18 m^{-3}');
            %     set(L1,'FontSize',10);
            %     title('Plasma equivalent resitance,at at \itT_{\rmgas}\rm=400 K，\itp_{\rmgas}\rm=0.3 Pa');
            %     xlabel('\itT_{\rme} \rm[eV]');
            %     ylabel('Plasma equivalent resitance,\itR_{\rmeq} \rm[\Omega]');
        end
    end
end

%% aid function
% 使用函数，主要是方便测试和复用。
% 弃用函数，在main中直接写表达式。因为在main中只使用一处，且基本不会复用，而写函数传参麻烦。

%美观目的的微调在origin中进行
function handle_fig=plot_parametric_1Y(Y1,name_Y1,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    semilogx(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,end),'-g','LineWidth',1.5);
else
    plot(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,end),'-g','LineWidth',1.5);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']); %\it斜体 \rm 正体 ^上标 _下标 {}组合
ylabel(name_Y1);
L1=legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
set(L1,'FontSize',10);
title([name_Y1 ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1与Y2在单纵轴下随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    semilogx(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,end),'-g','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,end),'--g','LineWidth',1.5);
else
    plot(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,end),'-g','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,end),'--g','LineWidth',1.5);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
ylabel([name_Y1 ' & ' name_Y2]);
L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
set(L1,'FontSize',10);
title([name_Y ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1与Y2双纵轴随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    yyaxis left;
    semilogx(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    semilogx(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
else
    yyaxis left;
    plot(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    plot(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
set(L1,'FontSize',10);
title([name_Y ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_1Ylog(Y1,name_Y1,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    loglog(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,end),'-g','LineWidth',1.5);
else
    semilogy(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,end),'-g','LineWidth',1.5);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']); %\it斜体 \rm 正体 ^上标 _下标 {}组合
ylabel(name_Y1);
L1=legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
set(L1,'FontSize',10);
title([name_Y1 ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_2Ylog(name_Y,Y1,name_Y1,Y2,name_Y2,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1与Y2在单对数纵轴下随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    loglog(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,end),'-g','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,end),'--g','LineWidth',1.5);
else
    semilogy(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,end),'-g','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,end),'--g','LineWidth',1.5);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
ylabel([name_Y1 ' & ' name_Y2]);
L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
set(L1,'FontSize',10);
title([name_Y ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_2Ylogaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1与Y2双纵轴随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    yyaxis left;
    loglog(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    loglog(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
else
    yyaxis left;
    semilogy(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    semilogy(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
set(L1,'FontSize',10);
title([name_Y ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function fig_background_transparent(gcf,gca)
% 使图片背景透明，用于复制到文档
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
L=legend('show');
set(L,'box','off')
end

% skin_depth(5.8e7,50,1)=0.0093 经测试，skin_depth可用
function skin_d=skin_depth(sigma,f,mu_r)
% calculate skin depth of metal
mu0=4*pi*1e-7;         %真空磁导率
skin_d=sqrt(1/(pi*f*mu_r*mu0*sigma));
end

% lambda(1,1e6,1)=300 经测试，lambda可用
function lambda_d=lambda(eps_r,f,mu_r)
% calculate wave length in low-loss medium or lossless medium
mu0=4*pi*1e-7;         %真空磁导率
eps0=8.85e-12;          %真空介电常数
medium_v=1/sqrt(mu0*mu_r*eps0*abs(eps_r));
lambda_d=medium_v/f;
end