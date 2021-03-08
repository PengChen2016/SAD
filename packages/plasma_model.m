function [ plasma ] = plasma_model( flag, plasma)
% ICP的等效电磁媒质模型 equivalent_medium_model_of_plasma
% 主要基于2014Cazzador、1995Vahedia，借鉴2018Jain

constants=get_constants();

disp('## ICP的等效电磁媒质模型')

% 矢量化
%以下求电子动量转移频率vm
Xsec=get_Xsec(false);
% 电子与中性粒子动量转移碰撞
plasma.kenp=get_k(plasma.Te, Xsec.enp); %速率系数
plasma.nu_enp=plasma.ng.*plasma.kenp; % 碰撞频率
lambda_coll_enp=plasma.ve./plasma.nu_enp; %平均自由程
% 电子与中性粒子电离碰撞
plasma.keniz=get_k(plasma.Te, Xsec.eniz);
plasma.nu_eniz=plasma.ng.*plasma.keniz;
lambda_coll_eniz=plasma.ve./plasma.nu_eniz;
% 电子与离子动量损失碰撞
lambda_e=(constants.eps0*constants.e*plasma.Te.^(3/2))./sqrt(plasma.ne)*...
    (4*pi*constants.mH/(constants.e^3)/(constants.me+constants.mH));%库仑算子中Λ
plasma.nu_eip=plasma.ne.*log(lambda_e)./sqrt((constants.e*plasma.Te).^3*constants.me)*...
    4*sqrt(2*pi)/3*((constants.e^2/4/pi/constants.eps0)^2);
lambda_coll_eip=plasma.ve./plasma.nu_eip;

plasma.nu_m=plasma.nu_enp+plasma.nu_eniz+plasma.nu_eip;
plasma.delta_m=get_plasma_skin_depth('as-medium-simplified',...
    plasma.f,plasma.nu_m,plasma.wpe,plasma.r);

%以下求解随机加热频率vst
if ~isempty(flag.stoc_model)
    plasma=stoc_model(flag.stoc_model, plasma);
end


%             %使用外部的vm、vst数据
%             if flag_single_independent_variable
%                 fprintf('%s = %.2e , ','ourcode，nu_m',nu_m);
%                 fprintf('%s = %.2e \n','vst',vst);
%                 switch flag.input_plasma
%                     case 'ELISE_base'
%                         % 2018Jain中间数据作为变压器输入，用于变压器模型benchmark
%                         disp('使用2018Jain的vm和vst')
%                         nu_m=5.6e6;
%                         vst=4.3e7;
%                     case 'small_source1_LZS'
%                         % 李增山-整体模型耦合解析电模型-2020.03.30\中间数据作为输入，用于解析电模型benchmark
%                         disp('使用李增山-整体模型耦合解析电模型-2020.03.30的vm和vst')
%                         nu_m=1.1238e7;
%                         vst=3.6537e7;
%                 end
%             end



%考虑两种加热机制后
veff(X1i,X2i)=vst(X1i,X2i)+plasma.nu_m;%有效碰撞频率

% 等离子体等效电磁媒质参数
epsp_r(X1i,X2i)=1-plasma.wpe^2/w_RF/(w_RF-1i*veff(X1i,X2i)); %等离子体复相对介电常数,即epsp
% 复介电常数epsc=constants.eps0*epsp_r
sigmap(X1i,X2i)=constants.eps0*plasma.wpe^2/(1i*w_RF+veff(X1i,X2i)); % 复电导率
mu_p=constants.mu0; %等离子体磁导率取为真空磁导率
sigma_medium(X1i,X2i)=constants.eps0*veff(X1i,X2i)*plasma.wpe^2/(w_RF^2+veff(X1i,X2i)^2); %与Re(sigmap)一致
eps_medium(X1i,X2i)=constants.eps0*(1-plasma.wpe^2/(w_RF^2+veff(X1i,X2i)^2)); %与Re(constants.eps0*epsp_r)一致
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
%         epsp_r_real(X1i,X2i)=1-plasma.wpe^2/(veff(X1i,X2i)^2+w_RF^2);         %复介电常数实部
%         epsp_r_imag(X1i,X2i)=-veff(X1i,X2i)*plasma.wpe^2/(veff(X1i,X2i)^2+w_RF^2)/w_RF;         %复介电常数虚部
%         sigma_from_epsp=-w_RF*constants.eps0*epsp_r_imag(X1i,X2i); %与sigmap_real一致
sigmap_real(X1i,X2i)=real(sigmap(X1i,X2i)); %即sigma
sigmap_imag(X1i,X2i)=imag(sigmap(X1i,X2i));
%         %检验表达式用-断点调试
%         sigmap_real(X1i,X2i)=constants.eps0*plasma.wpe^2*veff(X1i,X2i)/(veff(X1i,X2i)^2+w_RF^2);
%         sigmap_imag(X1i,X2i)=constants.eps0*plasma.wpe^2*(-w_RF)/(veff(X1i,X2i)^2+w_RF^2);
%         eps_r_from_sigmap(X1i,X2i)=sigmap_imag(X1i,X2i)/w_RF/constants.eps0+1; %与epsp_r_real一致
%         tandelta_from_sigmap(X1i,X2i)=sigmap_real(X1i,X2i)/(w_RF*constants.eps0*eps_r_from_sigmap(X1i,X2i)); %与tandelta一致

% 电磁波分析
k_p_waplasma.ve=w_RF*sqrt(mu_p*constants.eps0*epsp_r(X1i,X2i)); %wave number
%             gamma_wave=1i*k_p_waplasma.ve; %传播常数

% 等效集肤深度
skin_depth_eff(X1i,X2i)=delta_simplified_fun(veff(X1i,X2i),plasma.wpe);
if ~flag.sweep&&flag_output_plasma_model
    if r_plasma<3*skin_depth_eff(X1i,X2i)
        warning('δ>≈R ，电磁波穿透等离子体，电场基本均匀,集肤深度概念与随机加热模型不适用')
        pause
    end
end
%             alpha_wave=-imag(k_p_waplasma.ve); %衰减常数
%             skin_depth_eff(X1i,X2i)=1/alpha_wave; %一般介质中经典集肤深度
% 对比不同趋肤深度表达式
if ~flag.sweep&&flag_output_plasma_model
    % 对比不同适用情况的集肤深度表达式
    fprintf('%s = %.2e\n','EM经典δ',delta_fun(veff(X1i,X2i),plasma.wpe));
    disp('若wpe >> w,vc(碰撞频率)，集肤深度表达式可简化') %nu_m/vst<veff<<wpe
    fprintf('%s = %.2e , ','wpe',plasma.wpe);
    fprintf('%s = %.2e \n','veff',veff(X1i,X2i));
    fprintf('%s = %.2e \n','简化经典δ',delta_simplified_fun(veff(X1i,X2i),plasma.wpe));
    %                 fprintf('%s = %.2e \n','平面线圈ICP源考虑几何效应的简化δ',delta_geo_simplified_fun(veff(X1i,X2i),plasma.wpe));
    
    % 检验表达式：相同参数下，delta_fun结果应与后文skin_depth_general_medium一致
    %             % 使用良导体参数，验证一般介质集肤深度公式与良导体集肤深度公式，在良导体参数下一致
    %                             sigmap_real=constants.sigma_Cu;
    %                             f=50;
    %                             w_RF=2*pi*f;
    %                             epsp_r_real=1;
    %                             epsp_r=epsp_r_real-1i*constants.sigma_Cu/w_RF/constants.eps0;
    %                             tandelta=constants.sigma_Cu/w_RF/constants.eps0;
    
    skin_depth_good_conductor=skin_depth(sigmap_real(X1i,X2i),f,1); %良导体集肤深度
    k_p_waplasma.ve=w_RF*sqrt(mu_p*constants.eps0*epsp_r(X1i,X2i)); %wave number
    alpha_wave=-imag(k_p_waplasma.ve); %对于良导体，衰减常数≈相位常数
    beta_wave=real(k_p_waplasma.ve); %对于良导体，衰减常数≈相位常数
    %                 wavelength_wave=2*pi/beta_wave; %导体趋肤深度等于波长的2π倍
    %             % 检验表达式-断点调试
    %             % RF-ICP一般复介电常数实部为负值，因此表达式与常见表达式略有不同
    %             alpha_wave=w_RF*sqrt(-mu_p*constants.eps0*epsp_r_real*(sqrt(1+tandelta^2)+1)/2); %衰减常数
    skin_depth_general_medium=1/alpha_wave;
    fprintf('%s = %.2e , ','skin depth: 良导体表达式',skin_depth_good_conductor);
    fprintf('%s = %.2e \n','一般介质中表达式',skin_depth_general_medium);
    % 使用良导体参数，结果一致为skin_depth(5.8e7,50,1)=0.0093
    % 使用ELISE等体参数，集肤深度良导体公式结果略大于一般介质公式结果，可忽略差别
end

%             %使用外部的delta_eff数据
%             if ~flag.sweep&&flag_output_plasma_model
%                 fprintf('%s = %.2e , ','ourcode，veff',veff);
%                 fprintf('%s = %.2e \n','skin_depth_eff',skin_depth_eff);
%                 switch flag.input_plasma
%                     case 'ELISE_base'
%                         % 2018Jain中间数据作为变压器输入，用于变压器模型benchmark
%                         disp('使用2018Jain的δeff')
%                         skin_depth_eff=8.3e-3;
%                 end
%             end

%波长
beta_wave=real(k_p_waplasma.ve); %相位常数
wavelength_waplasma.ve=2*pi/beta_wave; %一般介质中波长
%             % 对比不同波长表达式
%             if ~flag.sweep&&flag_output_plasma_model
%                 %             % 使用无损介质参数，以验证一般介质中波长公式适用于无损介质
%                 %             f=1e6;
%                 %             w_RF=2*pi*f;
%                 %             epsp_r_real=1;
%                 %             epsp_r=epsp_r_real;
%                 %             tandelta=0;
%
%                 wavelength_lossless_medium=lambda(epsp_r_real(X1i,X2i),f,1); %无损介质中波长
%                 k_p_waplasma.ve=w_RF*sqrt(mu_p*constants.eps0*epsp_r(X1i,X2i)); %无损时 波数=相位常数
%                 alpha_wave=-imag(k_p_waplasma.ve); %无损时衰减常数=0
%                 beta_wave=real(k_p_waplasma.ve); %无损时 波数=相位常数
% %                 %检验表达式用-断点调试
% %                 % RF-ICP一般复介电常数实部为负值，因此表达式与常见表达式略有不同
% %                 beta_wave=w_RF*sqrt(-mu_p*constants.eps0*epsp_r_real*(sqrt(1+tandelta^2)-1)/2); %传播常数
%                 wavelength_general_medium=2*pi/beta_wave; %见前文
%                 fprintf('%s = %.2e , ','wavelength: 无损介质中',wavelength_lossless_medium);
%                 fprintf('%s = %constants.e \n','一般介质中',wavelength_general_medium);
%                 % 使用无损介质参数，结果一致为lambda(1,1e6,1)=300
%                 % 使用ELISE等离子体参数，波长一般介质公式结果约为无损介质公式结果的一半，考虑导电性后分布参数效应更明显
%             end
%             %测试时应在此设置断点

if ~flag.sweep&&flag_output_electric_model
    if  wavelength_waplasma.ve<3*r_plasma
        warning('λ＜≈R ，需要考虑相位变化,变压器模型不适用')
        pause
    end
end

%分析等效碰撞频率中的主导加热机制
vm_per_veff(X1i,X2i)=plasma.nu_m/veff(X1i,X2i);

% 分析带电粒子是否响应电磁场
wpe_per_w(X1i,X2i)=plasma.wpe/w_RF;
wpi_per_w(X1i,X2i)=plasma.wpi/w_RF;

%分析是否可以忽略heating model中位移电流，使用复电导率
veff_displacement_ratio(X1i,X2i)=plasma.wpe^2/(w_RF*sqrt(w_RF^2+veff(X1i,X2i)^2));

%分析是否频率非常低到可以使用直流电导率-忽略了复电导率虚部
veff_per_w(X1i,X2i)=veff(X1i,X2i)/w_RF;
sigmaeff_dc(X1i,X2i)=constants.eps0*plasma.wpe^2/veff(X1i,X2i);         %直流电导率



end