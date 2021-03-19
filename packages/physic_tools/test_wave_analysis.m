%% Main function to generate tests
function tests = test_wave_analysis
% test wave_analysis
tests = functiontests(localfunctions);
end

%% Test Functions
function test_get_skin_depth_plasma(testCase)
% test get_skin_depth_plasma


% skin_depth(5.8e7,50,1)=0.0093 
% lambda(1,1e6,1)=300
end

function test_skin_depth(testCase)
% test different get_skin_depth_plasma

skin_depth=sqrt(1/(pi*f*mu_r*mu0*sigma));

f=1e6*[ones(1,3);10*ones(1,3)];
vc=1e7*[ones(1,3);10*ones(1,3)];
wpe=1e9*[ones(1,3);10*ones(1,3)];
r_plasma=1;

get_plasma_skin_depth('as-medium',1e6,vc,wpe,r_plasma)
get_plasma_skin_depth('as-medium',f,vc,wpe,r_plasma)
get_plasma_skin_depth('as-medium-simplified-finite-radius',f,vc,wpe,r_plasma)

tolerance=1/100;
verify_equal=@(actual, expected) verifyEqual(testCase,actual,expected,'RelTol',tolerance);

% 根据ZP/ZC结果验算
verify_equal(kenp, [1.0505e-13*ones(1,3);9.1450e-14*ones(1,3)])
verify_equal(keniz, [1.3344e-21*ones(1,3);6.4445e-15*ones(1,3)])

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


% 对比不同趋肤深度表达式
if ~flag.sweep&&flag_output_plasma_model
    % 对比不同适用情况的集肤深度表达式
    fprintf('%s = %.2e\n','EM经典δ',delta_fun(plasma.nu_eff,plasma.wpe));
    disp('若wpe >> w,vc(碰撞频率)，集肤深度表达式可简化') %nu_m/plasma.nu_st<veff<<wpe
    fprintf('%s = %.2e , ','wpe',plasma.wpe);
    fprintf('%s = %.2e \n','veff',plasma.nu_eff);
    fprintf('%s = %.2e \n','简化经典δ',delta_simplified_fun(plasma.nu_eff,plasma.wpe));
    %                 fprintf('%s = %.2e \n','平面线圈ICP源考虑几何效应的简化δ',delta_geo_simplified_fun(plasma.nu_eff,plasma.wpe));
    
    % 检验表达式：相同参数下，delta_fun结果应与后文skin_depth_general_medium一致
    %             % 使用良导体参数，验证一般介质集肤深度公式与良导体集肤深度公式，在良导体参数下一致
    %                             sigmap_real=constants.sigma_Cu;
    %                             f=50;
    %                             w_RF=2*pi*f;
    %                             epsp_r_real=1;
    %                             epsp_r=epsp_r_real-1i*constants.sigma_Cu/w_RF/constants.eps0;
    %                             tandelta=constants.sigma_Cu/w_RF/constants.eps0;
    
    skin_depth_good_conductor=skin_depth(plasma.sigmap_real,f,1); %良导体集肤深度
    k_p_wave=plasma.w_RF*sqrt(plasma.mup*constants.eps0*plasma.epsp_r); %wave number
    alpha_wave=-imag(k_p_wave); %对于良导体，衰减常数≈相位常数
    beta_wave=real(k_p_wave); %对于良导体，衰减常数≈相位常数
    %                 wavelength_wave=2*pi/beta_wave; %导体趋肤深度等于波长的2π倍
    %             % 检验表达式-断点调试
    %             % RF-ICP一般复介电常数实部为负值，因此表达式与常见表达式略有不同
    %             alpha_wave=plasma.w_RF*sqrt(-plasma.mup*constants.eps0*epsp_r_real*(sqrt(1+tandelta^2)+1)/2); %衰减常数
    skin_depth_general_medium=1/alpha_wave;
    fprintf('%s = %.2e , ','skin depth: 良导体表达式',skin_depth_good_conductor);
    fprintf('%s = %.2e \n','一般介质中表达式',skin_depth_general_medium);
    % 使用良导体参数，结果一致为skin_depth(5.8e7,50,1)=0.0093
    % 使用ELISE等体参数，集肤深度良导体公式结果略大于一般介质公式结果，可忽略差别
    end
end

function test_get_wavelength(testCase)
% test get_wavelength
%             % 对比不同波长表达式
medium_v=1/sqrt(mu0*mu_r*eps0*eps_r);
    lambda_d=medium_v/f;

    lambda(1,1e6,1)
    
%             if ~flag.sweep&&flag_output_plasma_model
%                 %             % 使用无损介质参数，以验证一般介质中波长公式适用于无损介质
%                 %             f=1e6;
%                 %             w_RF=2*pi*f;
%                 %             epsp_r_real=1;
%                 %             epsp_r=epsp_r_real;
%                 %             tandelta=0;
%
%                 wavelength_lossless_medium=lambda(plasma.epsp_r_real,f,1); %无损介质中波长
%                 k_p_wave=w_RF*sqrt(plasma.mup*constants.eps0*plasma.epsp_r); %无损时 波数=相位常数
%                 alpha_wave=-imag(k_p_wave); %无损时衰减常数=0
%                 beta_wave=real(k_p_wave); %无损时 波数=相位常数
% %                 %检验表达式用-断点调试
% %                 % RF-ICP一般复介电常数实部为负值，因此表达式与常见表达式略有不同
% %                 beta_wave=w_RF*sqrt(-plasma.mup*constants.eps0*epsp_r_real*(sqrt(1+tandelta^2)-1)/2); %传播常数
%                 wavelength_general_medium=2*pi/beta_wave; %见前文
%                 fprintf('%s = %.2e , ','wavelength: 无损介质中',wavelength_lossless_medium);
%                 fprintf('%s = %constants.e \n','一般介质中',wavelength_general_medium);
%                 % 使用无损介质参数，结果一致为lambda(1,1e6,1)=300
%                 % 使用ELISE等离子体参数，波长一般介质公式结果约为无损介质公式结果的一半，考虑导电性后分布参数效应更明显
%             end
%             %测试时应在此设置断点

end



%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
% set a new path, for example
end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example
end

%% Optional fresh fixtures  
function setup(testCase)  % do not change function name
% open a figure, for example
end

function teardown(testCase)  % do not change function name
% close figure, for example
end
