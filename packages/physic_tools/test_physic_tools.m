%% Main function to generate tests
function tests = test_physic_tools
% test get_plasma_derived_parameters
tests = functiontests(localfunctions);
end

%% Test Functions
function test_basic(testCase)
% test basic
% the same as Nix_M
constants=get_constants();
tolerance=1/100;
verify_equal=@(actual, expected) verifyEqual(testCase,actual,expected,'RelTol',tolerance);

n0=4e17;
T0=2;
debye_length=get_debye_length( n0, T0 );
omega_pe= get_omega_pe( n0 );
omega_pHp= get_omega_pi( n0, 1 ,1 );
omega_pHn= get_omega_pi( n0,-1,1 );
vth_e=sqrt(-constants.q_m_ration_e*T0);
N_D=get_N_D(n0, T0);

verify_equal(debye_length, 1.66e-5)
verify_equal(omega_pe/2/pi, 5.679e9)
verify_equal(omega_pHp/2/pi, 1.325e8)
verify_equal(omega_pHn/2/pi, 1.325e8)
verify_equal(vth_e, 5.93e5)
verify_equal(N_D, n0*debye_length^3*4*pi/3)
end

function test_vectorization(testCase)
% test 基本函数支持矢量化
% the same as Nix_M
constants=get_constants();
tolerance=1/100;
verify_equal=@(actual, expected) verifyEqual(testCase,actual,expected,'RelTol',tolerance);

n0=4e17*ones(2,3);
T0=2*ones(2,3);
debye_length=get_debye_length( n0, T0 );
omega_pe= get_omega_pe( n0 );
omega_pHp= get_omega_pi( n0, 1 ,1 );
omega_pHn= get_omega_pi( n0,-1,1 );
vth_e=sqrt(-constants.q_m_ration_e*T0);
N_D=get_N_D(n0, T0);

verify_equal(debye_length, 1.66e-5*ones(2,3))
verify_equal(omega_pe/2/pi, 5.679e9*ones(2,3))
verify_equal(omega_pHp/2/pi, 1.325e8*ones(2,3))
verify_equal(omega_pHn/2/pi, 1.325e8*ones(2,3))
verify_equal(vth_e, 5.93e5*ones(2,3))
verify_equal(N_D, n0.*debye_length.^3*4*pi/3)
end

function test_get_Xsec_and_k(testCase)
% test 	get_Xsec and get_k
% Xsec截面的可靠性不在此处test，通过结果ν
% 与2014Cazzador中ν对比来test

Xsec=get_Xsec(true);
close
Te=[ones(1,3);10*ones(1,3)];
kenp=get_k(Te, Xsec.enp);
keniz=get_k(Te, Xsec.eniz);

tolerance=1/100;
verify_equal=@(actual, expected) verifyEqual(testCase,actual,expected,'RelTol',tolerance);

% 根据ZP/ZC结果验算
verify_equal(kenp, [1.0505e-13*ones(1,3);9.1450e-14*ones(1,3)])
verify_equal(keniz, [1.3344e-21*ones(1,3);6.4445e-15*ones(1,3)])
end

function test_get_skin_depth(testCase)
% test 	get_skin_depth and get_skin_depth_plasma

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
