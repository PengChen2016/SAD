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
vth_e=sqrt(-constants.q_m_ratio_e*T0);

verify_equal(debye_length, 1.66e-5)
verify_equal(omega_pe/2/pi, 5.679e9)
verify_equal(omega_pHp/2/pi, 1.325e8)
verify_equal(omega_pHn/2/pi, 1.325e8)
verify_equal(vth_e, 5.93e5)
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
vth_e=sqrt(-constants.q_m_ratio_e*T0);

verify_equal(debye_length, 1.66e-5*ones(2,3))
verify_equal(omega_pe/2/pi, 5.679e9*ones(2,3))
verify_equal(omega_pHp/2/pi, 1.325e8*ones(2,3))
verify_equal(omega_pHn/2/pi, 1.325e8*ones(2,3))
verify_equal(vth_e, 5.93e5*ones(2,3))
end

function test_get_Xsec_and_k(testCase)
% test 	formula of get_Xsec and get_k
Xsec=get_Xsec('1990Tawara', true);
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

function test_different_Xsec(testCase)
% test 	get_Xsec and get_k
% Xsec截面的可靠性不在此处test，通过结果ν
% 与2014Cazzador中ν对比来test
% get Xsec
Xsec1=get_Xsec('1990Tawara', false);
Xsec2=get_Xsec('Phelps-m', false);

% plot Xsec
plot_line_width=3;
gca_line_width=1;
font_size=15;

handle_fig=figure;
loglog(Xsec1.enp(:,1),Xsec1.enp(:,2),'-r','LineWidth',2*plot_line_width)
hold on
loglog(Xsec1.eniz(:,1),Xsec1.eniz(:,2),'--k','LineWidth',2*plot_line_width)
loglog(Xsec2.enp(:,1),Xsec2.enp(:,2),'-.m','LineWidth',plot_line_width)
loglog(Xsec2.eniz(:,1),Xsec2.eniz(:,2),'-.c','LineWidth',plot_line_width)

ylabel('Cross section (m^2)');
xlabel('Energy (eV)')
set(gca,'FontSize',font_size)
set(gca, 'LineWidth',gca_line_width)
%             title([name_Y ' \rmat \rm' now_str]);
grid on%显示网格
%     text(0.4*X1(1),0.5e5,'(b)','FontSize',font_size)

L1=legend('{\it\bf\sigma}_{en}^p-1990Tawara','{\it\bf\sigma}_{en}^{iz}-1990Tawara',...
    '{\it\bf\sigma}_{en}^p-Phelps','{\it\bf\sigma}_{en}^{iz}-Phelps');
set(L1,'FontSize',font_size);
set(L1,'location','southwest');
set(L1,'box','off')
set(L1,'AutoUpdate','off')
hold off

% get k
Te=logspace(-1,2,40);
kenp1=get_k(Te, Xsec1.enp);
keniz1=get_k(Te, Xsec1.eniz);
kenp2=get_k(Te, Xsec2.enp);
keniz2=get_k(Te, Xsec2.eniz);

% plot k
figure
loglog(Te,kenp1,'-r');
hold on
loglog(Te,keniz1,'-b');
loglog(Te,kenp2,'-.c');
loglog(Te,keniz2,'-.m');
ylabel('k');
xlabel('Te (eV)')
grid on
L1=legend('{\itk}_{en}^p-1990Tawara','{\itk}_{en}^{iz}-1990Tawara',...
    '{\itk}_{en}^p-Phelps','{\itk}_{en}^{iz}-Phelps');
set(L1,'FontSize',font_size);
set(L1,'location','southwest');
set(L1,'box','off')
set(L1,'AutoUpdate','off')
hold off

end
% test result: Xsec from 1990Tawara and Phelps are almost the same.
% details are recorded in .\others\Figures during code developing.pptx 

function test_get_Nagaoka(testCase)
% test 	get_Nagaoka
get_Nagaoka(nan);
hold on
a=0:1:11+0.005;
k=get_Nagaoka(a);
scatter(a,k)
legend('Nagaoka data','func get Nagaoka')
close
% 与1909Nagaoka - The inductance coefficients
% of solenoids中P31的表对比
tolerance=1/100;
verify_equal=@(actual, expected) verifyEqual(testCase,actual,expected,'RelTol',tolerance);
verify_equal(get_Nagaoka(0), 1)
verify_equal(get_Nagaoka(0.25), 0.901649)
verify_equal(get_Nagaoka(1), 0.688423)
verify_equal(get_Nagaoka(5), 0.319825)
verify_equal(get_Nagaoka(10), 0.203315)
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
