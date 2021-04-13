%% Main function to generate tests
function tests = test_ICP_heating_model
% test ICP_heating_model
tests = functiontests(localfunctions);
end

% test result: 20210406 @pengchen2016: 

% 1. Ohmic_heating_model is ok, according to
% test_nu_m_2018Jainb & test_nu_m_2014Cazzador.
% 2. stochastic_heating_model is ok, according to
% test_nu_st_2018Jainb

%% Test Functions
function test_nu_m_2018Jainb(testCase)
% compare ν_m with ν_m in the fig38 of 2018Jainb - Studies and experimental
% activities to qualify the behaviour of RF power circuits for Negative Ion
% Sources of Neutral Beam Injectors for ITER and fusion experiments

% get_input_plasma
constants=get_constants();
plasma.p=0.3;
plasma.Tg=1200;
plasma.ng=plasma.p./(constants.kB*plasma.Tg);
plasma.ne=5e18;
plasma.Te=logspace(-1,2,40);
plasma.ve=sqrt(8*plasma.Te*constants.e/(pi*constants.me)); 
plasma= Ohmic_heating_model( plasma, 'Phelps-m' );
plasma1= Ohmic_heating_model( plasma, '1990Tawara' );

% plot
figure
loglog(plasma.Te,plasma.nu_enp,'-r');
hold on
loglog(plasma.Te,plasma.nu_eniz,'--r');
loglog(plasma.Te,plasma.nu_eip,'-.r');
loglog(plasma.Te,plasma.nu_m,'-b');
loglog(plasma.Te,plasma1.nu_enp,'-.c','LineWidth',1.5);
loglog(plasma.Te,plasma1.nu_eniz,'-.m','LineWidth',1.5);
xlabel('{\itT}_e[eV]');
ylabel('collision frequency');
L1=legend('{\it\nu}_{enp}-Phelps-m','{\it\nu}_{eniz}-Phelps-m',...
    '{\it\nu}_{eip}','{\it\nu}_{m}',...
    '{\it\nu}_{enp}-1990Tawara','{\it\nu}_{eniz}-1990Tawara');
set(L1,'Location','best');
title('Ohmic heating model');
grid on%显示网格
axis([0.1,100,1e5,1e8]) %绘图显示范围，即[xmin,xmax,ymin,ymax]

% compare
% verifyEqual(testCase,plasma.nu_eip(plasma.Te==1),3.38e6,'RelTol',0.01);
end
% test result: 20210406 @pengchen2016:
% nu_eip of two code are the same.
% The differences of nu_enp and nu_eniz may be due to different cross
% sections used. 
% details: .\others\Figures during code developing.pptx

function test_nu_m_2014Cazzador(testCase)
% compare ν_m with ν_m in the fig3.1 of 2014Cazzador - Analytical and
% numerical models and first operations on the negative ion source NIO1

% get_input_plasma
constants=get_constants();
plasma.p=1;
plasma.Tg=400;
plasma.ng=plasma.p./(constants.kB*plasma.Tg);
alpha=1e-4; % ionization degree
plasma.ne=alpha*plasma.ng/(1-alpha)/2; % H2->2H+
% plasma.ne=alpha*plasma.ng/(1-alpha);
plasma.Te=logspace(-1,2,40);
plasma.ve=sqrt(8*plasma.Te*constants.e/(pi*constants.me)); 
plasma= Ohmic_heating_model( plasma, 'Phelps-m' );

% plot
figure
loglog(plasma.Te,plasma.nu_enp,'-r');
hold on
loglog(plasma.Te,plasma.nu_eniz,'--r');
loglog(plasma.Te,plasma.nu_eip,'-.r');
loglog(plasma.Te,plasma.nu_m,'-b');
xlabel('{\itT}_e[eV]');
ylabel('collision frequency');
L1=legend('{\it\nu}_{enp}','{\it\nu}_{eniz}','{\it\nu}_{eip}','{\it\nu}_{m}');
set(L1,'Location','best');
title('Ohmic heating model');
grid on%显示网格
axis([0.1,100,1e6,1e8]) %绘图显示范围，即[xmin,xmax,ymin,ymax]

% compare
% verifyEqual(testCase,plasma.nu_eip(plasma.Te==1),3.38e6,'RelTol',0.01);
end
% test result: 20210406 @pengchen2016:
% After code reviewed, there are still differences between nu_eniz and
% nu_eip of this code and of fig3.1 from 2014Cazzador. 
% The difference of nu_eniz may be explained by the difference of cross
% sections data source.
% Why the difference of nu_eip? Maybe the ne used is different.
% details: .\others\Figures during code developing.pptx

function test_nu_st_2018Jainb(testCase)
% compare ν_st with ν_st in the fig40 of 2018Jainb - Studies and experimental
% activities to qualify the behaviour of RF power circuits for Negative Ion
% Sources of Neutral Beam Injectors for ITER and fusion experiments

% get_input_plasma
flag.input_plasma='2018Jainb_ELISE_sweep_f';
input=get_input_data( flag );
plasma1=stochastic_heating_model('2018Jainb-simplify', input.plasma);
plasma2=stochastic_heating_model('Vahedi-simplify', input.plasma);
% plasma2=stochastic_heating_model('Cazzador-fit', plasma);

% plot
figure
loglog(plasma1.f(:,1),plasma1.nu_st(:,1),'-r');
hold on
loglog(plasma1.f(:,1),plasma1.nu_st(:,2),'-k');
loglog(plasma1.f(:,1),plasma1.nu_st(:,3),'-','Color',[0.5,0.5,1]);
loglog(plasma1.f(:,1),plasma1.nu_st(:,4),'-','Color',[0.5,0.5,0.5]);
h_plot5=loglog(plasma2.f(:,1),plasma2.nu_st(:,1),'-.c');
h_plot6=loglog(plasma2.f(:,1),plasma2.nu_st(:,2),'-.m');
h_plot7=loglog(plasma2.f(:,1),plasma2.nu_st(:,3),'-.y');
h_plot8=loglog(plasma2.f(:,1),plasma2.nu_st(:,4),'-.b');
xlabel('{\itf}[Hz]');
ylabel('{\it\nu}_{st}[Hz]');
grid on%显示网格
axis([1e2,1e10,1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
L1=legend('{\it\nu}_{st}-Jain, ne=5e16','{\it\nu}_{st}-Jain, ne=5e17',...
    '{\it\nu}_{st}-Jain, ne=5e18','{\it\nu}_{st}-Jain, ne=5e19');
set(L1,'Location','best');
axes2 = axes('position',get(gca,'position'),'visible','off');
L2=legend(axes2, [h_plot5,h_plot6,h_plot7,h_plot8],...
    '{\it\nu}_{st}-Vahedi, ne=5e16','{\it\nu}_{st}-Vahedi, ne=5e17',...
    '{\it\nu}_{st}-Vahedi, ne=5e18','{\it\nu}_{st}-Vahedi, ne=5e19');
% L2=legend(axes2, [h_plot5,h_plot6,h_plot7,h_plot8],...
%     '{\it\nu}_{st}-Cazzador, ne=5e16','{\it\nu}_{st}-Cazzador, ne=5e17',...
%     '{\it\nu}_{st}-Cazzador, ne=5e18','{\it\nu}_{st}-Cazzador, ne=5e19');
set(L2,'location','best');
end
% test result: 20210408 @pengchen2016:
% There are differences.
% details: .\others\Figures during code developing.pptx

function test_nu_2019Raunera(testCase)
% compare ν with ν in the fig7 of 2019Raunera - Influence of the excitation
% frequency on the RF power transfer efficiency of low pressure hydrogen
% ICPs

flag.input_plasma='2019Raunera_CHARLIE_sweep';
flag.stoc_model='Vahedi-simplify';
input=get_input_data( flag );
plasma=ICP_heating_model( flag, input.plasma);

% plot 1MHz
figure
loglog(plasma.p(:,1),plasma.nu_st(:,1),'--r');
hold on
loglog(plasma.p(:,1),plasma.nu_m(:,1),'--k');
loglog(plasma.p(:,1),plasma.nu_eff(:,1),'--b');
xlabel('{\itp}[Pa]');
ylabel('{\it\nu}[Hz]');
grid on%显示网格
axis([0.3,10,2e6,2e8])
L1=legend('{\it\nu}_{st}','{\it\nu}_{m}','{\it\nu}_{eff}');
set(L1,'location','best');

% plot 4MHz
figure
loglog(plasma.p(:,1),plasma.nu_st(:,2),'--r');
hold on
loglog(plasma.p(:,1),plasma.nu_m(:,2),'--k');
loglog(plasma.p(:,1),plasma.nu_eff(:,2),'--b');
xlabel('{\itp}[Pa]');
ylabel('{\it\nu}[Hz]');
grid on%显示网格
axis([0.3,10,2e6,2e8])
L1=legend('{\it\nu}_{st}','{\it\nu}_{m}','{\it\nu}_{eff}');
set(L1,'location','best');
end
% test result: almost the same.
% details: .\others\Figures during code developing.pptx

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
