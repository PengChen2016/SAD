%% Main function to generate tests
function tests = test_plasma_model
% test plasma_model
tests = functiontests(localfunctions);
end

% test result: 20210407 @pengchen2016: 
% Plasma model is ok, according to test_compare_with_v1.
% equivalent_EM_medium_model is ok, too.

%% Test Functions
function test_single_data_point(testCase)
% test single_data_point
flag = get_example_flag(0);
flag.electric_model='';
input=get_input_data(flag);

flag.output_plasma_model=true;
plasma=plasma_model(flag, input.plasma);

% 电模型需要的输入条件
verifyEqual(testCase,size(plasma.mu_r),[1,1])
verifyEqual(testCase,size(plasma.eps_r),[1,1])
verifyEqual(testCase,size(plasma.sigma),[1,1])

fprintf('\n\n')
end
% test result: ok.

function test_multi_data_point(testCase)
% test multi_data_point
flag = get_example_flag(2);
flag.output_plasma_model=true;
input=get_input_data(flag);
plasma=plasma_model(flag, input.plasma);

% 电模型需要的输入条件
verifyEqual(testCase,size(plasma.eps_r),size(plasma.ne))
verifyEqual(testCase,size(plasma.sigma),size(plasma.ne))
end
% test result: ok.

function test_compare_with_v1(testCase)
% compare with fig 4-6 for code v1 paper v1
flag=get_example_flag(2);
flag.skin_depth='as-medium-simplified';
input=get_input_data(flag);
plasma=plasma_model(flag, input.plasma);

% fig 5b
figure
loglog(plasma.ne(:,1),plasma.nu_enp(:,3),'--c','LineWidth',2);
hold on
loglog(plasma.ne(:,1),plasma.nu_eniz(:,3),'--m','LineWidth',2);
loglog(plasma.ne(:,1),plasma.nu_eip(:,3),'--y','LineWidth',2);
loglog(plasma.ne(:,1),plasma.nu_m(:,3),'--r','LineWidth',2);
axis([1e16,1e19,1e4,5e7])
ylabel('\nu');
xlabel('{\itn}_e')
grid on
legend('\it{\bf\nu}\rm_{en}^{(p)}','\it{\bf\nu}\rm_{en}^{(iz)}','\it{\bf\nu}\rm_{ei}^{(p)}')

% fig 5a
figure
loglog(plasma.ne(:,1),plasma.nu_m(:,3),'--c','LineWidth',2);
hold on
loglog(plasma.ne(:,1),plasma.nu_st(:,3),'--m','LineWidth',2);
axis([1e16,1e19,2e6,2e8])
ylabel('\nu');
xlabel('{\itn}_e')
grid on
legend('\it{\bf\nu}\rm_{m}','\it{\bf\nu}\rm_{st}')

% fig 4b
figure
yyaxis left
loglog(plasma.ne(:,1),-plasma.eps_r(:,1),'-.r','LineWidth',2);
axis([1e16,1e19,1e5,1e7])
name_Y2='-\it{\bf\epsilon}\rm_r';
ylabel(name_Y2);
yyaxis right
loglog(plasma.ne(:,1),plasma.sigma(:,1),'-.k','LineWidth',2);
axis([1e16,1e19,1e1,1e4])
hold on
yyaxis left
loglog(plasma.ne(:,1),-plasma.eps_r(:,3),'--m','LineWidth',2);
yyaxis right
loglog(plasma.ne(:,1),plasma.sigma(:,3),'--c','LineWidth',2);
hold on
name_Y1='\it{\bf\sigma}\rm';
ylabel([name_Y1 ' [S/m]']);
xlabel('{\itn}_e')
grid on
L1=legend([name_Y2 ', 5eV'],[name_Y2 ',15eV'],...
    [name_Y1 ', 5eV'],[name_Y1 ',15eV']);
set(L1,'Location','best')

% fig 6
figure
loglog(plasma.ne(:,1),plasma.skin_depth(:,3),'--c','LineWidth',2);
hold on
loglog(plasma.ne(:,1),plasma.wavelength(:,3),'--m','LineWidth',2);
axis([1e16,1e19,5e-3,1e0])
ylabel('length');
xlabel('{\itn}_e')
grid on
legend('\delta','\lambda')

% fig 7
figure
loglog(plasma.ne(:,1),plasma.wpe(:,3),'--c','LineWidth',2);
hold on
loglog(plasma.ne(:,1),plasma.wpi(:,3),'--m','LineWidth',2);
ylabel('\nu');
xlabel('{\itn}_e')
grid on
axis([1e16,1e19,1e6,1e12])
L1=legend('\it{\bf\omega}_{\rmpi}','\it{\bf\omega}_{\rmpe}');
set(L1,'location','northwest');
end
% test result: ok.
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
