%% Main function to generate tests
function tests = test_stoc_model
% test stoc_model
tests = functiontests(localfunctions);
end

%% Test Functions
function test_Cazzador_fit(testCase)
% test Cazzador_fit

vst_new_f_fit_fun1=@(x)vst_fun_Cazzador(x*w_RF_Cazzador^2/w_RF^2)*w_RF/w_RF_Cazzador; %测试确定，该变换表达式正确
% 测试确定，通过变换方法得到的结果合理

va_fun_Cazzador=@(X)vst_fun_Cazzador(K0*w_RF_Cazzador^2*X)/w_RF_Cazzador; %与频率无关的va-X关系
vst_fun_new_f=@(x)va_fun_Cazzador(x/K0/w_RF^2)*w_RF;

%                     test_temp=vst_new_f_fit_fun1(plasma.ne*plasma.Te);
% 测试驱动频率对Cazzador fit关系的影响
%         vst4_bug=vst_fit_Cazzador(plasma.ne*plasma.Te); %有较大差别

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