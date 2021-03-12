%% Main function to generate tests
function tests = test_analytical_EM_model
% test analytical_EM_model
tests = functiontests(localfunctions);
end

%% Test Functions
function test_small_source1_LZS(testCase)
% test small_source1_LZS

flag = get_flag( );
flag.input_plasma='small_source1_LZS';
flag.electric_model='analytical_base';
input=get_input_data( flag );
plasma=input.plasma;
% ICP heating model: 使用给定的nu_m、nu_st数据
% 李增山-整体模型耦合解析电模型-2020.03.30\中间数据作为输入，用于解析电模型benchmark
disp('使用李增山-整体模型耦合解析电模型-2020.03.30的vm和vst')
plasma.nu_m=1.1238e7;
plasma.nu_st=3.6537e7;
fprintf('%s = %.2e , ','[INFO] Results from LZS: ν_m= , ν_st= ',plasma.nu_m,plasma.nu_st);

flag.input_plasma='given_directly';
plasma=plasma_model(flag, plasma);

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