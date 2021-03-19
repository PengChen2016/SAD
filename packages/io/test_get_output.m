%% Main function to generate tests
function tests = test_get_output
% test get_output
tests = functiontests(localfunctions);
end

%% Test Functions
function test_get_output_json(testCase)
% test get_output_json
% 多级结构体 ok
flag = get_flag();
input=get_input_data(flag);
get_output_json( input, 'test_output-input');
% 多维矩阵 ok，但文件过大
flag.input_plasma='2020Chen_NIS_sweep_p';
input=get_input_data(flag);
plasma=plasma_model(flag, input.plasma);
get_output_json( plasma, 'test_output-plasma');
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
