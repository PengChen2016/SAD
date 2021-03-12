%% Main function to generate tests
function tests = test_transformer_model
% test transformer model
tests = functiontests(localfunctions);
end

%% Test Functions
function test_2018Jain(testCase)
% test 2018Jain

% 使用外部的nu_m、nu_st数据
% 2018Jain中间数据作为变压器输入，用于变压器模型benchmark
'ELISE_base'
plasma.nu_m=5.6e6;
plasma.nu_st=4.3e7;
fprintf('%s = %.2e , ','[INFO] Results from 2018Jain: ν_m= , ν_st= ',plasma.nu_m,plasma.nu_st);


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