%% Main function to generate tests
function tests = test_plasma_model
% test plasma_model
tests = functiontests(localfunctions);
end

%% Test Functions
function test_basic(testCase)
% test formula

flag = get_example_flag(0);
input=get_input_data(flag);
plasma=plasma_model(flag, input.plasma);

% 电模型需要的输入条件
size(plasma.eps_c_r)
size(plasma.sigma)
end

function test_formula(testCase)
% test formula
flag = get_example_flag(0);
input=get_input_data(flag);
plasma=plasma_model(flag, input.plasma);

%%%%%%% 介电常数/电导率
% plasma.sigma_medium=constants.eps0*plasma.nu_eff.*...
%     plasma.wpe.^2./(plasma.w_RF.^2+plasma.nu_eff.^2); %与Re(sigmap)一致
% plasma.eps_medium=constants.eps0*(1-...
%     plasma.wpe.^2./(plasma.w_RF.^2+plasma.nu_eff.^2)); %与Re(epsc)一致
% 
% plasma.sigma_p=constants.eps0*plasma.wpe.^2./...
%     (1i*plasma.w_RF+plasma.nu_eff); % 等离子体复电导率，即常见的sigma_p
% plasma.sigma_c=plasma.sigma_p+1i*plasma.w_RF*constants.eps0; 

%         epsp_r_real(X1i,X2i)=1-plasma.wpe^2/(plasma.nu_eff^2+plasma.w_RF^2);         %复介电常数实部
%         epsp_r_imag(X1i,X2i)=-plasma.nu_eff*plasma.wpe^2/(plasma.nu_eff^2+plasma.w_RF^2)/plasma.w_RF;         %复介电常数虚部
%         sigma_from_epsp=-plasma.w_RF*constants.eps0*epsp_r_imag(X1i,X2i); %与sigmap_real一致
%         %检验表达式用-断点调试
%         sigmap_real(X1i,X2i)=constants.eps0*plasma.wpe^2*plasma.nu_eff/(plasma.nu_eff^2+plasma.w_RF^2);
%         sigmap_imag(X1i,X2i)=constants.eps0*plasma.wpe^2*(-plasma.w_RF)/(plasma.nu_eff^2+plasma.w_RF^2);
%         eps_r_from_sigmap(X1i,X2i)=sigmap_imag(X1i,X2i)/plasma.w_RF/constants.eps0+1; %与epsp_r_real一致
%         tandelta_from_sigmap(X1i,X2i)=sigmap_real(X1i,X2i)/(plasma.w_RF*constants.eps0*eps_r_from_sigmap(X1i,X2i)); %与tandelta一致

end

function test_medium_approximation(testCase)
% test different medium_approximation


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
