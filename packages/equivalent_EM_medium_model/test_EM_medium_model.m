%% Main function to generate tests
function tests = test_EM_medium_model
% test EM_medium_model
tests = functiontests(localfunctions);
end

%% Test Functions
function test_sigma_2019Rauner(testCase)
% compare sigma with Re(sigma_p) in the fig4.2 of 2019Rauner - Efficiency
% of RF plasma generation for fusion relevant ion sources
plasma.nu_eff=logspace(6,8,21);
plasma.ne=1e17;
plasma.wpe=get_omega_pe(plasma.ne);
plasma.wpi=get_omega_pi(plasma.ne,1,1);
flag.input_plasma='';
plasma.r=inf;

plasma.w_RF=2*pi*1e6;
plasma1  = equivalent_EM_medium_model( flag, plasma);
plasma.w_RF=2*pi*2e6;
plasma2  = equivalent_EM_medium_model( flag, plasma);
plasma.w_RF=2*pi*4e6;
plasma3  = equivalent_EM_medium_model( flag, plasma);
plasma.w_RF=2*pi*8e6;
plasma4  = equivalent_EM_medium_model( flag, plasma);

figure
loglog(plasma.nu_eff, plasma1.sigma, '--r')
hold on
loglog(plasma.nu_eff, plasma2.sigma, '--g')
loglog(plasma.nu_eff, plasma3.sigma, '--b')
loglog(plasma.nu_eff, plasma4.sigma, '--k')

L1=legend('w_RF=1MHz');
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
xlabel('\nu [Hz]');
ylabel('\sigma [S/m]');
axis([1e6,2e8,1,400])
grid on

end

% 没必要实现该test
function test_formula(testCase)
% % % test formula
% % flag = get_example_flag(0);
% % flag.electric_model='';
% % input=get_input_data(flag);
% % plasma=plasma_model(flag, input.plasma);
% 
% %%%%%%% 介电常数/电导率
% % plasma.sigma_medium=constants.eps0*plasma.nu_eff.*...
% %     plasma.wpe.^2./(plasma.w_RF.^2+plasma.nu_eff.^2); %与Re(sigmap)一致
% % plasma.eps_medium=constants.eps0*(1-...
% %     plasma.wpe.^2./(plasma.w_RF.^2+plasma.nu_eff.^2)); %与Re(epsc)一致
% % 
% % plasma.sigma_p=constants.eps0*plasma.wpe.^2./...
% %     (1i*plasma.w_RF+plasma.nu_eff); % 等离子体复电导率，即常见的sigma_p
% % plasma.sigma_c=plasma.sigma_p+1i*plasma.w_RF*constants.eps0; 
% 
% %         epsp_r_real(X1i,X2i)=1-plasma.wpe^2/(plasma.nu_eff^2+plasma.w_RF^2);         %复介电常数实部
% %         epsp_r_imag(X1i,X2i)=-plasma.nu_eff*plasma.wpe^2/(plasma.nu_eff^2+plasma.w_RF^2)/plasma.w_RF;         %复介电常数虚部
% %         sigma_from_epsp=-plasma.w_RF*constants.eps0*epsp_r_imag(X1i,X2i); %与sigmap_real一致
% %         %检验表达式用-断点调试
% %         sigmap_real(X1i,X2i)=constants.eps0*plasma.wpe^2*plasma.nu_eff/(plasma.nu_eff^2+plasma.w_RF^2);
% %         sigmap_imag(X1i,X2i)=constants.eps0*plasma.wpe^2*(-plasma.w_RF)/(plasma.nu_eff^2+plasma.w_RF^2);
% %         eps_r_from_sigmap(X1i,X2i)=sigmap_imag(X1i,X2i)/plasma.w_RF/constants.eps0+1; %与epsp_r_real一致
% %         tandelta_from_sigmap(X1i,X2i)=sigmap_real(X1i,X2i)/(plasma.w_RF*constants.eps0*eps_r_from_sigmap(X1i,X2i)); %与tandelta一致
% % fprintf('\n\n')
% % end
% 
% % function test_medium_approximation(testCase)
% % % test different medium_approximation
% % 
% % 
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
