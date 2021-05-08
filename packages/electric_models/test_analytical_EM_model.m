%% Main function to generate tests
function tests = test_analytical_EM_model
% test analytical_EM_model
tests = functiontests(localfunctions);
end

%% Test Functions
function test_formula(testCase)
flag.type_Xsec='e-Ar-Biagi';
flag.input_plasma='2011Chabert';
flag.stoc_model='';
flag.medium_approximation='';
flag.skin_depth='as-medium-simplified';
flag.electric_model='analytical_base';
flag.input_geometry='2011Chabert';
flag.Rmetal='calculated-Rcoil-woplasma';
flag.Lcoil='calculated-Lcoil-woplasma';
input=get_input_data( flag );
input.plasma=plasma_model(flag, input.plasma);
source=electric_model( flag, input );

delta_PQ=source.PQplasma-source.PQ0;
disp('∫S|r=r_coil - ∫S|r=r_plasma: ')
disp(delta_PQ)
end
% result: r_plasma<r<r_coil范围内有 Q=- 2.7434e+02 W, 该值与等离子体参数无关
% Pabs in 2011Chabert=Pplasma in SAD code

% Pcp_part1=k_p_wave*r_p*besselj(1,k_p_wave*r_p)/...
%     (w_RF*constants.eps0*epsp_r*besselj(0,k_p_wave*r_p));
% Pcp_part2=Pcp_part1+w_RF*constants.mu0*(r_coil^2-r_p^2)/2;
% Pcp=1i*pi*N_equ^2*Icoil_rms^2*Pcp_part2/l_equ;
% % ERROR 该表达式与简化表达式结果不一致

%                 %检验表达式用-断点调试
%                 Pplasma1_part1=1i*k_p_wave*r_p*besselj(1,k_p_wave*r_p)/(w_RF*constants.eps0*epsp_r*besselj(0,k_p_wave*r_p));
%                 Pplasma1=pi*N_equ^2*I_coil^2*real(Pplasma1_part1)/l_equ; %与Pplasma一致，说明上下表达式自洽


%                 %检验表达式用-断点调试
% 与LZS实现结果一致
%                 Rp1_part1=1i*k_p_wave*r_p*besselj(1,k_p_wave*r_p)/(w_RF*constants.eps0*epsp_r*besselj(0,k_p_wave*r_p));
% 可能少了一个N^2
%                 Rp1=2*pi*real(Rp1_part1)/l_equ/(1/besselj(0,k_p_wave*r_p)-1)^2; %与Rp一致，说明上下表达式自洽

function test_E_H_2011Chabert(testCase)
% compare E/H with E/H in the fig7.3 of 2011Chabert - Physics of
% Radio-Frequency Plasmas
flag.type_Xsec='e-Ar-Biagi';
flag.input_plasma='2011Chabert';
flag.stoc_model='';
flag.medium_approximation='';
flag.skin_depth='as-medium-simplified';
flag.electric_model='analytical_base';
flag.input_geometry='2011Chabert';
flag.Rmetal='calculated-Rcoil-woplasma';
flag.Lcoil='calculated-Lcoil-woplasma';
input1=get_input_data( flag );
input2=input1;
input1.plasma=plasma_model(flag, input1.plasma);
source1=electric_model( flag, input1 );

Xsec=get_Xsec('e-Ar-Biagi',true);
axis([0.5,10,-inf,inf])
idx=Xsec.enp(:,1)>0.5 & Xsec.enp(:,1)<10; 
geometric_mean=trapz(Xsec.enp(idx,1),Xsec.enp(idx,2))/(10-0.5);
line([0.5,10],[geometric_mean,geometric_mean])
input2.plasma.nu_m=geometric_mean*input1.plasma.ng*input1.plasma.ve;
flag.input_plasma='given_directly';
input2.plasma=plasma_model(flag, input2.plasma);
source2=electric_model( flag, input2 );

num_ne=length(input1.plasma.ne);
r=0:5e-3:input1.geometry.r_chamber;
num_r=length(r);
normH1=zeros(num_r,num_ne);
normE1=normH1;
normH2=normH1;
normE2=normH1;
for i=1:num_r
    normH1(i,:)=abs(source1.emf.Hz_plasma(r(i)))./ abs(source1.emf.Hz_plasma(r(end)));
    normE1(i,:)=abs(source1.emf.Etheta_plasma(r(i)))./abs(source1.emf.Etheta_plasma(r(end)));
    normH2(i,:)=abs(source2.emf.Hz_plasma(r(i)))./ abs(source2.emf.Hz_plasma(r(end)));
    normE2(i,:)=abs(source2.emf.Etheta_plasma(r(i)))./abs(source2.emf.Etheta_plasma(r(end)));
end

idx1=input1.plasma.ne==1e15;
idx2=input1.plasma.ne==5e16;
idx3=input1.plasma.ne==1e18;
% idx=idx1 | idx2 | idx3; 
figure
plot(r,normH1(:,idx1),'-.r')
hold on
plot(r,normH1(:,idx2),'-.b')
plot(r,normH1(:,idx3),'-.y')
plot(r,normH2(:,idx1),'--c')
plot(r,normH2(:,idx2),'--m')
plot(r,normH2(:,idx3),'--g')
xlabel('r [m]')
ylabel('norm |H|')
grid on
axis([0,input1.geometry.r_chamber,0,1.2])
L1=legend('n_e=1e15','n_e=5e16','n_e=1e18');
set(L1,'location','best');
set(L1,'AutoUpdate','off');

figure
plot(r,normE1(:,idx1),'-.r')
hold on
plot(r,normE1(:,idx2),'-.b')
plot(r,normE1(:,idx3),'-.m')
xlabel('r [m]')
ylabel('norm |E|')
grid on
axis([0,input1.geometry.r_chamber,0,1.2])
L1=legend('n_e=1e15','n_e=5e16','n_e=1e18');
set(L1,'location','best');
set(L1,'AutoUpdate','off');
end

function test_analytical_EM_model_2011Chabert(testCase)
% compare results of analytical EM model with results in the fig7.6/7/8 of
% 2011Chabert - Physics of Radio-Frequency Plasmas
flag.type_Xsec='e-Ar-Biagi';
flag.input_plasma='2011Chabert';
flag.stoc_model='';
flag.medium_approximation='';
flag.skin_depth='collisionless';
flag.electric_model='analytical_base';
flag.input_geometry='2011Chabert';
flag.Rmetal='calculated-Rcoil-woplasma';
flag.Lcoil='calculated-Lcoil-woplasma';
input1=get_input_data( flag );
input2=input1;
input3=input1;
input4=input1;
input1.plasma=plasma_model(flag, input1.plasma);
source1=electric_model( flag, input1 );

flag.input_plasma='given_directly';
input2.plasma.nu_m=0.1*input2.plasma.w_RF;
input2.plasma=plasma_model(flag, input2.plasma);
source2=electric_model( flag, input2 );
input3.plasma.nu_m=input3.plasma.w_RF;
input3.plasma=plasma_model(flag, input3.plasma);
source3=electric_model( flag, input3 );
input4.plasma.nu_m=10*input4.plasma.w_RF;
input4.plasma=plasma_model(flag, input4.plasma);
source4=electric_model( flag, input4 );

% fig 7.6
figure
loglog(input1.plasma.ne, source2.Pplasma,'-.r')
hold on
loglog(input1.plasma.ne, source3.Pplasma,'-.b')
loglog(input1.plasma.ne, source4.Pplasma,'-.y')
xlabel('n_e [m^{-3}]')
ylabel('Power [W]')
grid on
axis([1e14,1e19,1e-1,5e2])
L1=legend('\nu/\omega=0.1','\nu/\omega=1','\nu/\omega=10');
set(L1,'location','best');
set(L1,'AutoUpdate','off');

% fig 7.7
figure
semilogx(input1.plasma.ne, 1e6*source3.Lind,'-.r')
hold on
semilogx(input1.plasma.ne, 1e6*source3.Lm,'-.b')
xlabel('n_e [m^{-3}]')
ylabel('L [\muH]')
grid on
axis([1e14,1e19,0,3])
L1=legend('L_{ind}','L_m');
set(L1,'location','best');
set(L1,'AutoUpdate','off');

% fig 7.8
Rp_conductor=2*pi*input3.geometry.r_plasma_eff./...
    (input3.plasma.sigma_dc.*input3.plasma.skin_depth*input3.geometry.l_plasma); 
skin_depth=get_plasma_skin_depth('as-medium',...
    input3.plasma.f,input3.plasma.nu_eff,input3.plasma.wpe,input3.plasma.r);
Rp_conductor_different_delta=2*pi*input3.geometry.r_plasma_eff./...
    (input3.plasma.sigma_dc.*skin_depth*input3.geometry.l_plasma); 
figure
yyaxis left
loglog(input1.plasma.ne, sqrt(2)*source3.Iplasma_rms,'-.r')
ylabel('I_p [A]')
axis([1e14,1e19,1e-2,1e2])
yyaxis right
loglog(input1.plasma.ne, source3.Rp,'-.b')
ylabel('R_p [\Omega]')
axis([1e14,1e19,1e-1,1e3])
hold on
loglog(input1.plasma.ne, source3.transformer.Rp1,'--c')
loglog(input1.plasma.ne, source3.transformer.Rp2,'--m')
loglog(input1.plasma.ne, Rp_conductor,'--y')
loglog(input1.plasma.ne, Rp_conductor_different_delta,'--g')
xlabel('n_e [m^{-3}]')
grid on
L1=legend('I_p','R_p-analytical','R_p-low pressure','R_p-high pressure','R_p-conductor','R_p-conductor, \delta use medium');
set(L1,'location','best');
set(L1,'AutoUpdate','off');

figure
semilogx(input1.plasma.ne, source3.PER./source3.Rp,'-.r')
hold on
line([1e14,1e19],[input1.geometry.N_coil^2,input1.geometry.N_coil^2]);
semilogx(input1.plasma.ne, source3.PER./source3.Rp,'-.r')
xlabel('n_e [m^{-3}]')
ylabel('PER/R_p')
grid on
end

% function test_small_source1_LZS(testCase)
% % test small_source1_LZS
% 
% flag = get_flag( );
% flag.input_plasma='small_source1_LZS';
% flag.electric_model='analytical_base';
% input=get_input_data( flag );
% plasma=input.plasma;
% % ICP heating model: 使用给定的nu_m、nu_st数据
% % 李增山-整体模型耦合解析电模型-2020.03.30\中间数据作为输入，用于解析电模型benchmark
% disp('使用李增山-整体模型耦合解析电模型-2020.03.30的vm和vst')
% plasma.nu_m=1.1238e7;
% plasma.nu_st=3.6537e7;
% fprintf('%s = %.2e , ','[INFO] Results from LZS: ν_m= , ν_st= ',plasma.nu_m,plasma.nu_st);
% 
% flag.input_plasma='given_directly';
% plasma=plasma_model(flag, plasma);
% 
% end


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