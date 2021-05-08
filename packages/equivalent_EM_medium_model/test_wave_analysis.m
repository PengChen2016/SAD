%% Main function to generate tests
function tests = test_wave_analysis
% test wave_analysis
tests = functiontests(localfunctions);
end

% test result: 20210413 @pengchen2016: 
% wave_analysis is ok, according to test_different_skin_depth_2018Jainb /
% test_different_wavelength_vary

%% Test Functions
function test_get_skin_depth_plasma_basic(testCase)
% test get_skin_depth_plasma basic
% single point
flag=get_example_flag(0);
input=get_input_data(flag);
plasma=ICP_heating_model(flag, input.plasma);
delta1=get_plasma_skin_depth('as-medium',...
    plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);
delta2=get_plasma_skin_depth('as-medium-simplified',...
    plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);
verifyEqual(testCase,delta1,delta2,'RelTol',1e-4);

% multi point
constants=get_constants();
flag=get_example_flag(2);
input=get_input_data(flag);
plasma=ICP_heating_model(flag, input.plasma);
delta=get_plasma_skin_depth('collisionless',...
    plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);
verifyEqual(testCase,delta,constants.c./plasma.wpe);
end
% test result: For typical parameters set, 'as-medium' is almost the same
% as 'as-medium-simplified'

function test_different_eps_c(testCase)
% test the wave of different eps_c
flag.input_plasma='2018Jainb_ELISE_sweep_f';
flag.stoc_model='2018Jainb-simplify';
flag.type_Xsec='e-H2-Phelps';
flag.skin_depth='';
input=get_input_data(flag);
input.plasma.r=0.14; % ELISE_case
plasma=ICP_heating_model(flag, input.plasma);
% 1. origin, - & -
plasma1=equivalent_EM_medium_model( flag, plasma);
% 2. modified, + & +
plasma2=plasma1;
plasma2.eps_c=-plasma1.eps_c;
plasma2=wave_analysis( plasma2, flag.skin_depth);
% 3. modified, + & -
plasma3=plasma1;
plasma3.eps_c=-real(plasma1.eps_c)+1i*imag(plasma1.eps_c);
plasma3=wave_analysis( plasma3, flag.skin_depth);
% 4. modified, - & +
plasma4=plasma1;
plasma4.eps_c=real(plasma1.eps_c)-1i*imag(plasma1.eps_c);
plasma4=wave_analysis( plasma4, flag.skin_depth);

idx1=find(plasma.ne(1,:)==5e18,1);
% skin_depth
figure
semilogx(plasma.f(:,1),plasma1.skin_depth(:,idx1),'-r');
hold on
semilogx(plasma.f(:,1),plasma2.skin_depth(:,idx1),'-.c');
semilogx(plasma.f(:,1),plasma3.skin_depth(:,idx1),'--m');
semilogx(plasma.f(:,1),plasma4.skin_depth(:,idx1),'.b');
ylabel('\delta')
grid on
xlabel('f')
legend('1. origin, - & -','2. + & +','3. + & -','4. - & +')

figure
loglog(plasma.f(:,1),plasma1.skin_depth(:,idx1),'-r');
hold on
loglog(plasma.f(:,1),plasma3.skin_depth(:,idx1),'--m');
ylabel('\delta')
grid on
xlabel('f')
legend('1. origin, - & -','3. + & -')

% wavelength
figure
loglog(plasma.f(:,1),plasma1.wavelength(:,idx1),'-r');
hold on
loglog(plasma.f(:,1),plasma2.wavelength(:,idx1),'-.c');
loglog(plasma.f(:,1),plasma3.wavelength(:,idx1),'--m');
loglog(plasma.f(:,1),plasma4.wavelength(:,idx1),'.b');
ylabel('\lambda')
grid on
xlabel('f')
legend('1. origin, - & -','2. + & +','3. + & -','4. - & +')
end

function test_different_skin_depth_2018Jainb(testCase)
% compare results from get_skin_depth_plasma() with results from fig.42 of
% 2018Jainb - Studies and experimental activities to qualify the behaviour
% of RF power circuits for Negative Ion Sources of Neutral Beam Injectors
% for ITER and fusion experiments
flag.input_plasma='2018Jainb_ELISE_sweep_f';
flag.stoc_model='2018Jainb-simplify';
flag.type_Xsec='e-H2-Phelps';
input=get_input_data(flag);
input.plasma.r=0.14; % ELISE_case
plasma=ICP_heating_model(flag, input.plasma);
% 1. result from wave number, lossy medium
flag.skin_depth='';
plasma=equivalent_EM_medium_model( flag, plasma);
delta1=plasma.skin_depth;
% 2. result of type 'as-medium' 
% is supposed to equal the result from wave number
delta2=get_plasma_skin_depth('as-medium',...
        plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);
% 3. result of type 'as-medium-simplified'
delta3=get_plasma_skin_depth('as-medium-simplified',...
        plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);
% 4. result of type 'as-medium-simplified'
delta4=get_plasma_skin_depth('collisionless',...
        plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);
% 5. result of type 'as-medium-simplified'
delta5=get_plasma_skin_depth('as-medium-simplified-finite-radius',...
        plasma.f,plasma.nu_eff,plasma.wpe,plasma.r);
% 6. result from sigma, good conductor
delta6=sqrt(2./(plasma.w_RF.*real(plasma.mu_c).*plasma.sigma));
% 7. result from wave number, lossy medium modified as capacitive dielectric
plasma2=plasma;
plasma2.eps_c=-real(plasma.eps_c)+1i*imag(plasma.eps_c);
plasma2=wave_analysis( plasma2, flag.skin_depth);
delta7=plasma2.skin_depth;
    
% ne=5e16
idx1=find(plasma.ne(1,:)==5e16);
figure
semilogx(plasma.f(:,1),delta1(:,idx1),'-.y','LineWidth',6);
hold on
semilogx(plasma.f(:,1),delta2(:,idx1),'-.c','LineWidth',5);
semilogx(plasma.f(:,1),delta3(:,idx1),'--m','LineWidth',2);
semilogx(plasma.f(:,1),delta4(:,idx1),'-.g');
semilogx(plasma.f(:,1),delta5(:,idx1),'-.k');
semilogx(plasma.f(:,1),delta6(:,idx1),'-.b');
semilogx(plasma.f(:,1),delta7(:,idx1),'--r','LineWidth',2);
ylabel('\delta_{skin}')
grid on
legend('default medium','as-medium','as-medium-simplified',...
    'collisionless','as-medium-simplified-finite-radius',...
    'as-good-conductor','medium modified as capacitive dielectric')
axis([1e4,1e10,0,35e-2])
xlabel('f')
title('2018Jainb ne=5e16')

figure
yyaxis left
loglog(plasma.f(:,1),-real(plasma.eps_c(:,idx1)),'-r');
yyaxis right
semilogx(plasma.f(:,1),-real(plasma.eps_c(:,idx1)),'--r');
hold on
yyaxis left
loglog(plasma.f(:,1),-imag(plasma.eps_c(:,idx1)),'-k');
yyaxis right
semilogx(plasma.f(:,1),-imag(plasma.eps_c(:,idx1)),'--k');
legend('lg(-Re(\epsilon_c))','lg(-Im(\epsilon_c))','-Re(\epsilon_c)','-Im(\epsilon_c)')
grid on
xlabel('f')

% ne=5e18
idx2=find(plasma.ne(1,:)==5e18);
figure
semilogx(plasma.f(:,1),delta1(:,idx2),'-.y','LineWidth',6);
hold on
semilogx(plasma.f(:,1),delta2(:,idx2),'-.c','LineWidth',5);
semilogx(plasma.f(:,1),delta3(:,idx2),'--m','LineWidth',2);
semilogx(plasma.f(:,1),delta4(:,idx2),'-.g');
semilogx(plasma.f(:,1),delta5(:,idx2),'-.k');
semilogx(plasma.f(:,1),delta6(:,idx2),'-.b');
semilogx(plasma.f(:,1),delta7(:,idx2),'--r','LineWidth',2);
grid on
ylabel('\delta_{skin}')
legend('default medium','as-medium','as-medium-simplified',...
    'collisionless','as-medium-simplified-finite-radius',...
    'as-good-conductor','medium modified as capacitive dielectric')
axis([1e4,1e10,0,6e-2])
xlabel('f')
title('2018Jainb ne=5e18')
end
% test result: get_skin_depth_plasma is ok. 
% details: .\others\Figures during code developing.pptx

function test_different_wavelength_vary(testCase)
% test different wavelength vary with different parameters
flag.input_plasma='2018Jainb_ELISE_sweep_f';
flag.stoc_model='2018Jainb-simplify';
flag.type_Xsec='e-H2-Phelps';
input=get_input_data(flag);
input.plasma.r=0.14; % ELISE_case
plasma=ICP_heating_model(flag, input.plasma);
% 1. result of lossy medium, default
flag.skin_depth='';
plasma1=equivalent_EM_medium_model( flag, plasma);
lambda1=plasma1.wavelength;
% 2. result of lossless dielectric
wave_speed=1./sqrt(real(plasma1.mu_c).*(-real(plasma1.eps_c)));
lambda2=wave_speed./plasma1.f;
% 3. result of lossy medium, modified as capacitive dielectric
plasma2=plasma;
plasma2.eps_c=-real(plasma1.eps_c)+1i*imag(plasma1.eps_c);
plasma2.mu_c=plasma1.mu_c;
plasma2=wave_analysis( plasma2, flag.skin_depth);
lambda3=plasma2.wavelength;

% f=1e6, sweep ne
idx1=find(plasma.f(:,1)==1e6);
figure
yyaxis left
loglog(plasma.ne(1,:),lambda1(idx1,:),'-r');
yyaxis right
loglog(plasma.ne(1,:),real(plasma1.eps_c(idx1,:))./imag(plasma1.eps_c(idx1,:)),'-.c');
hold on
yyaxis left
loglog(plasma.ne(1,:),lambda2(idx1,:),'-.k');
loglog(plasma.ne(1,:),lambda3(idx1,:),'--m');
ylabel('\lambda')
grid on
xlabel('n_{e}')
legend('as-lossy-medium','as-lossless-dielectric',...
    'as-lossy-medium, modified as capacitive dielectric',...
    'lg(Re(\epsilon_c)/Im(\epsilon_c))')

% ne=5e18, sweep f
idx2=find(plasma.ne(1,:)==5e18);
figure
yyaxis left
loglog(plasma.f(:,1),lambda1(:,idx2),'-r');
ylabel('\lambda')
yyaxis right
loglog(plasma.f(:,1),real(plasma1.eps_c(:,idx2))./imag(plasma1.eps_c(:,idx2)),'--c');
hold on
yyaxis left
loglog(plasma.f(:,1),lambda2(:,idx2),'--k');
loglog(plasma.f(:,1),lambda3(:,idx2),'-.m');
grid on
xlabel('f')
legend('as-lossy-medium','as-lossless-dielectric',...
    'as-lossy-medium, modified as capacitive dielectric',...
    'lg(Re(\epsilon_c)/Im(\epsilon_c))')
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
