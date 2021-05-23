%% Main function to generate tests
function tests = test_transformer_model
% test transformer model
tests = functiontests(localfunctions);
end

%% Test Functions
function test_PER_2018Jainb(testCase)
% compare PER with PER in the fig45 of 2018Jainb - Studies and experimental
% activities to qualify the behaviour of RF power circuits for Negative Ion
% Sources of Neutral Beam Injectors for ITER and fusion experiments

flag.input_plasma='2018Jainb_ELISE_sweep_f';
flag.type_Xsec='e-H2-Phelps';
flag.stoc_model='2018Jainb-simplify';
flag.medium_approximation='';
flag.skin_depth='as-medium-simplified-finite-radius';
flag.output_plasma_model=false;
flag.electric_model='transformer-2018Jainb';
flag.input_geometry='ELISE_base';
flag.Rmetal='measured-Rmetal-woplasma';
flag.Lcoil='calculated-Lcoil-woplasma';
flag.output_electric_model=false;
input1=get_input_data( flag );
input1.plasma=plasma_model(flag, input1.plasma);
p1=input1.plasma;
s1=electric_model( flag, input1 );
% PC
flag.stoc_model='Vahedi-simplify';
flag.skin_depth='as-medium';
flag.electric_model='transformer-base';
flag.Lcoil='measured-Lcoil-woplasma';
input2=get_input_data( flag );
input2.plasma=plasma_model(flag, input2.plasma);
p2=input2.plasma;
s2=electric_model( flag, input2 );

% plot
figure
legend_text={};
loglog(p1.f(:,1),s1.PER(:,1),'-.y');
legend_text{end+1}='Jain, ne=5e16';
hold on
loglog(p1.f(:,1),s1.PER(:,2),'-.b');
legend_text{end+1}='Jain, ne=5e17';
loglog(p1.f(:,1),s1.PER(:,3),'-.k');
legend_text{end+1}='Jain, ne=5e18';
loglog(p1.f(:,1),s1.PER(:,4),'-.m');
legend_text{end+1}='Jain, ne=5e19';
axis([5e5,5e7,3e0,2e2])
% PC
marker_indices=1:20:141;
idx=p2.ne(1,:)==5e16;
loglog(p2.f(:,1),s2.PER(:,idx),'-oy','MarkerIndices',marker_indices);
legend_text{end+1}='PC, ne=5e16';
idx=p2.ne(1,:)==5e17;
loglog(p2.f(:,1),s2.PER(:,idx),'-ob','MarkerIndices',marker_indices);
legend_text{end+1}='PC, ne=5e17';
idx=p2.ne(1,:)==5e18;
loglog(p2.f(:,1),s2.PER(:,idx),'-ok','MarkerIndices',marker_indices);
legend_text{end+1}='PC, ne=5e18';
idx=p2.ne(1,:)==5e19;
loglog(p2.f(:,1),s2.PER(:,idx),'-om','MarkerIndices',marker_indices);
legend_text{end+1}='PC, ne=5e19';
L1=legend(legend_text(1:5));
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
xlabel('f [Hz]');
ylabel('PER [\Omega]');
axis([-inf,inf,-inf,inf])
grid on

% plot
figure
legend_text={};
semilogx(p1.f(:,1),s1.PER(:,1),'-.y');
legend_text{end+1}='Jain, ne=5e16';
hold on
semilogx(p1.f(:,1),s1.PER(:,2),'-.b');
legend_text{end+1}='Jain, ne=5e17';
semilogx(p1.f(:,1),s1.PER(:,3),'-.k');
legend_text{end+1}='Jain, ne=5e18';
semilogx(p1.f(:,1),s1.PER(:,4),'-.m');
legend_text{end+1}='Jain, ne=5e19';
axis([5e5,5e7,3e0,2e2])
% PC
marker_indices=1:20:141;
idx=p2.ne(1,:)==5e16;
semilogx(p2.f(:,1),s2.PER(:,idx),'-oy','MarkerIndices',marker_indices);
legend_text{end+1}='PC, ne=5e16';
idx=p2.ne(1,:)==5e17;
semilogx(p2.f(:,1),s2.PER(:,idx),'-ob','MarkerIndices',marker_indices);
legend_text{end+1}='PC, ne=5e17';
idx=p2.ne(1,:)==5e18;
semilogx(p2.f(:,1),s2.PER(:,idx),'-ok','MarkerIndices',marker_indices);
legend_text{end+1}='PC, ne=5e18';
idx=p2.ne(1,:)==5e19;
semilogx(p2.f(:,1),s2.PER(:,idx),'-om','MarkerIndices',marker_indices);
legend_text{end+1}='PC, ne=5e19';
L1=legend(legend_text(1:5));
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
xlabel('f [Hz]');
ylabel('PER [\Omega]');
axis([-inf,inf,-inf,inf])
grid on

% typical
idx=p1.ne(1,:)==5e18 & p1.f(:,1)==1e6 ;
s1.PER(idx)
idx=p2.ne(1,:)==5e18 & p2.f(:,1)==1e6 ;
s2.PER(idx)

p2.sigma(idx)
p2.eps_r(idx)

end
% test result: PER for 2018Jain case in SAD is smaller than 
% details: .\others\Figures during code developing.pptx

function test_Rp_PER_Ar(testCase)
% compare the effect of Rp on PER in different models.
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
flag.input_plasma='given_directly';
input1.plasma.nu_m=input1.plasma.w_RF;
input1.plasma=plasma_model(flag, input1.plasma);
source1=electric_model( flag, input1 );

flag.input_plasma='2011Chabert';
flag.electric_model='transformer-base';
input2=get_input_data( flag );
flag.input_plasma='given_directly';
input2.plasma.nu_m=input2.plasma.w_RF;
input2.plasma=plasma_model(flag, input2.plasma);
source2=electric_model( flag, input2 );

flag.input_plasma='2011Chabert';
flag.electric_model='transformer-2018Jainb';
flag.skin_depth='as-medium-simplified-finite-radius'; 
input3=get_input_data( flag );
flag.input_plasma='given_directly';
input3.plasma.nu_m=input3.plasma.w_RF;
input3.plasma=plasma_model(flag, input3.plasma);
source3=electric_model( flag, input3 );

% Rp, PER
figure
yyaxis left
loglog(input1.plasma.ne, source1.Rp,'-b')
ylabel('R_p [\Omega]')
axis([1e14,1e19,6e-2,3e3])
yyaxis right
loglog(input1.plasma.ne, source1.PER,'--b')
ylabel('PER [\Omega]')
axis([1e14,1e19,7e-2,4e2])
hold on
yyaxis left
loglog(input1.plasma.ne, source1.transformer.Rp1,'-c')
loglog(input1.plasma.ne, source1.transformer.Rp2,'-m')
loglog(input1.plasma.ne, source2.transformer.Rp,'-.y')
loglog(input1.plasma.ne, source3.transformer.Rp,'-.g')
yyaxis right
loglog(input1.plasma.ne, source2.PER,'--y')
loglog(input1.plasma.ne, source3.PER,'--g')
xlabel('n_e [m^{-3}]')
grid on
L1=legend('Rp-a','Rp-low \nu-a','Rp-high \nu-a','Rp-t','Rp-tj','PER-a','PER-t','PER-tj');
set(L1,'location','best');
set(L1,'AutoUpdate','off');
end

% TODO: NIS_sweep

% function test_Rp_PER_H2(testCase)
% % compare the effect of Rp on PER in different models.
% flag.type_Xsec='e-Ar-Biagi';
% flag.input_plasma='2011Chabert';
% flag.stoc_model='';
% flag.medium_approximation='';
% flag.skin_depth='as-medium-simplified'; 
% flag.electric_model='analytical_base';
% flag.input_geometry='2011Chabert';
% flag.Rmetal='calculated-Rcoil-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';
% input1=get_input_data( flag );
% flag.input_plasma='given_directly';
% input1.plasma.nu_m=input1.plasma.w_RF;
% input1.plasma=plasma_model(flag, input1.plasma);
% source1=electric_model( flag, input1 );
% 
% flag.input_plasma='2011Chabert';
% flag.electric_model='transformer-base';
% input2=get_input_data( flag );
% flag.input_plasma='given_directly';
% input2.plasma.nu_m=input2.plasma.w_RF;
% input2.plasma=plasma_model(flag, input2.plasma);
% source2=electric_model( flag, input2 );
% 
% flag.input_plasma='2011Chabert';
% flag.electric_model='transformer-2018Jainb';
% flag.skin_depth='as-medium-simplified-finite-radius'; 
% input3=get_input_data( flag );
% flag.input_plasma='given_directly';
% input3.plasma.nu_m=input3.plasma.w_RF;
% input3.plasma=plasma_model(flag, input3.plasma);
% source3=electric_model( flag, input3 );
% 
% % Rp, PER
% figure
% yyaxis left
% loglog(input1.plasma.ne, source1.Rp,'-b')
% ylabel('R_p [\Omega]')
% axis([1e14,1e19,6e-2,3e3])
% yyaxis right
% loglog(input1.plasma.ne, source1.PER,'--b')
% ylabel('PER [\Omega]')
% axis([1e14,1e19,7e-2,4e2])
% hold on
% yyaxis left
% loglog(input1.plasma.ne, source1.transformer.Rp1,'-c')
% loglog(input1.plasma.ne, source1.transformer.Rp2,'-m')
% loglog(input1.plasma.ne, source2.transformer.Rp,'-.y')
% loglog(input1.plasma.ne, source3.transformer.Rp,'-.g')
% yyaxis right
% loglog(input1.plasma.ne, source2.PER,'--y')
% loglog(input1.plasma.ne, source3.PER,'--g')
% xlabel('n_e [m^{-3}]')
% grid on
% L1=legend('Rp-a','Rp-low \nu-a','Rp-high \nu-a','Rp-t','Rp-tj','PER-a','PER-t','PER-tj');
% set(L1,'location','best');
% set(L1,'AutoUpdate','off');
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