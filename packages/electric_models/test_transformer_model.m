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

end
% test result: PER for 2018Jain case in DSA is smaller than 
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