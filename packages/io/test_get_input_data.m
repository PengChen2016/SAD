%% Main function to generate tests
function tests = test_get_input_data
% test get_input_data
tests = functiontests(localfunctions);
end

%% Test Functions
function test_get_plasma_single(testCase)
% test single point
plasma=get_input_plasma( '2018Jainb_ELISE_typical' );
verifyEqual(testCase,plasma.size,1)
end

function test_get_plasma_multi_decoupled_2D(testCase)
% test multi point with parameters decoupled
plasma=get_input_plasma( '2020Chen_NIS_sweep1' );
dne=16:1:19;ne=10.^dne;
Te=5:5:25;
verifyTrue(testCase,isempty(plasma.ne(plasma.ne(end,:)~=ne(end))))
verifyTrue(testCase,isempty(plasma.Te(plasma.Te(:,end-1)~=Te(end-1))))
end

function test_get_plasma_multi_decoupled_3D(testCase)
% test multi point with parameters decoupled
plasma=get_input_plasma( '2020Chen_NIS_sweep_p' );
p=0.3:0.3:10;
dne=16:0.1:19;
ne=10.^dne;
Te=1:4:21;
verifyEqual(testCase,plasma.size,[length(p),1,1,...
length(ne),length(Te)])

verifyTrue(testCase,isempty(plasma.p(  plasma.p(end,:,:)~=p(end)  )))
verifyTrue(testCase,isempty(plasma.Te(plasma.Te(:,:,end-1)~=Te(end-1))))
end

function test_get_plasma_multi_coupled(testCase)
% test multi point with parameters coupled
plasma=get_input_plasma( '2019Raunera_CHARLIE_sweep' );
verifyTrue(testCase,isempty(plasma.p(plasma.p(end,:,:)~=10)))
verifyTrue(testCase,isempty(plasma.f(plasma.f(:,end,:)~=4e6)))
verifyEqual(testCase,plasma.size,[size(plasma.ne),1])
end

function test_get_external_formula(testCase)
% test get_input_external formula
input_plasma=get_input_plasma( '2018Jainb_ELISE_typical' );
input_geo=get_input_geometry( 'HUST_small_driver_ZL' );
flag.electric_model='not empty';
flag.Rmetal='measured-Rmetal-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
external=get_input_external( flag, input_geo, input_plasma.w_RF );
% 可见，公式计算值与测量值有较大区别
verifyEqual(testCase,external.Rmetal,external.Rmetal_ex)
verifyEqual(testCase,external.Lcoil,external.Lcoil_ex)
end

function test_get_external_other_case(testCase)
% test get_input_external_other_case
input_plasma=get_input_plasma( '2018Jainb_ELISE_typical' );
input_geo=get_input_geometry( 'ELISE_base' );
flag.electric_model='not empty';
flag.Rmetal='measured-Rmetal-woplasma';
flag.Lcoil='calculated-Lcoil-woplasma';
external=get_input_external( flag, input_geo, input_plasma.w_RF );
verifyEqual(testCase,external.Rmetal,external.Rcoil_ex)
verifyEqual(testCase,external.Lcoil,external.Lcoil_th)
% 可见，公式计算值与测量值有较大区别
end

function test_get_data_base(testCase)
% test get_input_data base
flag.input_plasma='2018Jainb_ELISE_typical';
flag.input_geometry='ELISE_base';
flag.electric_model='not empty';
flag.Rmetal='measured-Rmetal-woplasma';
flag.Lcoil='calculated-Lcoil-woplasma';
input=get_input_data( flag );
verifyEqual(testCase,input.geometry.r_plasma_eff,input.geometry.r_chamber)
end

function test_get_data_CHARLIE_single(testCase)
% test get_input_data CHARLIE_single case
flag.input_plasma='CHARLIE_10Pa_4MHz_520W';
flag.input_geometry='CHARLIE_base';
flag.Rmetal='measured-Rmetal-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
flag.electric_model='not empty';
input=get_input_data( flag );

verifyEqual(testCase,input.plasma.ne,3.7e17)
verifyEqual(testCase,input.external.Rmetal,0.213659)
end

function test_CHARLIE_raza_sweep(testCase)
flag.input_plasma='CHARLIE_raza_sweep';
input=get_input_data( flag );
actual_plasma=input.plasma;

excepted_plasma=get_input_plasma( '2019Raunera_CHARLIE_sweep' );
ratio_origin2goal.ne_r=nonuniform_dist.get_ne_r([0,45.5])/nonuniform_dist.get_ne_r(10);
excepted_plasma.ne= excepted_plasma.ne*ratio_origin2goal.ne_r;
constants=get_constants();
excepted_plasma.wpe=get_omega_pe(excepted_plasma.ne); 
excepted_plasma.wpi=get_omega_pi(excepted_plasma.ne,1,1); %离子等离子体频率
excepted_plasma.ve=sqrt(8*excepted_plasma.Te*constants.e/(pi*constants.me));

verifyEqual(testCase,actual_plasma,excepted_plasma)
end
% test result: ok  20210420

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