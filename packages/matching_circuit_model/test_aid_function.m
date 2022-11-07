% disp(table(runtests('test_aid_function.m')))
%% Main function to generate tests
function tests = test_aid_function
% test aid_function
tests = functiontests(localfunctions);
end

%% Test Functions
function test_get_impedance_mode1(testCase)
% 2021Jain eP-210607-01 双驱串并联融合\期刊论文\Records of input and output data processing.xlsx
% RD [Ω]	LD [μH]	XD [Ω]	|ZD| [Ω]	cosφ	φ [°]
% 1.9	19.25	121.0 	121.0 	0.0157 	89.10 
f=1e6;
RD=1.9;
LD=19.25e-6;
XD=121;
ZD=RD+1i*XD;
phi=89.1;
% test type Z and keep other field of Z struct
Z.cm=ZD;
Z.other='keep';
Z=get_impedance('Z',Z,f);
verifyEqual(testCase,Z.L,LD,'RelTol',1e-3);
assert(strcmp(Z.other,'keep'))
clear Z
% test type RL
Z.R=RD;
Z.L=LD;
Z=get_impedance('R_L',Z,f);
verifyEqual(testCase,Z.phase,phi,'AbsTol',1e-3);
clear Z
% test type ap
Z.abs=abs(ZD);
Z.phase=phi;
Z=get_impedance('a_p',Z);
verifyEqual(testCase,Z.R,RD,'AbsTol',1e-3);
clear Z
% test type aR
Z.abs=abs(ZD);
Z.R=RD;
Z=get_impedance('aR',Z);
verifyEqual(testCase,Z.cm,ZD,'AbsTol',1e-3);
clear Z
end

function test_get_impedance_mode2(testCase)
% test type C2X and X2C
f=1e6;
w=2*pi*f;
C0=1e-9;
X0=-1/(w*C0);
X=get_impedance( 'C', C0, f );
verifyEqual(testCase,X,X0,'AbsTol',1e-9);
C=get_impedance( 'x2c', X, f );
verifyEqual(testCase,C,C0,'RelTol',1e-5);
% test type s
Z=[2+2i, 2+2i];
verifyEqual(testCase,get_impedance( 's', Z ),4+4i,'RelTol',1e-5);
% test type p
verifyEqual(testCase,get_impedance( 'p', Z ),1+1i,'RelTol',1e-5);
end

function test_get_impedance_RC(testCase)
f=1e6;
ZCs.R=0.192;
ZCs.C=1.42e-9;
C0=ZCs.C;
% 直接用 mode1-RC
ZCs=get_impedance('RC',ZCs,f);
C1=ZCs.C;
Z1=ZCs.cm;
% 先用mode2-C2X，然后mode1-RX
ZCs.X=get_impedance('C2X',ZCs.C,f);
ZCs=get_impedance('RX',ZCs,f);
C2=ZCs.C;
Z2=ZCs.cm;
verifyEqual(testCase,C1,C0,'RelTol',1e-5);
verifyEqual(testCase,C2,C0,'RelTol',1e-5);
verifyEqual(testCase,Z1,Z2,'RelTol',1e-5);
end

function test_get_UIP(testCase)
% test get_UIP
% 2021Jain eP-210607-01 双驱串并联融合\期刊论文\Records of input and output data processing.xlsx
PD=1e5;
ZD=3.8+118.06i;
phase_ZD=phase(ZD)*180/pi;
ID_rms=162.2214211;
UD_rms=19161.94967;
QD=3106869.787;
% test type UI
s.U.abs=sqrt(2)*UD_rms;
s.I.abs=sqrt(2)*ID_rms;
s=get_UIP( 'UI', s, phase_ZD );
verifyEqual(testCase,s.P,PD,'relTol',1e-3);
clear s
% test type IZ
s1.Z.cm=ZD;
s1.I.abs=sqrt(2)*ID_rms;
s1=get_UIP( 'IZ', s1);
verifyEqual(testCase,s1.Q,QD,'relTol',1e-3);
% test type PZ
s.Z.cm=ZD;
s.P=PD;
s=get_UIP( 'PZ', s, 360);
verifyEqual(testCase,s.I.rms,ID_rms,'relTol',1e-3);
verifyEqual(testCase,s.U.cm,s1.U.cm,'relTol',1e-3);
clear s
end

function test_get_RF_mode1(testCase)
% test type Z
% eP-210603-02 基于小源结合实验的研究\实验记录\S20210729-小源调试激发\S20220325 6匝单螺旋线圈\S20220414 6匝单螺旋线圈2MHz氢气试激发\匹配单元-离线冷测\
% 网分测得 2MHz：46.69+j1.0109 Ω，|S11|=-28.925dB
s.Z.cm=46.69+1.0109i;
s.Z0=50;
f=2e6;
s=get_RF( 'Z', s, f );
verifyEqual(testCase,s.S11dB,-28.925,'RelTol',1e-3);

% TODO: other types

% TODO: UIP

end

function test_get_RF_mode2(testCase)
% eP-210607-01 双驱串并联融合.docx 仅abs 换算表
s.gamma.abs=0.5;
s=get_RF( 'gamma_abs', s );
verifyEqual(testCase,s.VSWR,3,'RelTol',1e-3);
s.VSWR=2;
s=get_RF( 'VSWR', s );
verifyEqual(testCase,s.gamma.abs,1/3,'RelTol',1e-3);
verifyEqual(testCase,s.S11dB,-4.77*2,'RelTol',1e-3);
s.S11dB=-20;
s=get_RF( 'S11dB', s );
verifyEqual(testCase,s.gamma.abs,0.1,'RelTol',1e-3);
end

function test_get_dB(testCase)
% 
assert(get_dB( 'P',1)==0)
verifyEqual(testCase,get_dB( '',0.5),-3.01,'RelTol',1e-3);
assert(get_dB( 'dBm',1)==30)
assert(get_dB( 'U',1)==0)

assert(get_dB( 'dBmi',30)==1)
verifyEqual(testCase,get_dB( 'i',-3.01),0.5,'RelTol',1e-3);
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
