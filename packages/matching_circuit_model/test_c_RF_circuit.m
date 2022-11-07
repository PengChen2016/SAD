% disp(table(runtests('test_c_RF_circuit.m')))
%% Main function to generate tests
function tests = test_c_RF_circuit
% test c_RF_circuit
tests = functiontests(localfunctions);
end

%% Test Functions
function test_Constructor(testCase)
% test Constructor method
% 默认f
c0=c_RF_circuit;
assert(c0.D.U.cm==0) % 多层引用
assert(c0.s.Z0==50)
% 输入f
f=2e6;
c0=c_RF_circuit(f);
assert(c0.w==2*pi*f)
end

function test_init_source_port(testCase)
% test 使用内部s和内部P
c0=c_RF_circuit;
c0.s.Z.cm=50;
c0.s.P=1000;
c0=c0.init_source_port('Z');
verifyEqual(testCase,c0.s.I.cm,sqrt(2*1000/50))

% test 使用外部s和无P
s.Z.cm=50;
s.P=0;
c0=c0.init_source_port('Z',s);
verifyEqual(testCase,c0.s.I.cm,sqrt(2*1/50))

% test 使用外部s和外部P
s.Z.cm=50;
c0=c0.init_source_port('Z',s,100, 75);

verifyEqual(testCase,c0.s.I.cm,sqrt(2*100/50))
end

function test_init_ZD(testCase)
% test init_ZD
% 2021Jain eP-210607-01 双驱串并联融合\期刊论文\Records of input and output data processing.xlsx
% RD [Ω]	LD [μH]	XD [Ω]	|ZD| [Ω]	cosφ	φ [°]
% 1.9	19.25	121.0 	121.0 	0.0157 	89.10 
RD=1.9;
LD=19.25e-6;
XD=121;
ZD=RD+1i*XD;
phi=89.1;

% test 单驱+使用内部Z
c0=c_RF_circuit;
c0.D.Z.cm=ZD;
c0=c0.init_ZD('Z');
verifyEqual(testCase,c0.D.Z.phase,phi,'AbsTol',1e-3);
c0=c_RF_circuit;
c0.D.Z.R=RD;
c0.D.Z.X=XD;
c0=c0.init_ZD('R_X',1);
verifyEqual(testCase,c0.D.Z.phase,phi,'AbsTol',1e-3);
% test 单驱+使用外部Z
c0=c_RF_circuit;
Z.cm=ZD;
c0=c0.init_ZD('Z',1,Z);
verifyEqual(testCase,c0.D.Z.phase,phi,'AbsTol',1e-3);
% test 多驱s+使用内部Z
c0=c_RF_circuit;
c0.Dlist(1).Z.cm=ZD;
c0.Dlist(2).Z.cm=ZD;
c0=c0.init_ZD('Z',2);
verifyEqual(testCase,c0.D.Z.cm,2*ZD,'RelTol',1e-3);
% test 多驱p+使用外部Z
c0=c_RF_circuit;
Z(1).R=RD;
Z(1).L=LD;
Z(2).R=RD;
Z(2).L=LD;
f=2e6;
w=2*pi*f;
ZD=RD+1i*w*LD;
c0=c0.init_ZD('R_L',2,Z,0,f,'p');
verifyEqual(testCase,c0.D.Z.cm,ZD/2,'RelTol',1e-3);
verifyEqual(testCase,c0.Dlist(1).Z.cm,ZD,'RelTol',1e-3);
verifyEqual(testCase,c0.Dlist(2).Z.phase,phase(ZD)*180/pi,'AbsTol',1e-3);
% test 多驱p+使用外部Z+ZD_leadline
ZD_leadline.cm=0.2;
c0=c0.init_ZD('R_L',2,Z,ZD_leadline,f,'p');
verifyEqual(testCase,c0.D.Z.cm,ZD/2+0.2,'RelTol',1e-3);
c0=c0.init_ZD('R_L',2,Z,0.1+0.1i,f,'p');
verifyEqual(testCase,c0.D.Z.cm,ZD/2+0.1+0.1i,'RelTol',1e-3);
end

function test_init_matching_t(testCase)
% test init_matching_t
% 小源-周磊数据
L1 = (92.28)*10^-6;             %变压器原边电感
L2 = (10.38)*10^-6;              %变压器副边线圈电感
k =0.9826 ;                     %互感耦合系数  
M = k*(L1*L2)^0.5;              %互感

% test ''
c0=c_RF_circuit;
c0.m.t.n=3;
c0=c0.init_matching_t('');
assert(c0.m.t.n==1)
% test ideal
c0=c0.init_matching_t('i',3);
assert(c0.m.t.n==3)
% test mutual
c0=c0.init_matching_t('m',M,L1,L2);
verifyEqual(testCase,c0.m.t.k,k,'RelTol',1e-3);
end

function test_init_matching_C(testCase)
% test init_matching_C
f=1e6;
c0=c_RF_circuit(f);
% test single C
C1=1e-9;
C=c0.init_C(C1);
assert(C1==get_impedance( 'X2C', f, C.Z.X ))
% test multi C with UI_withstand
C2=2e-9;
Cs=[C1 C1 C1];
Cp=[C1 C2];
Cs_UI_withstand=[1e3 100; 2e3 100; 2e3 50];
Cp_UI_withstand=[2e3 100; 1e3 50];
c0=c0.init_matching_C( 'gamma', Cs, Cp, Cs_UI_withstand, Cp_UI_withstand);
assert(c0.m.p.U_withstand==1e3)
verifyEqual(testCase,c0.m.s.I_withstand,150)
verifyEqual(testCase,c0.m.s.Z.C,3*C1,'RelTol',1e-5)
end

%% test get_C/get_Zs/get_ZD/get_all
% 对于 get_C/get_Zs/get_ZD，前三种情况 与Multisim结果一致，ok 2022/05/06 16:35:17
% 第四种情况与Multisim尚未对比

% 对于get_all：% 正反能对起来，ok 2022/06/09 14:55:47

function test_2009Zamengo(testCase)
% 2009Zamengo: t.model='';
f=1e6;
Cs0=1.5e-9; % 完美匹配时电容
Cp0=10e-9;

c0=c_RF_circuit(f);
% init matching t
c0=c0.init_matching_t('');
% -------------- get ZD
c0.s.Z.cm=c0.s.Z0;
c0=c0.init_matching_C( 'gamma', Cs0, Cp0);
c0=c0.get_ZD();
fprintf('ZD = %.2f Ω+ %.2f μH\n',c0.D.Z.R,c0.D.Z.L*1e6)
% ZD = 4.60 Ω+ 19.19 μH，在2009Zamengo范围内，且与2018Maistrello单激励器Z接近
verifyEqual(testCase,c0.D.Z.R,2.3*2,'RelTol',1e-3)
RD0=c0.D.Z.R;
LD0=c0.D.Z.L;
% -------------- get C
% init ZD
Z.R=RD0;
Z.L=LD0;
c0=c0.init_ZD('RL',1,Z);
c0.s.Z0=50;
c0=c0.get_C();
verifyEqual(testCase,c0.m.s.Z.C,Cs0,'RelTol',1e-3)
verifyEqual(testCase,c0.m.p.Z.C,Cp0,'RelTol',1e-3)
% -------------- get Zs
c0=c0.get_Zs();
verifyEqual(testCase,c0.s.Z.cm,50,'RelTol',1e-3)
% -------------- get all D2s
c0.D.P=1e3;
c0=c0.get_all('PD');
verifyEqual(testCase,c0.s.Z.cm,50,'RelTol',1e-3)
Iratio=c0.s.I.cm/c0.D.I.cm;
% -------------- get all D2s
c0=c0.get_all('Ps');
verifyEqual(testCase,c0.s.I.cm/c0.D.I.cm,Iratio,'RelTol',1e-3)
verifyEqual(testCase,c0.D.P,1e3,'RelTol',1e-3)
end
% result: get_Zs/get_ZD/get_C互相能对起来,ok 2022/05/06 14:25:52

function test_2016Yue_tIdeal(testCase)
% 2016岳海昆: t.model='ideal'
% 数据源：大源，1MHz – 匹配设计 2016岳海昆 ch5.2.2
f=1e6;
RD0=1.8;
LD0=20.5e-6;
Cs0=1.26e-9; % 完美匹配时电容
Cp0=4.6e-9;

c0=c_RF_circuit(f);
% init matching t
c0=c0.init_matching_t('ideal',3);
% -------------- get C
% init ZD
Z.R=RD0;
Z.L=LD0;
c0=c0.init_ZD('RL',1,Z);
c0.s.Z0=50;
c0=c0.get_C();
verifyEqual(testCase,c0.m.s.Z.C,Cs0,'RelTol',1e-2)
verifyEqual(testCase,c0.m.p.Z.C,Cp0,'RelTol',1e-2)
Cs_precise=c0.m.s.Z.C;
Cp_precise=c0.m.p.Z.C;
% -------------- get Zs
% ------- 使用精准电容值
c0=c0.get_Zs();
verifyEqual(testCase,c0.s.Z.cm,50,'RelTol',1e-3)
% ------- 使用四舍五入后电容值
c0=c0.init_matching_C( 'gamma', Cs0, Cp0);
c0=c0.get_Zs();
verifyEqual(testCase,c0.s.Z.cm,50,'RelTol',1e-1)
fprintf('当Cs相差%.2e, Cp相差%.2e, 求得Zs会相差%.2e+j %.2e (abs=%.2f)\n',...
    Cs_precise-c0.m.s.Z.C,Cp_precise-c0.m.p.Z.C,50-c0.s.Z.R,c0.s.Z.X,abs(50-c0.s.Z.cm))
% -------------- get ZD
% ------- 使用四舍五入后电容值
c0.s.Z.cm=c0.s.Z0;
c0=c0.get_ZD();
verifyEqual(testCase,c0.D.Z.R,RD0,'RelTol',1e-2)
verifyEqual(testCase,c0.D.Z.L,LD0,'RelTol',1e-2)
RD1=c0.D.Z.R;
XD1=c0.D.Z.X;
ZD1=c0.D.Z.cm;
% ------- 使用精准电容值
c0=c0.init_matching_C( 'gamma', Cs_precise, Cp_precise);
c0=c0.get_ZD();
verifyEqual(testCase,c0.D.Z.R,RD0,'RelTol',1e-3)
verifyEqual(testCase,c0.D.Z.L,LD0,'RelTol',1e-3)
fprintf('当Cs相差%.2e, Cp相差%.2e, 求得ZD会相差%.2e+j %.2e (abs=%.2f)\n',...
    c0.m.s.Z.C-Cs0, c0.m.p.Z.C-Cp0, c0.D.Z.R-RD1, c0.D.Z.X-XD1, abs(c0.D.Z.cm-ZD1))
% -------------- get all D2s
c0.D.P=1e3;
c0=c0.get_all('PD');
verifyEqual(testCase,c0.s.Z.cm,50,'RelTol',1e-3)
Iratio=c0.s.I.cm/c0.D.I.cm;
% -------------- get all D2s
c0=c0.get_all('Ps');
verifyEqual(testCase,c0.s.I.cm/c0.D.I.cm,Iratio,'RelTol',1e-3)
verifyEqual(testCase,c0.D.P,1e3,'RelTol',1e-3)
end
% result: get_Zs/get_ZD/get_C与2016岳海昆计算结果一致,ok 2022/05/06 14:25:52

function test_2016Yue_tMutual(testCase)
% 2016岳海昆: t.model='mutual'
f=1e6;
c0=c_RF_circuit(f);
% init ZD
RD0=1.8;
LD0=20.5e-6;
Z.R=RD0;
Z.L=LD0;
c0=c0.init_ZD('RL',1,Z);
% init matching t
L1 = 135.4e-6;             %变压器原边电感
L2 = 16.8e-6;              %变压器副边线圈电感
k =0.96 ;                     %互感耦合系数  
M = k*(L1*L2)^0.5;              %互感
c0=c0.init_matching_t('m',M,L1,L2);
verifyEqual(testCase,c0.m.t.M,45.8e-6,'RelTol',1e-3);
% -------------- get C
c0.s.Z0=50;
c0=c0.get_C();
% 用get到的C，去get Zs与ZD，若get Zs与ZD ok，说明get C也ok
% -------------- get Zs
c0=c0.get_Zs();
verifyEqual(testCase,c0.s.Z.cm,50,'RelTol',1e-3)
% -------------- get ZD
c0.s.Z.cm=c0.s.Z0;
c0=c0.get_ZD();
verifyEqual(testCase,c0.D.Z.R,RD0,'RelTol',1e-2)
verifyEqual(testCase,c0.D.Z.L,LD0,'RelTol',1e-2)
% -------------- get all D2s
c0.D.P=1e3;
c0=c0.get_all('PD');
verifyEqual(testCase,c0.s.Z.cm,50,'RelTol',1e-3)
Iratio1=c0.s.I.cm/c0.D.I.cm;
Uscm1=c0.s.U.cm;
% -------------- get all s2D
c0=c0.get_all('Ps');
Iratio2=c0.s.I.cm/c0.D.I.cm;
verifyEqual(testCase,Iratio2,Iratio1,'RelTol',1e-3)
verifyEqual(testCase,c0.D.P,1e3,'RelTol',1e-3)
% -------------- get all D2s
c0=c0.get_all('IDrms');
Uscm2=c0.s.U.cm;
verifyEqual(testCase,Uscm2,Uscm1,'RelTol',1e-3)
end
% result: get_Zs/get_ZD/get_C互相能对起来,ok 2022/05/06 14:25:52
% get_all的四种type均ok。

function test_2022Chen_HI_tMutual_R(testCase)
% m.t.model='mutual'
% m.R1>0 || m.R2>0
f=1e6;
c0=c_RF_circuit(f);
% init ZD
RD0=0.382;
LD0=3.16e-6;
Z.R=RD0;
Z.L=LD0;
c0=c0.init_ZD('RL',1,Z);
% init matching t
L1 = 45.307e-6;             %变压器原边电感
L2 = 10.646e-6;              %变压器副边线圈电感
k =0.9826 ;                     %互感耦合系数  
M = k*(L1*L2)^0.5;              %互感
c0=c0.init_matching_t('m',M,L1,L2);

c0=c0.init_matching_R(0.125, 0);

% -------------- get C
c0.s.Z0=50;
c0=c0.get_C();
% 与MMA-匹配公式推导.nb结果比较
verifyEqual(testCase,c0.m.s.Z.X,-20.0423,'RelTol',1e-2)
verifyEqual(testCase,c0.m.p.Z.X,-9.38968,'RelTol',1e-2)
% -------------- get Zs
c0=c0.get_Zs();
verifyEqual(testCase,c0.s.Z.cm,50,'RelTol',1e-3)
% -------------- get ZD
c0.s.Z.cm=c0.s.Z0;
c0=c0.get_ZD();
verifyEqual(testCase,c0.D.Z.R,RD0,'RelTol',1e-2)
verifyEqual(testCase,c0.D.Z.L,LD0,'RelTol',1e-2)
% -------------- get all D2s
c0.D.P=1e3;
c0=c0.get_all('PD');
verifyEqual(testCase,c0.s.Z.cm,50,'RelTol',1e-3)
Iratio=c0.s.I.cm/c0.D.I.cm;
% -------------- get all s2D
c0=c0.get_all('Ps');
verifyEqual(testCase,c0.s.I.cm/c0.D.I.cm,Iratio,'RelTol',1e-3)
verifyEqual(testCase,c0.D.P,1e3,'RelTol',1e-3)

% 检验Icoil表达式
Icoil_temp1=c0.w^2*(c0.m.t.L1 - c0.s.Z.L - c0.m.p.Z.C*(c0.m.R1*c0.s.Z.R + c0.w^2*c0.m.t.L1*c0.s.Z.L))^2;
Icoil_temp2=(c0.s.Z.R-c0.m.R1- c0.w^2*c0.m.p.Z.C*(c0.s.Z.L*c0.m.R1+c0.m.t.L1*c0.s.Z.R))^2;
Icoil=sqrt((Icoil_temp1+Icoil_temp2)*c0.s.P/c0.s.Z.R)/(c0.w*c0.m.t.M);
verifyEqual(testCase,Icoil,c0.D.I.rms,'RelTol',1e-3)
end
% result: get_Zs/get_ZD/get_C互相能对起来,get C与MMA结果一致，ok  2022/06/08 18:44:29

%% test multiD
function test_get_PUI_multiD(testCase)
% test 多驱p+使用外部Z
f=1e6;
c0=c_RF_circuit(f);
ZD_list(1).R=3.8/2;
ZD_list(1).L=18.79e-6/2;
ZD_list(2).R=1.9/2;
ZD_list(2).L=19.25e-6/2;
c0=c0.init_ZD('RL',2,ZD_list,0,f,'p');

c0.D.P=1e5;
c0=c0.get_PUI_multiD('PZ');

verifyEqual(testCase,c0.Dlist(1).I.rms,188.78499,'RelTol',1e-3);
verifyEqual(testCase,c0.Dlist(2).U.rms,11149.84831,'RelTol',1e-3);
verifyEqual(testCase,c0.Dlist(1).P,67715.56767,'RelTol',1e-3);
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
