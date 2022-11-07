% 小源6匝单 匹配电路模型
% S20211223-放电前匹配电路调试\
clear
addpath(genpath('../../packages'))

%% 分析实测数据
%--------------------- 实验数据分析
% % 1MHz 2022陈鹏 测量
% f=1e6;
% disp('1MHz')
% disp('考虑匹配单元中电阻')
% Z_sCsL2=4.42675+47.24955351i;
% Z_pCpL1=0.06419-9.862428697i;
% Z_Cs.cm=0.93351-19.64407778i;
% Z_Cp.cm=0.049415-9.532185742i;
% Z_L2.cm=Z_sCsL2-Z_Cs.cm;
% Z_L1.cm=1/(1/Z_pCpL1-1/Z_Cp.cm);
% Z_Cs=get_impedance('Z',Z_Cs,f);
% Z_Cp=get_impedance('Z',Z_Cp,f);
% Z_L1=get_impedance('Z',Z_L1,f);
% Z_L2=get_impedance('Z',Z_L2,f);
% disp('忽略匹配单元中电阻')
% Z_sCsL2=imag(Z_sCsL2)*1i;
% Z_pCpL1=imag(Z_pCpL1)*1i;
% Z_Cs.cm=imag(Z_Cs.cm)*1i;
% Z_Cp.cm=imag(Z_Cp.cm)*1i;
% Z_L2.cm=Z_sCsL2-Z_Cs.cm;
% Z_L1.cm=1/(1/Z_pCpL1-1/Z_Cp.cm);
% Z_Cs=get_impedance('Z',Z_Cs,f)
% Z_Cp=get_impedance('Z',Z_Cp,f)
% Z_L1=get_impedance('Z',Z_L1,f)
% Z_L2=get_impedance('Z',Z_L2,f)

% 2MHz 2022陈鹏 测量
% f=2e6;
% disp('2MHz')
% disp('考虑匹配单元中电阻')
% Z_sCsL2=258.71+241.7405281i;
% Z_pCpL1=0.296925-5.103180871i;
% Z_Cs.cm=316.345-38.56459899i;
% Z_Cp.cm=0.186725-5.077311934i;
% Z_L2.cm=Z_sCsL2-Z_Cs.cm;
% Z_L1.cm=1/(1/Z_pCpL1-1/Z_Cp.cm);
% Z_Cs=get_impedance('Z',Z_Cs,f)
% Z_Cp=get_impedance('Z',Z_Cp,f)
% Z_L1=get_impedance('Z',Z_L1,f)
% Z_L2=get_impedance('Z',Z_L2,f)
% disp('忽略匹配单元中电阻')
% Z_sCsL2=imag(Z_sCsL2)*1i;
% Z_pCpL1=imag(Z_pCpL1)*1i;
% Z_Cs.cm=imag(Z_Cs.cm)*1i;
% Z_Cp.cm=imag(Z_Cp.cm)*1i;
% Z_L2.cm=Z_sCsL2-Z_Cs.cm;
% Z_L1.cm=1/(1/Z_pCpL1-1/Z_Cp.cm);
% Z_Cs=get_impedance('Z',Z_Cs,f)
% Z_Cp=get_impedance('Z',Z_Cp,f)
% Z_L1=get_impedance('Z',Z_L1,f)
% Z_L2=get_impedance('Z',Z_L2,f)

% 2021周磊 阻分 直接测量 RF变 应为1MHz
t_ZL.L1 = 92.28e-6;             %变压器原边电感
t_ZL.L2 = 10.38e-6;              %变压器副边线圈电感
t_ZL.k =0.9826 ;                     %互感耦合系数  
t_ZL.M = t_ZL.k*sqrt(t_ZL.L1*t_ZL.L2);              %互感

%% 1MHz 简化匹配电路模型
% 考虑引线电感，忽略RM
fprintf('---------- 1MHz\n')
f=1e6;
%--------------------- 初始化
% 2022陈鹏 测量
ZD0.R=0.382;
ZD0.L=3.16e-6;
Zs0=28.858 +21.987i;
Cs0=8.102e-9; 
Cp0=16.697e-9;

% 2022陈鹏 阻分 间接测量 RF变
t_CP.L1 = 45.307e-6;             %变压器原边电感
t_CP.L2 = 10.646e-6;              %变压器副边线圈电感
t_CP.k =t_ZL.k  ;                     % 未测，使用周磊结果
t_CP.M = t_CP.k*sqrt(t_CP.L1*t_CP.L2);              %互感

c_ZL=c_RF_circuit(f);
% init matching t
c_ZL=c_ZL.init_matching_t('m',t_ZL.M,t_ZL.L1,t_ZL.L2);
% init matching C
c_ZL=c_ZL.init_matching_C( 'gamma', Cs0, Cp0);
% init Zs
c_ZL.s.Z.cm=Zs0;

c_CP=c_ZL.init_matching_t('m',t_CP.M,t_CP.L1,t_CP.L2);
%--------------------- get ZD
fprintf('ZD(阻分) = %.2f Ω+ %.2f μH\n',ZD0.R,ZD0.L*1e6)
c_ZL=c_ZL.get_ZD();
fprintf('ZD(t_ZL) = %.2f Ω+ %.2f μH\n',c_ZL.D.Z.R,c_ZL.D.Z.L*1e6)
c_CP=c_CP.get_ZD();
fprintf('ZD(t_CP) = %.2f Ω+ %.2f μH\n',c_CP.D.Z.R,c_CP.D.Z.L*1e6)

%--------------------- get ZD considering R1 R2
ZD0=get_impedance('RL',ZD0,f);
fprintf('ZD(阻分) = %.2f + %.2f Ω\n',ZD0.R,ZD0.X)
Rline0=0.03011;

c_CP.s.Z.cm=Zs0;
R1=0.025;
R2=0.025;
c_new=get_ZD_R12( c_CP, R1, R2 );
fprintf('ZD(t_CP, R1R2) = %.2f + %.2f Ω\n',c_new.D.Z.R,c_new.D.Z.X)
fprintf('Rline = %.4f Ω\n',c_new.D.Rline)
% 与MMA计算结果一致

R2_list=0:0.005:0.03;
R1_list=flip(5*R2_list);
len=length(R2_list);
Rline=zeros(len,1);
ZD=zeros(len,1);
for i=1:len
    c_new=get_ZD_R12( c_new, R1_list(i), R2_list(i) );
    ZD(i)=c_new.D.Z.cm;
    Rline(i)=c_new.D.Rline;
end

h_fig=figure;
yyaxis left
plot(R2_list, real(ZD),'-.b')
hold on
plot(R2_list, Rline,'-.m')
ylabel('R [\Omega]')
axis([0,max(R2_list),0,0.5]) 
yyaxis right
plot(R2_list, imag(ZD)/(2*pi),'-.c')
ylabel('L [\muH]')
axis([0,max(R2_list),0,3.2]) 
xlabel('R2 [\Omega]')
grid on
L1=legend('RD','Rline','LD');
set(L1,'location','best');
set(L1,'box','off')
set(L1,'AutoUpdate','off')
yyaxis left
line([0,max(R2_list)],[ZD0.R,ZD0.R],'linestyle','--','color','b','LineWidth',1);
line([0,max(R2_list)],[Rline0,Rline0],'linestyle','--','color','m','LineWidth',1);
yyaxis right
line([0,max(R2_list)],1e6*[ZD0.L,ZD0.L],'linestyle','--','color','c','LineWidth',1);

%--------------------- get Zs and compare with Multisim
c_CP=c_CP.init_ZD('RL',1,ZD0);
c_CP=c_CP.get_Zs();
fprintf('Zs = %.2f +∠ %.2f \n',c_CP.s.Z.abs,c_CP.s.Z.phase)
fprintf('Zs = %.2f + %.2f Ω\n',c_CP.s.Z.R,c_CP.s.Z.X)
fprintf('ZD = %.2f +∠ %.2f \n',c_CP.D.Z.abs,c_CP.D.Z.phase)


%% 2MHz 简化匹配电路模型
% 考虑引线电感，忽略RM
fprintf('---------- 2MHz\n')
%--------------------- 初始化
f=2e6;
% 2022陈鹏 测量
ZD0.R=0.296;
ZD0.L=3.09e-6;
ZD0=get_impedance('RL',ZD0,f);
Cs0=2.063e-9; 
Cp0=15.673e-9;
Zs=34.010 - 25.258i;

% 2022陈鹏 阻分 间接测量 RF变
t_CP.L1 = 79.705e-6;             %变压器原边电感
t_CP.L2 = 22.306e-6;              %变压器副边线圈电感
t_CP.k =t_ZL.k  ;                     % 未测，使用周磊结果
t_CP.M = t_CP.k*sqrt(t_CP.L1*t_CP.L2);              %互感

c_ZL=c_RF_circuit(f);
% init matching t
c_ZL=c_ZL.init_matching_t('m',t_ZL.M,t_ZL.L1,t_ZL.L2);
% init matching C
c_ZL=c_ZL.init_matching_C( 'gamma', Cs0, Cp0);
% init Zs
c_ZL.s.Z.cm=Zs;

c_CP=c_ZL.init_matching_t('m',t_CP.M,t_CP.L1,t_CP.L2);
%--------------------- get ZD
fprintf('ZD(阻分) = %.2f Ω+ %.2f μH\n',ZD0.R,ZD0.L*1e6)
c_ZL=c_ZL.get_ZD();
fprintf('ZD(t_ZL) = %.2f Ω+ %.2f μH\n',c_ZL.D.Z.R,c_ZL.D.Z.L*1e6)
c_CP=c_CP.get_ZD();
fprintf('ZD(t_CP) = %.2f Ω+ %.2f μH\n',c_CP.D.Z.R,c_CP.D.Z.L*1e6)

%--------------------- get Zs and compare with Multisim
c_CP=c_CP.init_ZD('RL',1,ZD0);
c_CP=c_CP.get_Zs();
fprintf('Zs = %.2f +∠ %.2f \n',c_CP.s.Z.abs,c_CP.s.Z.phase)
fprintf('Zs = %.2f + %.2f Ω\n',c_CP.s.Z.R,c_CP.s.Z.X)
fprintf('ZD = %.2f +∠ %.2f \n',c_CP.D.Z.abs,c_CP.D.Z.phase)

%--------------------- get ZD 假设k
fprintf('ZD(阻分) = %.2f Ω+ %.2f μH\n',ZD0.R,ZD0.L*1e6)
t_CP.k =1.1  ; 
t_CP.M = t_CP.k*sqrt(t_CP.L1*t_CP.L2);
c_CP=c_CP.init_matching_t('m',t_CP.M,t_CP.L1,t_CP.L2);
c_CP.s.Z.cm=Zs;
c_CP=c_CP.get_ZD();
fprintf('ZD(设k=%.3f) = %.2f Ω+ %.2f μH\n',t_CP.k,c_CP.D.Z.R,c_CP.D.Z.L*1e6)
% 20220509 详见MMA-匹配公式推导.nb，可知2MHz下ZM推得ZD偏小，主要不是
% 因为测得k偏小，而是因为其他原因

%% 结果
%{
---------- 1MHz
ZD(阻分) = 0.38 Ω+ 3.16 μH
ZD(t_ZL) = 0.16 Ω+ 2.91 μH
ZD(t_CP) = 0.34 Ω+ 3.05 μH
---------- 2MHz
ZD(阻分) = 0.30 Ω+ 3.09 μH
ZD(t_ZL) = 0.06 Ω+ 2.76 μH
ZD(t_CP) = 0.15 Ω+ 2.42 μH

分析：
1MHz时ZD(阻分) ≈ZD(t_CP) ，说明使用正确测量参数时，阻分直接测量ZD与 使用电路模型 推算ZD一致
2022陈鹏 阻分 间接测量 RF变 效果 优于 2021周磊 阻分 直接测量 RF变。可能与布线影响有关。
2021周磊 阻分 直接测量 RF变 得到的耦合系数k适用于1MHz情况，不适用于2MHz情况
%}

%% aid function
function obj_c=get_ZD_R12( obj_c, R1, R2 )
% 已知Zs与匹配电路模型，计算ZD
Z_sCp=-get_impedance('parallel',[obj_c.m.p.Z.cm,-obj_c.s.Z.cm]);

ZL1=1i*obj_c.w*obj_c.m.t.L1;
ZL2=1i*obj_c.w*obj_c.m.t.L2;
ZM=1i*obj_c.w*obj_c.m.t.M;
Z2=ZM^2/(ZL1+R1-Z_sCp)-ZL2-R2;

obj_c.D.Z.cm=Z2 - obj_c.m.s.Z.cm;
obj_c.D.Z=get_impedance( 'Z', obj_c.D.Z, obj_c.f );

obj_c.D.Rline=(4.82073e7 - abs((8763.01 -7939.95i) - (28.858 + 31.5192i)*R1)...
    *abs((579.956 +6711.57i) + (2108.43 - 1930.42i)*R1...
    + (8763.01 - 7939.95i)*R2 - (28.858 + 31.5192i)*R1*R2)...
    *cos(phase((0. - 66.8936i) + ( 18384.9 + 7.10543e-14i)...
    /((1.4358 - 276.707i) - R1) -R2)))/abs((8763.01 - 7939.95i) - (28.858 + 31.5192i)*R1)^2;
end