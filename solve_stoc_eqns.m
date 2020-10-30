%% Abstract
% solve stoc eqns
%% History
% v190111 对Φ(α)的拟合自左晨，赵鹏
% v200424 基于2014Cazzador，

%% TODO

%% Introduction
% 见 1995Vahedia的ln

close all
clear
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
disp(now_str)
%物理常量
e=1.6e-19;            %电子电量
me=9.1e-31;           %电子质量
mi=1.67e-27;          %氢离子质量
eps0=8.85e-12;          %真空介电常数
c=3e8;                %真空光速[m/s]
kB=1.381e-23;         %玻尔兹曼常数 J/K
mu0=4*pi*1e-7;         %真空磁导率
sigma_Cu=5.8e7;       %铜的电导率

%% 输入
dne=15:1:19;             %指数
ne=10.^dne;                 %电子密度[m^-3]
Te=5:5:25;                 %电子温度[eV]
num_ne=length(ne);
num_Te=length(Te);
f=10^6;               %驱动频率，单位Hz
w_RF=2*pi*f;

% Φ(α)
ei_Vahedi=@(alpha)-ei(-alpha);
%     %检验表达式
%     syms x alpha;
%     ei_Vahedi(i)=double(int(exp(-x)/x,x,alpha,inf));
jphi_fun=@(alpha)(exp(alpha)*(1+alpha)*ei_Vahedi(alpha)-1)/pi; % vstoc与α之间的关系函数

%% 拟合Φ(α)
da=-4:0.2:2; %指数
alpha_arr=10.^da; %表征电子穿过等离子体趋肤层耗费时间与射频电磁场周期之比值的参量，其值越小，随机加热更加显著。
num_alpha=length(alpha_arr);
for ai=1:num_alpha
    jphi_origin(ai)=jphi_fun(alpha_arr(ai));
    %分段拟合
    if alpha_arr(ai)<0.03
       jphi_fit_Vahedi(ai)=pi/2;%1995Vahedia拟合
       jphi_fit1_ZP(ai)=-(log(alpha_arr(ai))+1.58)/pi; %赵鹏/左晨拟合
    else
        if alpha_arr(ai)<10
            jphi_fit_Vahedi(ai)=1/8/pi/alpha_arr(ai);
            jphi_fit1_ZP(ai)=jphi_fit_Vahedi(ai);
        else
            jphi_fit_Vahedi(ai)=1/pi/(alpha_arr(ai)^2);
            jphi_fit1_ZP(ai)=jphi_fit_Vahedi(ai);
        end
    end

    %赵鹏/左晨拟合    
    if alpha_arr(ai)<0.1 %第一段与第二段分段界限与fit1不同
        jphi_fit2_ZP(ai)=0.4*alpha_arr(ai)^(-0.2); %赵鹏/左晨拟合
    else
        if alpha_arr(ai)<10
            jphi_fit2_ZP(ai)=1/8/pi/alpha_arr(ai); %第二段表达式与fit1相同
        else
            jphi_fit2_ZP(ai)=1/pi/(alpha_arr(ai)^2);
        end
    end    
    
    %2014Cazzador
    if alpha_arr(ai)<1
        jphi_fit_Cazzador(ai)=0.47*alpha_arr(ai)^(-0.18);
    else
        jphi_fit_Cazzador(ai)=1/pi/(alpha_arr(ai)^2);
    end
end

% loglog(alpha_arr,jphi_origin,...
%     alpha_arr,jphi_fit_Vahedi,alpha_arr,jphi_fit1_ZP,alpha_arr,jphi_fit2_ZP,...
%     alpha_arr,jphi_fit_Cazzador,'-','LineWidth',1.5);
% L1=legend('\itI\rm(\alpha)','Vahedi fit','ZP fit1','ZP fit2','Cazzador fit','Location','southWest');
loglog(alpha_arr,jphi_origin,...
    alpha_arr,jphi_fit_Vahedi,'-','LineWidth',1.5);
L1=legend('\itI\rm(\alpha)','Vahedi fit','Location','southWest');
set(L1,'FontSize',12);
xlabel('\alpha');
ylabel('\itI\rm(\alpha)');
title('vstoc与α之间的关系函数曲线Φ(α)的拟合')
%结果见result_fit_Phi_alpha_200424.jpg，个人喜欢ZP fit2，无波折且贴近

%% 数值求解stoc方程组
% eqn1
% 隐函数绘图失败
% f_eqn1=@(alpha,va)2*sqrt(pi*alpha)*jphi_fun(alpha)-va/(1+va*va);
% fimplicit(f_eqn1,[0 0.5 0 5])
T_eqn1=@(alpha)2*sqrt(pi*alpha)*jphi_fun(alpha);
vaH=@(alpha)(1+sqrt(1-4*T_eqn1(alpha)*T_eqn1(alpha)))/2/T_eqn1(alpha);
vaL=@(alpha)(1-sqrt(1-4*T_eqn1(alpha)*T_eqn1(alpha)))/2/T_eqn1(alpha);
alpha_eqn1=[0:0.0005:0.003 (0.003+0.001):0.001:0.01 0.02:0.02:0.7]; % 根据vahedi模型算出来的alpha范围
num_alpha=length(alpha_eqn1);
for ai=1:num_alpha
    vaH_eqn1(ai)=vaH(alpha_eqn1(ai));
    vaL_eqn1(ai)=vaL(alpha_eqn1(ai));
end
figure
plot(alpha_eqn1,vaH_eqn1,'-b','LineWidth',1);
hold on
plot(alpha_eqn1,vaL_eqn1,'-b','LineWidth',1);
hold on
xlabel('\it\alpha');
ylabel('\it\nu_{\rma}=\it\nu_{\rmst}/\omega');
axis([0,0.7,0,7]) %绘图显示范围，即[xmin,xmax,ymin,ymax]

% eqn2
K0=me*me/(mu0*e^3); %与频率无关 %Te取eV，不需要再乘以e
% % 判断X范围
% for nei=1:num_ne
%     for Tei=1:num_Te
%         x(nei,Tei)=ne(nei)*Te(Tei);
%         X(nei,Tei)=x(nei,Tei)/K0/w_RF^2;
%     end
% end
% disp(X)

dX=(0:1:20)/5; %指数
X=10.^dX;
num_X=length(X);
for Xi=1:num_X
    alpha_va_fun2=@(va)(1+va*va)/(1+sqrt(1+va*va))/X(Xi);
    va_eqn2=0:0.05:7;
    num_va=length(va_eqn2);
    for vai=1:num_va
        alpha_eqn2(Xi,vai)=alpha_va_fun2(va_eqn2(vai));
    end
    plot(alpha_eqn2(Xi,:),va_eqn2,'-r','LineWidth',1);
    hold on
%     [x,y]=ginput(2);%光标移到交点处点击专鼠标，获得属两个交点坐标
%     str=sprintf("(%f,%f)",x,y);
%     text(x,y,str); %在图上标记
end

L1=legend('eqn1-H','eqn2-L','eqn2 at X_i(n_e,T_e,\omega)');
set(L1,'FontSize',10);
title(['numerical root of stoc eqns' ' \rmat \rm' now_str]);
grid on%显示网格

%% 拟合va-X
% xlswrite('Result\vst-alpha,X数值根.xlsx',X',1,'A3')
% 从excel表导入数据
% rootH=vstalphaX
% rootL=vstalphaX
% X_rootH=table2array(rootH(:,1));
% va_rootH=table2array(rootH(:,3));
% X_rootL=table2array(rootL(:,1));
% va_rootL=table2array(rootL(:,3));

load('result/rootH_L.mat')

figure
loglog(X_rootH,va_rootH,'s');
hold on
loglog(X_rootL,va_rootL,'o');
hold on
xlabel('\itX=\rm(\itv_{\rme}\itT_{\rme}\rm)/\rm(K_0\it\omega^{\rm2}\rm)');
ylabel('\it\nu_{\rma}=\nu_{\rmst}/\omega');
% axis([0,0.7,0,7]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
title(['numerical root of stoc eqns' ' \rmat \rm' now_str]);
grid on%显示网格
% 以上关系，是与频率无关的。如果某频率下可以取到该α/X范围的值，则知道va的值
line([1.6e2,1.6e2],[va_rootL(1),va_rootH(1)],'linestyle','-.','linewidth',1.5,'color','k');
line([3.1e3,3.1e3],[va_rootL(1),va_rootH(1)],'linestyle','-.','linewidth',1.5,'color','k');
L1=legend('H','L','1e17,5eV@1MHz','1e18,15eV@1MHz');
set(L1,'location','northwest');

% %2014Cazzador拟合效果检验
% f_Cazzador=2.1e6;               %驱动频率，单位Hz
% w_RF_Cazzador=2*pi*f_Cazzador;
% vst_fit_Cazzador=@(x)10^(-244.1+48.1*log10(x)-3.467*log10(x)^2+0.1113*log10(x)^3-0.001336*log10(x)^4);
% va_fit_Cazzador=@(X)vst_fit_Cazzador(K0*w_RF_Cazzador^2*X)/w_RF_Cazzador;
% % 验证了计算流程：即将频率无关的通用关系，处理成语频率相关的结果。
% for Xi=1:num_X
%     va_arr_fit_Cazzador(Xi)=va_fit_Cazzador(X(Xi));
% end
% loglog(X,va_arr_fit_Cazzador,'-');

% 陈鹏200424拟合
% X_rootHL=[flip(X_rootL(end-3:end));flip(X_rootH(1:end-4))];
% va_rootHL=[flip(va_rootL(end-3:end));flip(va_rootH(1:end-4))];
% loglog(X_rootHL,va_rootHL,'-');
A0=-0.500323937504147;
A1=0.391557506163918;
A2=0.141238121518155;
A3=-0.0891994841475049;
A4=0.0125156938251928;
va_fit_PC=@(X)10^(A0+A1*log10(X)+A2*log10(X)^2+A3*log10(X)^3+A4*log10(X)^4);
for Xi=1:num_X
    va_arr_fit_PC(Xi)=va_fit_PC(X(Xi));
end
loglog(X,va_arr_fit_PC,'-');
L1=legend('H','L','1e17,5eV@1MHz','1e18,15eV@1MHz','Cazzador fit','PC fit');
L1=legend('H','L','1e17,5eV@1MHz','1e18,15eV@1MHz','PC fit');

% %vst-x
% figure
% loglog(X_rootH*K0*w_RF_Cazzador*w_RF_Cazzador,va_rootH*w_RF_Cazzador,'s');
% hold on
% loglog(X_rootL*K0*w_RF_Cazzador*w_RF_Cazzador,va_rootL*w_RF_Cazzador,'o');
% hold on
% xlabel('\itv_{\rme}\itT_{\rme}');
% ylabel('\nu_{\rmst}');
% % axis([0,0.7,0,7]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
% L1=legend('H','L');
% set(L1,'FontSize',10);
% title(['numerical root of stoc eqns' ' \rmat \rm' now_str]);
% grid on%显示网格


