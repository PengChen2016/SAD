%% Main function to generate tests
function tests = test_stoc_model
% test stoc_model
tests = functiontests(localfunctions);
end

% test result: 20210407 @pengchen2016: 
% Three stoc model is ok, according to test_Cazzador_fit /
% test_reproduce_2018Jainb / test_compare_with_v1.

%% Test Functions
function test_Cazzador_fit(testCase)
% test Cazzador_fit
% basic
flag=get_example_flag(0);
input=get_input_data( flag );
flag.stoc_model='Cazzador-fit';
plasma=stochastic_heating_model(flag, input.plasma);
verifyEqual(testCase,size(plasma.nu_st),size(plasma.ne))

% formula
constants=get_constants();
f_Cazzador=2.1e6; %Cazzador拟合式使用的驱动频率，单位Hz
w_RF_Cazzador=2*pi*f_Cazzador;
nu_st_fun_Cazzador=@(x)10.^(-244.1+48.1*log10(x)-3.467*log10(x).^2+...
    0.1113*log10(x).^3-0.001336*log10(x).^4); % 2.1MHz下的nu_st-x关系，x=neTe
% 驱动频率对Cazzador fit关系的影响
nu_st_wrong=nu_st_fun_Cazzador(plasma.ne*plasma.Te); %有较大差别

% 与频率无关的va-X关系
K0=constants.me*constants.me/(constants.mu0*constants.e^3); %Te取eV，不需要再乘以e
va_fun_Cazzador=@(X)nu_st_fun_Cazzador(K0*w_RF_Cazzador^2*X)...
    /w_RF_Cazzador; 
nu_st_fun_new_f2=@(ne, Te, w_RF)va_fun_Cazzador(ne*Te/K0/w_RF^2)*w_RF;
nu_st_right=nu_st_fun_new_f2(plasma.ne, plasma.Te, plasma.w_RF);

verifyEqual(testCase,nu_st_right,plasma.nu_st)
end
% test result: ok. The transformation of Cazzador-fit expression is
% correct.

function test_Vahedi_simplify(testCase)
% basic
flag=get_example_flag(0);
input=get_input_data( flag );
flag.stoc_model='Vahedi-simplify';
plasma=stochastic_heating_model(flag, input.plasma);
verifyEqual(testCase,size(plasma.nu_st),size(plasma.ne))
end
% test result: ok. Code for this case run.

function test_compare_with_v1(testCase)
% compare with fig S1 for code v1 paper v1
flag=get_example_flag(2);
input=get_input_data( flag );
flag.stoc_model='Cazzador-fit';
plasma_Cazzador_fit=stochastic_heating_model(...
    flag, input.plasma);
flag.stoc_model='Vahedi-simplify';
plasma_Vahedi_simplify=stochastic_heating_model(...
    flag, input.plasma);

ne=input.plasma.ne(:,1);
% compare nu_st and delta_st
nu_st1_15eV=plasma_Cazzador_fit.nu_st(:,3);
nu_st2_15eV=plasma_Vahedi_simplify.nu_st(:,3);
delta_st1_15eV=plasma_Cazzador_fit.delta_st(:,3);
delta_st2_15eV=plasma_Vahedi_simplify.delta_st(:,3);

% plot
figure
yyaxis left
loglog(ne,nu_st1_15eV,'--c');
axis([1e16,1e19,4e6,1e8])
yyaxis right
loglog(ne,delta_st1_15eV,'--m');
axis([1e16,1e19,1e-3,3e-1])
hold on % make the figure fixxed
yyaxis left
loglog(ne,nu_st2_15eV,'-.c');
ylabel('{\it\nu}_{st}');
yyaxis right
loglog(ne,delta_st2_15eV,'-.m');
ylabel('{\it\delta}_{st}');

xlabel('{\itn}_e')
grid on
end
% test result: ok.
% details: .\others\Figures during code developing.pptx

function test_reproduce_2018Jainb(testCase)
% reproduce nu_st in the fig40 of 2018Jainb - Studies and experimental
% activities to qualify the behaviour of RF power circuits for Negative Ion
% Sources of Neutral Beam Injectors for ITER and fusion experiments
flag.input_plasma='2018Jainb_ELISE_sweep_f';
input=get_input_data( flag );
plasma= Ohmic_heating_model( input.plasma, 'e-H2-Phelps' );

% plasma1=plasma;
% constants=get_constants();
% alpha_st_fun=@(delta_st,Ve,w_RF) 4*delta_st.^2.*w_RF.^2/pi./(Ve.^2);
% % pre-processing
% if 1==plasma1.size
%     size_mat=[1,1];
% else
%     size_mat=plasma1.size(plasma1.size>1);
% end
% plasma1.alpha_st=zeros(size_mat);
% plasma1.delta_st=zeros(size_mat);
% plasma1.nu_st=zeros(size_mat);
% % For vectorization parallel
% temp_ve=sqrt(plasma1.Te*constants.e/constants.me); %2018Jainb
% if isscalar(plasma1.Te)
%     temp_ve=repmat(temp_ve,size_mat);
% end
% if isscalar(plasma1.ne)
%     temp_wpe=repmat(plasma1.wpe,size_mat);
% else
%     temp_wpe=plasma1.wpe;
% end
% if isscalar(plasma1.w_RF)
%     temp_f=repmat(plasma1.f,size_mat);
%     temp_w_RF=repmat(plasma1.w_RF,size_mat);
% else
%     temp_f=plasma1.f;
%     temp_w_RF=plasma1.w_RF;
% end
% 
% if ~isfield(plasma1,'nu_m') || isempty(plasma1.nu_m)
%     warning('nu_m is needed.')
%     plasma1=Ohmic_heating_model(plasma1, 'Phelps-m');
% end
% 
% plasma1.delta_st=get_plasma_skin_depth('as-medium-simplified-finite-radius',...
%     temp_f,plasma1.nu_m,temp_wpe,plasma1.r);
% plasma1.alpha_st=alpha_st_fun(plasma1.delta_st,temp_ve,temp_w_RF);
% %         % case 1: Vahedi-simplify 3 phase. Not used.
% %         % α <= 0.03
% %         idx1=plasma.alpha_st<=0.03;
% %         plasma.nu_st(idx1)=temp_ve(idx1)./(2*pi*plasma.delta_st(idx1));
% %         % 0.03 < α <= 10
% %         idx2=(plasma.alpha_st>0.03) & (plasma.alpha_st<=10);
% %         plasma.nu_st(idx2)=temp_ve(idx2)./plasma.delta_st(idx2)/4;
% %         % α > 10
% %         idx3=plasma.alpha_st>10;
% %         plasma.nu_st(idx3)=pi/4./(temp_w_RF(idx3).^2)...
% %             .*(temp_ve(idx3)./plasma.delta_st(idx3)).^3;
% %         % case 2: Vahedi-simplify 2 phase. Not used.
% %         % α << 1
% %         idx1=plasma.alpha_st<=1;
% %         plasma.nu_st(idx1)=temp_ve(idx1)./(2*pi*plasma.delta_st(idx1));
% %         % α >> 1
% %         idx2=plasma.alpha_st>1;
% %         plasma.nu_st(idx2)=pi/4./(temp_w_RF(idx2).^2)...
% %             .*(temp_ve(idx2)./plasma.delta_st(idx2)).^3;
% % case 3: Cazzador-simplify 2 phase. Not used.
% x=plasma1.ne.*plasma1.Te;
% % α << 1
% idx1=plasma1.alpha_st<=1;
% A=0.47; % 2014Cazzador fit
% B=-0.18; % 2014Cazzador Result in small difference.
% %         B=0.18; % 2018Jain Result in big difference.
% plasma1.nu_st(idx1)=((temp_w_RF(idx1)/A/sqrt(4*pi)).*...
%     (constants.mu0*constants.e^3*x(idx1)./temp_w_RF(idx1)/constants.me^2).^(B+1/2))...
%     .^(1/(B+3/2));
% % α >> 1
% idx2=plasma1.alpha_st>1;
% % 2018Jain ok
% plasma1.nu_st(idx2)=pi/4./(temp_w_RF(idx2).^2).*...
%     (8*constants.mu0*constants.e^3*x(idx2)/pi/constants.me^2).^(3/2);
% % 2004Cazzador Result in big difference.
% %         plasma.nu_st(idx2)=pi/4./(temp_w_RF(idx2).^2).*...
% %             (8*constants.mu0*constants.e^3*x(idx2).*temp_w_RF(idx2)/pi/constants.me^2).^(3/2);
% case 4: Vahedi/Cazzador-simplify 2 phase. Used.
flag.stoc_model='2018Jainb-simplify';
plasma2=stochastic_heating_model(flag, plasma);
% plot
figure
loglog(plasma2.f(:,1),plasma2.nu_st(:,1),'-.c');
hold on
loglog(plasma2.f(:,1),plasma2.nu_st(:,2),'-.m');
loglog(plasma2.f(:,1),plasma2.nu_st(:,3),'-.y');
loglog(plasma2.f(:,1),plasma2.nu_st(:,4),'-.b');

xlabel('{\itf}[Hz]');
ylabel('collision frequency');
grid on%显示网格
axis([1e2,1e10,1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
end
% test result: ok.
% details: .\others\Figures during code developing.pptx

function test_compare_models(testCase)
% test compare different stoc models
flag=get_example_flag(2);
input=get_input_data( flag );
flag.stoc_model='Cazzador-fit';
plasma_Cazzador_fit=stochastic_heating_model(...
    flag, input.plasma);
flag.stoc_model='Vahedi-simplify';
plasma_Vahedi_simplify=stochastic_heating_model(...
    flag, input.plasma);
% flag.stoc_model='Cazzador-simplify';

ne=input.plasma.ne(:,1);
% compare nu_st
nu_st1_5eV=plasma_Cazzador_fit.nu_st(:,1);
nu_st1_15eV=plasma_Cazzador_fit.nu_st(:,3);
nu_st1_25eV=plasma_Cazzador_fit.nu_st(:,5);
nu_st2_5eV=plasma_Vahedi_simplify.nu_st(:,1);
nu_st2_15eV=plasma_Vahedi_simplify.nu_st(:,3);
nu_st2_25eV=plasma_Vahedi_simplify.nu_st(:,5);

% plot
figure
h_plot1=semilogx(ne,nu_st1_5eV,'.c','LineWidth',6);
hold on
h_plot2=semilogx(ne,nu_st1_15eV,'.m','LineWidth',6);
h_plot3=semilogx(ne,nu_st1_25eV,'.k','LineWidth',6);
h_plot4=semilogx(ne,nu_st2_5eV,'-.c');
semilogx(ne,nu_st2_15eV,'-.m')
semilogx(ne,nu_st2_25eV,'-.k')
ylabel('{\it\nu}_{st}');
xlabel('{\itn}_e')
grid on
L1=legend([h_plot1,h_plot2,h_plot3],{'Te=5eV, Cazzador fit','Te=15eV','Te=25eV'});
set(L1,'location','best')
set(L1,'AutoUpdate','off')
% 额外图例：第二个坐标轴中处理
axes2 = axes('position',get(gca,'position'),'visible','off');
L2=legend(axes2, [h_plot1,h_plot4],{'Te=5eV, Cazzador fit','Vahedi simplify'});
set(L2,'location','best')
set(L2,'AutoUpdate','off')

% compare alpha_st
alpha_st1_5eV=plasma_Cazzador_fit.alpha_st(:,1);
alpha_st1_15eV=plasma_Cazzador_fit.alpha_st(:,3);
alpha_st1_25eV=plasma_Cazzador_fit.alpha_st(:,5);
alpha_st2_5eV=plasma_Vahedi_simplify.alpha_st(:,1);
alpha_st2_15eV=plasma_Vahedi_simplify.alpha_st(:,3);
alpha_st2_25eV=plasma_Vahedi_simplify.alpha_st(:,5);

% plot
figure
h_plot1=semilogx(ne,alpha_st1_5eV,'--c');
hold on
h_plot2=semilogx(ne,alpha_st1_15eV,'--m');
h_plot3=semilogx(ne,alpha_st1_25eV,'--k');
h_plot4=semilogx(ne,alpha_st2_5eV,'-.c');
semilogx(ne,alpha_st2_15eV,'-.m')
semilogx(ne,alpha_st2_25eV,'-.k')
ylabel('{\it\alpha}_{st}');
xlabel('{\itn}_e')
grid on
L1=legend([h_plot1,h_plot2,h_plot3],{'Te=5eV, Cazzador fit','Te=15eV','Te=25eV'});
set(L1,'location','best')
set(L1,'AutoUpdate','off')
% 额外图例：第二个坐标轴中处理
axes2 = axes('position',get(gca,'position'),'visible','off');
L2=legend(axes2, [h_plot1,h_plot4],{'Te=5eV, Cazzador fit','Vahedi simplify'});
set(L2,'location','best')
set(L2,'AutoUpdate','off')

% compare delta_st
delta_st1_5eV=plasma_Cazzador_fit.delta_st(:,1);
delta_st1_15eV=plasma_Cazzador_fit.delta_st(:,3);
delta_st1_25eV=plasma_Cazzador_fit.delta_st(:,5);
delta_st2_5eV=plasma_Vahedi_simplify.delta_st(:,1);
delta_st2_15eV=plasma_Vahedi_simplify.delta_st(:,3);
delta_st2_25eV=plasma_Vahedi_simplify.delta_st(:,5);

% plot
figure
h_plot1=semilogx(ne,delta_st1_5eV,'--c');
hold on
h_plot2=semilogx(ne,delta_st1_15eV,'--m');
h_plot3=semilogx(ne,delta_st1_25eV,'--k');
h_plot4=semilogx(ne,delta_st2_5eV,'-.c');
semilogx(ne,delta_st2_15eV,'-.m')
semilogx(ne,delta_st2_25eV,'-.k')
ylabel('{\it\delta}_{st}');
xlabel('{\itn}_e')
grid on
L1=legend([h_plot1,h_plot2,h_plot3],{'Te=5eV, Cazzador fit','Te=15eV','Te=25eV'});
set(L1,'location','best')
set(L1,'AutoUpdate','off')
% 额外图例：第二个坐标轴中处理
axes2 = axes('position',get(gca,'position'),'visible','off');
L2=legend(axes2, [h_plot1,h_plot4],{'Te=5eV, Cazzador fit','Vahedi simplify'});
set(L2,'location','best')
set(L2,'AutoUpdate','off')

% 对比stoc模型参数
x1='ne';
x2='Te';
y1='alpha_st';
y2='delta_st';

X1=plasma_Cazzador_fit.(x1);
Y1_C=plasma_Cazzador_fit.(y1);
Y2_C=plasma_Cazzador_fit.(y2);
Y1_V=plasma_Vahedi_simplify.(y1);
Y2_V=plasma_Vahedi_simplify.(y2);
% 可视化
h1=figure;
yyaxis left;
X2=1;h_plot1=semilogx(X1(:,X2),Y1_C(:,X2),'-r');hold on;h_plot4=semilogx(X1(:,X2),Y1_V(:,X2),'--r');
X2=3;h_plot2=semilogx(X1(:,X2),Y1_C(:,X2),'-g');semilogx(X1(:,X2),Y1_V(:,X2),'--g');
X2=5;h_plot3=semilogx(X1(:,X2),Y1_C(:,X2),'-b');semilogx(X1(:,X2),Y1_V(:,X2),'--b');
ylabel('\alpha_{st}');
yyaxis right;
X2=1;h_plot5=semilogx(X1(:,X2),Y2_C(:,X2),'-or');hold on;h_plot6=semilogx(X1(:,X2),Y2_V(:,X2),'--or');
X2=3;semilogx(X1(:,X2),Y2_C(:,X2),'-og');semilogx(X1(:,X2),Y2_V(:,X2),'--og');
X2=5;semilogx(X1(:,X2),Y2_C(:,X2),'-ob');semilogx(X1(:,X2),Y2_V(:,X2),'--ob');
ylabel('\delta_{st}');
xlabel(x1)
grid on
yyaxis left;
line([X1(1),X1(end)],[0.03,0.03],'linestyle','-.','color','k');
% 图例
L1=legend([h_plot1,h_plot2,h_plot3],{'Te=5eV, Cazzador fit, \alpha_{st}','Te=15eV','Te=25eV'});
set(L1,'location','best')
set(L1,'AutoUpdate','off')
% 额外图例：第二个坐标轴中处理
axes2 = axes('position',get(gca,'position'),'visible','off');
L2=legend(axes2, [h_plot1,h_plot4],{'Te=5eV, Cazzador fit, \alpha_{st}','Vahedi simplify'});
set(L2,'location','best')
set(L2,'AutoUpdate','off')
axes3 = axes('position',get(gca,'position'),'visible','off');
L3=legend(axes3, [h_plot1,h_plot5],{'Te=5eV, Cazzador fit, \alpha_{st}','\delta_{st}'});
set(L3,'location','best')
set(L3,'AutoUpdate','off')
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