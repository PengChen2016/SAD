%% Main function to generate tests
function tests = test_stoc_model
% test stoc_model
tests = functiontests(localfunctions);
end

%% Test Functions
function test_Cazzador_fit(testCase)
% test Cazzador_fit
% basic
flag=get_example_flag(0);
input=get_input_data( flag );
plasma=stochastic_heating_model(flag.stoc_model, input.plasma);
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
% test passed: 变换前后等价
end

function test_Vahedi_simplify(testCase)
flag=get_example_flag(0);
input=get_input_data( flag );
plasma=stochastic_heating_model('Vahedi-simplify', input.plasma);
verifyEqual(testCase,size(plasma.nu_st),size(plasma.ne))
end

function test_compare_models(testCase)
% test compare different stoc models
flag=get_example_flag(2);
input=get_input_data( flag );
plasma_Cazzador_fit=stochastic_heating_model(...
    'Cazzador-fit', input.plasma);
plasma_Vahedi_simplify=stochastic_heating_model(...
    'Vahedi-simplify', input.plasma);
% flag.stoc_model='Cazzador-simplify';

verifyEqual(testCase,size(plasma_Cazzador_fit.nu_st),...
    size(plasma_Cazzador_fit.ne))
verifyEqual(testCase,size(plasma_Vahedi_simplify.nu_st),...
    size(plasma_Vahedi_simplify.ne))

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
L1=legend([h_plot1,h_plot2,h_plot3],{'Te=1eV, Cazzador fit, \alpha_{st}','Te=3eV','Te=5eV'});
set(L1,'location','best')
set(L1,'AutoUpdate','off')
% 额外图例：第二个坐标轴中处理
axes2 = axes('position',get(gca,'position'),'visible','off');
L2=legend(axes2, [h_plot1,h_plot4],{'Te=1eV, Cazzador fit, \alpha_{st}','Vahedi simplify'});
set(L2,'location','best')
set(L2,'AutoUpdate','off')
axes3 = axes('position',get(gca,'position'),'visible','off');
L3=legend(axes3, [h_plot1,h_plot5],{'Te=1eV, Cazzador fit, \alpha_{st}','\delta_{st}'});
set(L3,'location','best')
set(L3,'AutoUpdate','off')

% vst简化表达式与Cazzador fit表达式对比
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