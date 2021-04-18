% fit spatially nonuniform ne/Te of plasma in CHARLIE

%% Abstract
% 对CHARLIE轴向、径向电子密度分布做拟合与拟合表达式的测试
%% History
% 拟合ok v200419 by PengChen
% 绘图用于论文 v200914 by PengChen

%% TODO
% 分段取积分平均值
% 如何在Matlab中做拟合？

%% Introduction
close all
clear
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
% now_str=datestr(now,'yyyy-mm-dd');
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
%%%% 控制位
flag_output_for_paper=true;

%% axial electron density data from 2019Raunner
x_experiment1=[70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200];
y_experiment1=[0.196224, 0.23341, 0.286613, 0.340961, 0.405034, 0.476545, 0.553776, 0.639588, 0.727689, 0.816362, 0.891304, 0.947941, 0.98913, 1.00229];
x_fit1_point=[20:20:60,70:10:200];
x_fit1_line=0:1:200;
x_fit1_rauner=20:20:60;
y_fit1_rauner=[0.0446224,0.0926773,0.15389];

% % fit10
% doing_str='fit10: x0=200, read half_gamma from figure';
% x0=200;
% half_gamma0=76.36; %手测半高宽，可见拟合结果误差较大
% L_max=Lorentzian([x0,half_gamma0],x0);
% k0=max(y_experiment1)/L_max;
% y_fit10_point= k0*Lorentzian([x0,half_gamma0],x_fit1_point);
% y_fit10_line=k0*Lorentzian([x0,half_gamma0],x_fit1_line);
% fprintf('%s: %s = %.2e ,%s = %.2e ,%s = %.2e\n',doing_str,'k',k0,'x0',x0,'half_gamma',half_gamma0);
% 
% % fit11
% doing_str='fit11: x0=200, fit k and half_gamma in Matlab';
% x1=x0;
% func_fit1 = @(coff_vec,x_vec) coff_vec(1)*Lorentzian([x1,coff_vec(2)],x_vec);
% opts = statset('nlinfit');
% opts.RobustWgtFun = 'bisquare';
% coff_vec0=[half_gamma0,k0];
% [coff_vec1,~]= nlinfit(x_experiment1,y_experiment1,func_fit1,coff_vec0,opts);
% %f:符号函数句柄,如果是以m文件的形式调用的时候,需要加@.
% % 拟合参数的标准是(f-y)^2取最小值
% %a:最开始预估的值(预拟合的未知参数的估计值)。
% k1=coff_vec1(1);
% half_gamma1=coff_vec1(2);
% fprintf('%s: %s = %.2e ,%s = %.2e ,%s = %.2e\n',doing_str,'k',k1,'x0',x1,'half_gamma',half_gamma1);
% y_fit11_point= func_fit1(coff_vec1,x_fit1_point);
% y_fit11_line=func_fit1(coff_vec1,x_fit1_line);
% 
% % fit12
% doing_str='fit12: fit k, half_gamma in 1stOpt';
% x2=x0;
% k2=234.98;
% half_gamma2=72.88; 
% fprintf('%s: %s = %.2e ,%s = %.2e ,%s = %.2e\n',doing_str,'k',k2,'x0',x2,'half_gamma',half_gamma2);
% y_fit12_point= k2*Lorentzian([x2,half_gamma2],x_fit1_point);
% y_fit12_line=k2*Lorentzian([x2,half_gamma2],x_fit1_line);
% 
% % fit13
% doing_str='fit13: fit half_gamma in 1stOpt';
% x3=x0;
% k3=k0;
% half_gamma3=74.73; 
% y_fit13_point= k3*Lorentzian([x3,half_gamma3],x_fit1_point);
% y_fit13_line=k3*Lorentzian([x3,half_gamma3],x_fit1_line);
% fprintf('%s: %s = %.2e ,%s = %.2e ,%s = %.2e\n',doing_str,'k',k3,'x0',x3,'half_gamma',half_gamma3);

% fit14
doing_str='fit14: 多项式拟合 in 1stOpt';
y_poly=@(x)-6.2861e-3+4.0769e-3.*x-6.6033e-5.*x.^2+8.4313e-7*x.^3-2.4400e-9*x.^4;
y_fit14_point= y_poly(x_fit1_point);
y_fit14_line=y_poly(x_fit1_line);

intID_y_poly=@(x)-6.2861e-3*x+4.0769e-3.*x^2/2-6.6033e-5.*x.^3/3+8.4313e-7*x.^4/4-2.4400e-9*x.^5/5;
y_poly_mean=@(x_down,x_up) (intID_y_poly(x_up)-intID_y_poly(x_down))/(x_up-x_down);
y_mean=y_poly_mean(0,200);





if ~flag_output_for_paper
%% 后处理
figure
scatter(x_experiment1,y_experiment1,'o')
hold on
% plot(x_fit1_line,y_fit10_line,'r-')
% hold on
% plot(x_fit1_line,y_fit11_line,'b-')
% hold on
% plot(x_fit1_line,y_fit12_line,'y-')
% hold on
% plot(x_fit1_line,y_fit13_line,'k-')
% hold on
plot(x_fit1_line,y_fit14_line,'g-')
hold on
scatter(x_fit1_rauner,y_fit1_rauner,'+')
% 图像格式控制
title(['fit of ne(z) of CHARLIE \rmat \rm' now_str])
xlabel('z')
ylabel('normalized ne')
% legend('experiment','fit10','fit11','fit12','fit13','fit14','Location','best') %图例。位置为最佳，可以不写，有需要再手调。
legend('experiment','CP-fit','Rauner fit','Location','best') %图例。位置为最佳，可以不写，有需要再手调。
% axis([-2*pi,2*pi,-8,2]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
% set(gca,'xTick',-2*pi:pi/2:2*pi) %设置x轴标签点
% set(gca,'xTickLabel',{'-2pi','','-pi','','0','','pi','','2pi'}) %设置x轴标签值
% text(1.5,-2,'y1在原点处不连续','color','r') %在指定点添加文本
line([0,200],[y_mean,y_mean],'linestyle','--');

end


%% lateral intensity data from 2019Raunner
x_experiment2=[45.5, 42, 38.5, 35, 31.5, 28, 24.5, 21, 17.5, 14, 10.5, 7, 3.5, 0];
y_experiment2=[0.00146771, 0.255382, 0.412427, 0.515656, 0.618885, 0.708904, 0.774462, 0.824364, 0.876712, 0.926125, 0.957926, 0.973581, 0.991194, 0.998532];

x_fit2_point=x_experiment2;
x_fit2_line=0:0.5:45.5;
R=45.5;

% fit20
% 归一化Hβ谱线强度 对 轴向z=100mm，侧壁截弦 的平均
doing_str='fit20: 4阶多项式拟合 in 1stOpt';
a0=0.9920;
a1=5.4601e-3;
a2=-1.1959e-3;
a3=3.7981e-5;
a4=-5.4331e-7;
y_poly=@(x)a0+a1.*x+a2.*x.^2+a3*x.^3+a4*x.^4;
y_fit20_point= y_poly(x_fit2_point);
y_fit20_line=y_poly(x_fit2_line);

% Abel逆变换 见MMA Abel反演.nb 
f_ne=@(r)(sqrt(R^2-r^2)*(12*a2+9*a3*R+8*a4*(2*r^2+R^2))+(6*a1+9*a3*r^2)*log((R+sqrt(R^2-r^2))/r))/6;
f_ne(0); %Inf
f_ne(R); %0
f_ne(R/2); % -0.0347
f_ne(2); % -0.0376
% 不可接受

if ~flag_output_for_paper
%% 后处理
figure
scatter(x_experiment2,y_experiment2,'o')
hold on
plot(x_fit2_line,y_fit20_line,'r-')
% 图像格式控制
title(['fit of ne(r) of CHARLIE \rmat \rm' now_str])
xlabel('r')
ylabel('normalized ne')
legend('experiment','CP-fit','Location','best')
% axis([-2*pi,2*pi,-8,2]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
% set(gca,'xTick',-2*pi:pi/2:2*pi) %设置x轴标签点
% set(gca,'xTickLabel',{'-2pi','','-pi','','0','','pi','','2pi'}) %设置x轴标签值
% text(1.5,-2,'y1在原点处不连续','color','r') %在指定点添加文本
% line([500,2500],[1,1],'linestyle','--');
end


%% radial density data from 2017Raunera-fig4
x_experiment3=[45, 40, 35, 30, 25, 20, 15, 10, 5, 0];
y_experiment3=[[0.04, 0.18, 0.39, 0.56, 0.7, 0.8, 0.87, 0.9, 0.94, 1]];

x_fit3_point=x_experiment3;
x_fit3_line=0:0.5:45.5;
R=45.5;

% fit30
doing_str='fit30: Bessel fit in 1stOpt';
% 零阶贝塞尔

a0=0.9920;
a1=5.4601e-3;
a2=-1.1959e-3;
a3=3.7981e-5;
a4=-5.4331e-7;
y_poly=@(x)a0+a1.*x+a2.*x.^2+a3*x.^3+a4*x.^4;
y_fit20_point= y_poly(x_fit2_point);
y_fit20_line=y_poly(x_fit2_line);

% fit31
doing_str='fit31: 4阶多项式拟合 in 1stOpt';
a0=0.9920;
a1=5.4601e-3;
a2=-1.1959e-3;
a3=3.7981e-5;
a4=-5.4331e-7;
y_poly=@(x)a0+a1.*x+a2.*x.^2+a3*x.^3+a4*x.^4;
y_fit20_point= y_poly(x_fit2_point);
y_fit20_line=y_poly(x_fit2_line);



%% 论文绘图
if flag_output_for_paper
plot_line_width=3;
gca_line_width=1;
marker_size=8;
font_size=15;

% 轴向
x_plot=[0, 20, 40, 60, 70, 80, 90, 100, 110, 120, 130, 140, 149.99, 150, 160, 170, 180, 190, 200];
y_nonuniform_plot=[0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9];
figure
h1=scatter(x_axial_competion(x_experiment1),y_axial_competion(y_experiment1),'o','MarkerEdgeColor','r','LineWidth',plot_line_width);
hold on
% scatter(x_fit1_rauner,y_fit1_rauner,'+')
h2=plot(x_axial_competion(x_fit1_line),y_axial_competion(y_fit14_line),'-b','LineWidth',plot_line_width);
h3=plot(x_axial_competion(x_plot),y_axial_competion(y_nonuniform_plot),'--r','LineWidth',plot_line_width);
h4=line([-200,200],[1,1],'Color','k','linestyle','-.','LineWidth',plot_line_width);
% line([0,200],[y_mean,y_mean],'linestyle','--');
% y_uniform_plot=ones(1,length(x_plot));
% plot(x_flip_and_connect(x_plot),y_flip_and_connect(y_uniform_plot),'-.k','LineWidth',plot_line_width)
% plot(x_flip_and_connect([x_fit1_rauner x_experiment1]),y_flip_and_connect([y_fit1_rauner y_experiment1]),'-b','LineWidth',plot_line_width)
% 为了既占住第一个图例，又不被覆盖，重画一遍
scatter(x_axial_competion(x_experiment1),y_axial_competion(y_experiment1),'o','MarkerEdgeColor','r','LineWidth',plot_line_width)
% 图像格式控制
% title(['ne(z) of CHARLIE \rmat \rm' now_str])
axis([-200,200,0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
xlabel('\itz \rm[mm]')
ylabel('Normalized \itn_{\rme}')
set(gca,'FontSize',font_size)
set(gca, 'LineWidth',gca_line_width)
grid on%显示网格
text(-250,-0.12,'(c)','color','k','FontSize',font_size) %在指定点添加文本

L1=legend([h1,h2],'Experiment','Fit'); %图例。位置为最佳，可以不写，有需要再手调。
set(L1,'FontSize',font_size);
set(L1,'location','south');
set(L1,'Orientation','horizontal');
set(L1,'box','off')
set(L1,'AutoUpdate','off')

axes2 = axes('position',get(gca,'position'),'visible','off');
L2=legend(axes2,[h3,h4], 'Nonuniform case','Uniform case');
set(L2,'FontSize',font_size);
set(L2,'location','south');
set(L2,'box','off')
set(L2,'Orientation','horizontal');
                    
% 径向
x_plot=[45.5, 45.4999, 40, 35, 30.0001, 30, 25, 20, 15.0001, 15, 10, 5, 0];
y_nonuniform_plot=[0, 0.360763542, 0.360763542, 0.360763542, 0.360763542, 0.7961105, 0.7961105, 0.7961105, 0.7961105, 0.9694716, 0.9694716, 0.9694716, 0.9694716];
figure
scatter(x_experiment2,y_experiment2,'o','MarkerEdgeColor','r','LineWidth',plot_line_width)
hold on 
plot(x_fit2_line,y_fit20_line,'-b','LineWidth',plot_line_width)
plot(x_plot,y_nonuniform_plot,'--r','LineWidth',plot_line_width)
line([0,45.5],[0.47,0.47],'Color','k','linestyle','-.','LineWidth',plot_line_width)
% 为了既占住第一个图例，又不被覆盖，重画一遍
scatter(x_experiment2,y_experiment2,'o','MarkerEdgeColor','r','LineWidth',plot_line_width)
% 图像格式控制
% title(['fit of ne(r) of CHARLIE \rmat \rm' now_str])
axis([0,45.5,0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
xlabel('{\it\bf\rho} [mm]')
ylabel('Normalized \itn_{\rme}')
set(gca,'FontSize',font_size)
set(gca, 'LineWidth',gca_line_width)
grid on%显示网格
text(-5,-0.12,'(b)','color','k','FontSize',font_size) %在指定点添加文本

L1=legend('PIC/MCC','Fit','Nonuniform case','Uniform case');
set(L1,'FontSize',font_size);
set(L1,'location','southwest');
% set(L1,'Orientation','horizontal');
set(L1,'box','off')

end

%% aid function
% 使用函数，主要是方便测试和复用。
% 不使用函数，主要是懒得传参

function out_arr=y_radial_competion(in_arr)
% 将行向量in_arr与flip(in_arr)拼接，去除中央的那一个重复值
temp_arr=flip(in_arr);
out_arr=[in_arr temp_arr(2:end)];
end

function out_arr=x_axial_competion(in_arr)
% 将行向量in_arr,偏置-200，然后与flip(in_arr)拼接，去除中央的那一个重复值
temp_arr=flip(in_arr);
out_arr=[in_arr-200 200-temp_arr(2:end)];
end

function out_arr=y_axial_competion(in_arr)
% 将行向量in_arr与flip(in_arr)拼接，去除中央的那一个重复值
temp_arr=flip(in_arr);
out_arr=[in_arr temp_arr(2:end)];
end