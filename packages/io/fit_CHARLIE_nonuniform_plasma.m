% fit spatially nonuniform ne/Te of plasma in CHARLIE

%% Introduction
close all
clear
now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
disp(now_str)

z_fit_point=[0,70:10:200];
z_fit_line=0:1:200;
%% 1. axial norm ne from 2019Rauner-fig7.2
z_ne=[200, 130, 120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0];
ne_z=[0, 0.196224, 0.23341, 0.286613, 0.340961, 0.405034, 0.476545, 0.553776, 0.639588, 0.727689, 0.816362, 0.891304, 0.947941, 0.98913, 1.00229];

% fit1 by Rauner
z_fit1_rauner=180:-20:140;
ne_fit1_rauner=[0.0446224,0.0926773,0.15389];

% fit10 k*Lorentzian([z0,half_gamma],z)
doing_str='fit10: Lorentzian, read half_gamma from figure';
z0=0;
half_gamma0=76.36; %手测半高宽，可见拟合结果误差较大
L_max=Lorentzian([z0,half_gamma0],z0);
k0=max(ne_z)/L_max;
y_fit10_point= k0*Lorentzian([z0,half_gamma0],z_fit_point);
y_fit10_line=k0*Lorentzian([z0,half_gamma0],z_fit_line);
fprintf('%s: %s = %.2f ,%s = %.2f ,%s = %.2f\n',doing_str,'k',k0,'z0',z0,'half_gamma',half_gamma0);

% fit11 k*Lorentzian([z0,half_gamma],z)
doing_str='fit11: Lorentzian, fit k and half_gamma in Matlab';
func_fit1 = @(coff_vec,x_vec) coff_vec(1)*Lorentzian([z0,coff_vec(2)],x_vec);
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
coff_vec0=[half_gamma0,k0];
[coff_vec1,~]= nlinfit(z_ne,ne_z,func_fit1,coff_vec0,opts);
%f:符号函数句柄,如果是以m文件的形式调用的时候,需要加@.
% 拟合参数的标准是(f-y)^2取最小值
%a:最开始预估的值(预拟合的未知参数的估计值)。
k1=coff_vec1(1);
half_gamma1=coff_vec1(2);
y_fit11_point= func_fit1(coff_vec1,z_fit_point);
y_fit11_line=func_fit1(coff_vec1,z_fit_line);
fprintf('%s: %s = %.2e ,%s = %.2e ,%s = %.2e\n',doing_str,'k',k1,'z0',z0,'half_gamma',half_gamma1);

% fit12 k*Lorentzian([z0,half_gamma],z)
doing_str='fit12: Lorentzian, fit k, half_gamma in 1stOpt';
k2=234.98;
half_gamma2=72.88; 
fprintf('%s: %s = %.2e ,%s = %.2e ,%s = %.2e\n',doing_str,'k',k2,'z0',z0,'half_gamma',half_gamma2);
y_fit12_point= k2*Lorentzian([z0,half_gamma2],z_fit_point);
y_fit12_line=k2*Lorentzian([z0,half_gamma2],z_fit_line);

% fit13 k*Lorentzian([z0,half_gamma],z)
doing_str='fit13: Lorentzian, fit half_gamma in 1stOpt';
k3=k0;
half_gamma3=74.73; 
y_fit13_point= k3*Lorentzian([z0,half_gamma3],z_fit_point);
y_fit13_line=k3*Lorentzian([z0,half_gamma3],z_fit_line);
fprintf('%s: %s = %.2e ,%s = %.2e ,%s = %.2e\n',doing_str,'k',k3,'z0',z0,'half_gamma',half_gamma3);

% fit14 Σx^n
doing_str='fit14: polyfit in 1stOpt'; % older version 
y_poly14=@(x)-6.2861e-3+4.0769e-3.*x-6.6033e-5.*x.^2+8.4313e-7*x.^3-2.4400e-9*x.^4;
y_fit14_point= y_poly14(200-z_fit_point);
y_fit14_line=y_poly14(200-z_fit_line);
y_mean14 = integral(y_poly14,200,0)/-200;

% fit15 Σx^n
doing_str='fit15: polyfit in Matlab cftool';
y_poly15=@(x)1.003+0.0002894*x-0.0001769*x.^2+1.397e-06*x.^3-3.226e-09*x.^4;
y_fit15_point= y_poly15(z_fit_point);
y_fit15_line=y_poly15(z_fit_line);
y_mean15 = integral(y_poly15,0,200)/200;

% 对比
figure
scatter(z_ne,ne_z,'o')
hold on
plot(z_fit_line,y_fit10_line,'r-','LineWidth',6)
plot(z_fit_line,y_fit11_line,'b-.','LineWidth',5)
plot(z_fit_line,y_fit12_line,'y-.','LineWidth',4)
plot(z_fit_line,y_fit13_line,'k--')
plot(z_fit_line,y_fit14_line,'g-.')
plot(z_fit_line,y_fit14_line,'c--')
scatter(z_fit1_rauner,ne_fit1_rauner,'+')
% 图像格式控制
title(['fit of ne(z) of CHARLIE \rmat \rm' now_str])
xlabel('z')
ylabel('norm ne')
grid on%显示网格
line([0,200],[y_mean14,y_mean14],'linestyle',':','color','b');
line([0,200],[y_mean15,y_mean15],'linestyle',':','color','r');
legend('experiment','fit10','fit11','fit12','fit13','fit14','fit15','Rauner',...
    num2str(y_mean14),num2str(y_mean15),'Location','best') %图例。位置为最佳，可以不写，有需要再手调。

rms_e14=norm(y_fit14_point-ne_z); %2.4413
rms_e15=norm(y_fit15_point-ne_z); %2.3844

% 综上，使用fit15 多项式拟合 
% 详细对比拟合是没有必要的，直接cftool比较就好
get_ne_z=@(z) 1.003+0.0002894*z-0.0001769*z.^2+1.397e-06*z.^3-3.226e-09*z.^4;

%% 2. axial norm Te from 2019Rauner-fig7.2
z_Te=[130, 120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0];
Te_z=[0.750324414, 0.773111028, 0.780272704, 0.801105955, 0.833982816, 0.866860659, 0.89778514, 0.930337587, 0.954425788, 0.968097363, 0.986001062, 0.988930615, 0.992184582, 1];
% 对比：见cftool
get_Te_z=@(z) 672.1*Lorentzian([0,214.4],z);

%% 3. radial norm ne from 2017Raunera-fig4
r_ne=[45.5, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0];
ne_r=[0, 0.035130017, 0.180627079, 0.393114122, 0.5617439, 0.699048337, 0.798758892, 0.865696707, 0.899862222, 0.939810461, 1];
% 对比：见cftool
get_ne_r=@(r) 1.002-0.01921*r+0.001517*r.^2-6.698e-05*r.^3+7.112e-07*r.^4;

%% 4. radial norm Te from 2017Raunera-fig4
r_Te=[45, 40, 35, 30, 25, 20, 15, 10, 5, 0];
Te_r=[0.973151812, 1, 0.971381704, 0.94632216, 0.91439163, 0.879948658, 0.854767708, 0.837248269, 0.83310566, 0.829764696];
% 对比：见cftool
get_Te_r=@(r) 0.8296+0.0007247*r-8.095e-05*r.^2+1.27e-05*r.^3-2.147e-07*r.^4;