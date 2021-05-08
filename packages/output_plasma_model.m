function output_plasma_model(flag,plasma)
% output information of plasma model
% mainly for single point now

fprintf('[INFO] Results from plasma model.\n');
if 1==plasma.size
    % single point
    disp('等离子体参数')
    fprintf('%s = %.1e m^-3 , ','ne',plasma.ne);
    fprintf('%s = %.1f eV\n','Te',plasma.Te);
    fprintf('%s = %.2e Hz, ','f',plasma.f);
    fprintf('%s = %.2e Pa, ','p',plasma.p);
    fprintf('%s = %.2e K, ','Tg',plasma.Tg)
    fprintf('%s = %.2e m^-3\n','ng',plasma.ng);
    disp('特征频率')
    fprintf('%s = %.2e rad/s, ','ω_RF',plasma.w_RF);
    fprintf('%s = %.2e rad/s, ','ω_pe',plasma.wpe);
    fprintf('%s = %.2e rad/s\n','ω_pi',plasma.wpi);
    disp('特征时间')
    fprintf('%s = %.2e s\n','RF周期T',2*pi/plasma.w_RF);
else
    % multi point
    disp('等离子体参数')
    fprintf('%s = %.1e ~ %.1e m^-3 , ','ne',min(plasma.ne(:)),max(plasma.ne(:)));
    fprintf('%s = %.1f ~ %.1f eV\n','Te',min(plasma.Te(:)),max(plasma.Te(:)));
    fprintf('%s = %.2e ~ %.2e Hz, ','f',min(plasma.f(:)),max(plasma.f(:)));
    fprintf('%s = %.1f ~ %.1f Pa, ','p',min(plasma.p(:)),max(plasma.p(:)));
    fprintf('%s = %.1f ~ %.1f K, ','Tg',min(plasma.Tg(:)),max(plasma.Tg(:)));
    fprintf('%s = %.2e ~ %.2e m^-3\n','ng',min(plasma.ng(:)),max(plasma.ng(:)));
    disp('特征频率')
    fprintf('%s = %.2e ~ %.2e rad/s, ','ω_RF',min(plasma.w_RF(:)),max(plasma.w_RF(:)));
    fprintf('%s = %.2e ~ %.2e rad/s, ','ω_pe',min(plasma.wpe(:)),max(plasma.wpe(:)));
    fprintf('%s = %.2e ~ %.2e rad/s\n','ω_pi',min(plasma.wpi(:)),max(plasma.wpi(:)));
    disp('特征时间')
    fprintf('%s = %.2e ~ %.2e s\n','RF周期T',min(2*pi./plasma.w_RF(:)),max(2*pi./plasma.w_RF(:)));
end
% child model
output_ICP_heating_model( flag, plasma );
output_equivalent_EM_medium_model( flag, plasma )
end

%% 留存
% 部分图需要手动调整图形大小以避免legend覆盖曲线，因此不自动保存
%             save_path='d:\School\DoctorProgram\eP-项目笔记文件夹\eP-190821-01激励器FEM模型\200221期刊论文\';
%             saveas(gcf,[save_path name_Y '.svg'],'svg')

%字符数组不能够存储不同长度的字符串，因此用元胞
%     X_var={'ne';'Te';'p'};
%     name_X_var={'\itn\rm_e';'\itT\rm_e';'\itp'};
%     unit_X_var={'m^{-3}';'eV';'Pa'};
%
%     output.idx_X1=1; %绘图时第一自变量，即横轴
%     output.X1=plasma.ne;
%     output.idx_X2='T_e'; %绘图时第二自变量，即legend
%     output.no_mid_X2=3;
