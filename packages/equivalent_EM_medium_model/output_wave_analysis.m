function output_wave_analysis( plasma )
% 对wave analysis的输出与可视化
if 1==plasma.size
%% 单数据点：输出文本
    disp('特征长度')
    fprintf('%s = %.2e \n','debye length',get_debye_length(plasma.ne,plasma.Te));
    fprintf('%s = %.2e , ','effective_skin_depth',plasma.skin_depth);
    fprintf('%s = %.2e \n','λ_wave',plasma.wavelength);
    fprintf('%s = %.2e \n','r_plasma',plasma.r);
    disp('')
else
%% 多数据点：输出图像
%             Y1=skin_depth_eff;
%             name_Y1='\it\delta_{\rmeff}\rm[m]';
%             handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%             title([name_Y1 ' \rmat \rm' now_str ', if \it\delta_{\rmeff}>\itr_{\rmchamber},电磁波穿透等离子体']);
%             line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle','-.','linewidth',1.5,'color','k');
%             legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end},'\itr_{\rmchamber}');
%             
%             Y1=wavelength_wave;
%             name_Y1='\it\lambda_{\rmwave}\rm[m]';
%             handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%             title([name_Y1 ' \rmat \rm' now_str ', if \it\lambda_{\rmwave}<\itr_{\rmchamber},各处不同相位']);
%             line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle','-.','linewidth',1.5,'color','k');
%             legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end},'\itr_{\rmchamber}');
%             pause
end