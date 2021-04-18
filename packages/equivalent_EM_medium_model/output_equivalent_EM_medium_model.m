function output_equivalent_EM_medium_model( flag, plasma )
% 对equivalent EM_medium model的输出与可视化
if 1==plasma.size
    %% 单数据点：输出文本
    
    disp('复介电常数')
    disp(['eps_c_r = ' num2str(plasma.eps_c_r,'%.2e')])
    fprintf('%s = %.2e , ','eps_r',plasma.eps_r);
    fprintf('%s = %.2e \n','tandelta',plasma.tan_delta);
    
    % 复电导率
    fprintf('%s = %.2e \n','位移电流占比',plasma.ratio_displacement_current);
    if plasma.ratio_displacement_current<0.01
        disp('ω_pe^2/(ω_RF*sqrt(ω_RF^2+ν^2)) >> 1, 可忽略等离子体位移电流，ε_p近似等效为σ_p')
        disp(['σ_p = ' num2str(plasma.sigma_p,'%.2e')])
    end
    
    % 直流电导率
    fprintf('%s = %.2e \n','veff/w_RF',plasma.nu_eff/plasma.w_RF);
    if plasma.nu_eff/plasma.w_RF>10
        disp('ν/ω_RF >> 1, 可近似ν+jω≈ν, σ_p近似为σ_dc')
        fprintf('%s = %.2e \n','σ_dc',plasma.sigma_dc);
    else
        fprintf('%s = %.2e \n','ν/ω_RF, 不应使用σ_dc\n',plasma.nu_eff/plasma.w_RF);
        if isfield(flag,'medium_approximation') && strcmp(flag.medium_approximation,'sigma_dc')
            warning('The used medium_approximation sigma_dc is not suitable.')
        end
    end
    
    % 带电粒子响应RF电磁场
    fprintf('%s = %.2e , ','ω_pe/ω_RF',plasma.wpe2wRF);
    fprintf('%s = %.2e \n','ω_pi/ω_RF',plasma.wpi2wRF);
    if plasma.wpe2wRF>1
        disp('ω_pe/ω_RF > 1, 电子运动响应RF电场')
    end
    if plasma.wpi2wRF<1
        disp('ω_pi/ω_RF < 1, 离子运动跟不上RF电场变化，近似为静止')
    elseif plasma.wpe/plasma.wpi>10
        fprintf('%s = %.2e >>1,  可忽略离子动力学\n','ω_pe/ω_pi',plasma.wpe/plasma.wpi)
    else
        fprintf('[WARN] %s = %.2e <10,  不可忽略离子动力学，现有电磁描述不适用\n',...
            'ω_pe/ω_pi',plasma.wpe/plasma.wpi)
    end
    disp('')
else
    %% 多数据点：输出图像

    
%     % TODO：除Te外，全部为SI单位，待添加进图
%     %                     name_Y='等离子体等效电磁参数';
%     %                     Y1=sigma_medium;
%     %                     name_Y1='\it\sigma_{\rmeff}';
%     %                     Y2=-w_RF*eps_medium;
%     %                     name_Y2='\it-\omega\epsilon_{\rmeff}';
%     %                     handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
%     %                         X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%     %                     % 200529 试过logy的图，不怎么好看，只是看上去近了一点
%     %                     L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
%     %                         [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
%     %                     set(L1,'location','northwest');
%     
%     name_Y='等离子体等效电磁参数';
%     Y1=sigma_medium;
%     name_Y1='\it\sigma_{\rmeff}';
%     Y2=-eps_medium/constants.eps0;
%     name_Y2='\it-\epsilon_{\rmr-eff}';
%     handle_fig=plot_parametric_2Ylogaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
%         X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%     L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
%         [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
%     set(L1,'location','northwest');
%     %         ne_center=5e18;
%     %         ne_skin=ne_center*0.55;
%     %         ne_axi_ave=ne_center*0.42;
%     %         line([,X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
%     
%     Y1=sigma_medium./(-w_RF*eps_medium);
%     name_Y1='\it\sigma_{\rmeff}/\it\omega|\epsilon_{\rmeff}|';
%     handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%     title([name_Y1 ' \rmat \rm' now_str]);
%     line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
%     line([X1(1),X1(end)],[10,10],'linestyle','-.','linewidth',1.5,'color','k');
%     L1=legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
%     set(L1,'location','northwest');
%     
%     pause
%     
%     %                 name_Y='有损介质参数';
%     %                 Y1=-epsp_r_real;
%     %                 name_Y1='\it-\epsilon_{\rmr}';
%     %                 Y2=-tandelta;
%     %                 name_Y2='\rm-tan\it\delta';
%     %                 handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
%     %                     X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%     %                 L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
%     %             [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
%     %                 set(L1,'location','northwest');
%     %                 pause
%     
%     %         name_Y='复电导率\it\sigma_{\rmp}';
%     %         Y1=sigmap_real;
%     %         name_Y1='\rmRe(\it\sigma_{\rmp}\rm)[S/m]';
%     %         Y2=-sigmap_imag;
%     %         name_Y2='-\rmIm(\it\sigma_{\rmp}\rm) = - \it\omega\epsilon[S/m]';
%     %         handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
%     %             X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%     %         name_Y2='-\rmIm(\it\sigma_{\rmp}\rm)[S/m]';
%     %         handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
%     %             X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%     %         pause
%     
%     name_Y='带电粒子对电磁场的响应';
%     Y1=wpi_per_w;
%     name_Y1='\it\omega_{\rmpi}/\it\omega_{\rmRF}';
%     Y2=wpe_per_w;
%     name_Y2='\it\omega_{\rmpe}/\it\omega_{\rmRF}';
%     handle_fig=plot_parametric_2Ylog(name_Y,Y1,name_Y1,Y2,name_Y2,...
%         X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%     line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
%     L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
%         [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
%         '\it\omega_{\rmpi}/\it\omega_{\rmRF}\rm=1');
%     set(L1,'location','northwest');
%     pause
%     
%     
%     if flag.data_map
%         % 等离子体等效电磁参数
%         figure
%         yyaxis left
%         semilogx(X1,sigma_medium,'-r','LineWidth',1.5);
%         ylabel('\it\sigma_{\rmeff}');
%         yyaxis right
%         semilogx(X1,-eps_medium/constants.eps0,'-b','LineWidth',1.5);
%         ylabel('\it-\epsilon_{\rmr-eff}');
%         xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
%         title(['等离子体等效电磁参数 \rmat \rm' now_str]);
%         L1=legend('\it\sigma_{\rmeff}','\it-\epsilon_{\rmr-eff}');
%         set(L1,'FontSize',10);
%         %                     set(L1,'location','southeast');
%         %         set(L1,'box','off')
%         %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
%         grid on%显示网格
%         %             fig_background_transparent(gcf,gca)  %临时使用
%         
%         figure
%         semilogx(X1,sigma_medium./(-w_RF*eps_medium),'-r','LineWidth',1.5);
%         ylabel('\it\sigma_{\rmeff}/\it\omega|\epsilon_{\rmeff}|');
%         xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
%         title(['\it\sigma_{\rmeff}/\it\omega|\epsilon_{\rmeff}|' ' \rmat \rm' now_str]);
%         line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
%         line([X1(1),X1(end)],[10,10],'linestyle','-.','linewidth',1.5,'color','k');
%         %                     L1=legend('\it\sigma_{\rmeff}','\it-\epsilon_{\rmr-eff}');
%         %                     set(L1,'FontSize',10);
%         %                     set(L1,'location','southeast');
%         %         set(L1,'box','off')
%         %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
%         grid on%显示网格
%         pause
%     end
end

output_wave_analysis( plasma )

end