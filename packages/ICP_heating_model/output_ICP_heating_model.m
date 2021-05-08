function output_ICP_heating_model( flag,plasma )
% 对ICP heating model的输出与可视化
if isfield(flag,'stoc_model') && ~isempty(flag.stoc_model)
    flag_use_stoc_model=true;
else
    flag_use_stoc_model=false;
end

if 1==plasma.size
    %% single point
    fprintf('%s = %.2e , ','ν_m',plasma.nu_m);
    if flag_use_stoc_model
        fprintf('%s = %.2e , ','ν_stoc',plasma.nu_st);
    end
    fprintf('%s = %.2e \n','ν_eff',plasma.nu_eff);
    
    disp('碰撞频率')
    fprintf('%s = %.2e , ','ν_enp',plasma.nu_enp);
    fprintf('%s = %.2e , ','ν_eniz',plasma.nu_eniz);
    fprintf('%s = %.2e \n','ν_eip',plasma.nu_eip);
    %         fprintf('%s = %.2e , ','λenp',lambda_coll_enp);
    %         fprintf('%s = %.2e , ','λeniz',lambda_coll_eniz);
    %         fprintf('%s = %.2e \n','λeip',lambda_coll_eip);
    
    % 主导加热机制
    if plasma.nu_m2nu_eff<0.1
        disp('ν_m/ν_eff < 0.1, 加热机制以随机加热为主.')
    elseif plasma.nu_m2nu_eff<=0.9
        disp('需同时考虑欧姆加热和随机加热.')
    elseif flag_use_stoc_model
        disp('ν_stoc/ν_eff > 0.9, 加热机制以欧姆加热为主.')
    end
    if flag_use_stoc_model
        fprintf('%s = %.2e , ','趋肤层渡越时间τ',plasma.delta_st/plasma.ve);
    end
    
    % 带电粒子响应RF电磁场
    if flag_use_stoc_model && plasma.wpi2wRF>=1
        fprintf('%s = %.2e \n','ω_pi/ω_RF',plasma.wpi2wRF);
        disp('[WARN] 1995Vahedia的stoc model要求离子不响应电场, 但当前等离子体参数集并不满足')
    end
    disp('')
else
    %% multi point
%     % 对比不同加热机制
%     if ~isempty(flag.stoc_model)
%         name_Y='\nu_c';
%         Y1=plasma.nu_m;
%         name_Y1='\nu_{\rmm}\rm[rad/s]'; %可能是Hz
%         Y2=plasma.nu_st;
%         name_Y2='\nu_{\rmst}\rm[rad/s]';
%         
%         
%         handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
%             X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%         hold on
%         line([X1(1),X1(end)],[w_RF,w_RF],'linestyle','-.','linewidth',1.5,'color','k');
%         xlabel([X_var '\rm[' unit_X_var ']']);
%         ylabel(name_Y);
%         L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
%             [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
%         title([name_Y ' \rmat \rm' now_str]);
%         %         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
%         
%         
%         L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
%             [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
%             '\it\omega_{\rmRF}');
%         set(L1,'location','northwest');
%         set(L1,'box','off')
%     end
%     
%     Y1=veff_per_w;
%     name_Y1='\it\nu_{\rmeff}/\it\omega_{\rmRF}';
%     handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%     title([name_Y1 ' \rmat \rm' now_str ', if>>1, \sigma_{\rmDC}可用']);
%     line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
%     line([X1(1),X1(end)],[10,10],'linestyle','-.','linewidth',1.5,'color','k');
%     legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
%     pause
%     
%     % 欧姆加热不同种类碰撞频率-for 聚变负源典型参数
%     figure
%     loglog(X1,venp(:,1),'-r','LineWidth',1.5);
%     hold on
%     loglog(X1,veniz(:,1),'--r','LineWidth',1.5);
%     hold on
%     loglog(X1,veip(:,1),'-.r','LineWidth',1.5);
%     hold on
%     loglog(X1,veip(:,no_mid_X2),'-.b','LineWidth',1.5);
%     hold on
%     loglog(X1,veip(:,end),'-.g','LineWidth',1.5);
%     xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
%     ylabel('碰撞频率');
%     L1=legend('venp','veniz',...
%         ['veip at ' legend_X2{1}],['veip at ' legend_X2{no_mid_X2}],['veip at ' legend_X2{end}]);
%     set(L1,'FontSize',10);
%     title(['不同种类碰撞频率 \rmat \rm' now_str]);
%     axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
%     grid on%显示网格
    
end

end