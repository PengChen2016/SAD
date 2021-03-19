function output_plasma_model(flag, plasma)
% 对plasma model的输出与可视化

if 1==plasma.size
    fprintf('%s = %.2e , ','[INFO] Results from used ICP heating model: ν_m= , ν_st= \n',plasma.nu_m,plasma.nu_st);
    %用于单点计算
    disp('等离子体参数')
    fprintf('%s = %.2.e , ','ne',plasma.ne);
    fprintf('%s = %.2e \n','Te',plasma.Te);
    fprintf('%s = %.2e ,','f',plasma.f);
    fprintf('%s = %.2e , ','p',plasma.p);
    fprintf('%s = %.2e , ','Tg',plasma.Tg)
    fprintf('%s = %.2e \n','ng',plasma.ng);
    disp('特征频率')
    fprintf('%s = %.2e , ','ν_m',plasma.nu_m);
    fprintf('%s = %.2e , ','vst',plasma.nu_st);
    fprintf('%s = %.2e \n','veff',plasma.nu_eff);
    fprintf('%s = %.2e , ','ω_RF',plasma.w_RF);
    fprintf('%s = %.2e , ','ωpe',plasma.wpe);
    fprintf('%s = %.2e \n','ωpi',plasma.wpi);
    disp('模型适用条件') %参考利伯曼教材
    disp('if wp/w_RF > 1, 带电粒子响应RF电场')
    fprintf('%s = %.2e , ','wpe/w_RF',plasma.wpe2wRF);
    fprintf('%s = %.2e \n','wpi/w_RF',plasma.wpi2wRF);
    if wpi_per_w>1
        disp('1995Vahedia的ICP heating model要求wpi<w_RF,但当前等离子体参数集并不满足')
    end
    disp('if vm/veff << 1, 加热机制以随机加热为主.')
    fprintf('%s = %.2e \n','vm/veff',plasma.nu_m2nu_eff);
    disp('')
    disp('复介电常数')
    disp(['eps_c_r = ' num2str(plasma.eps_c_r,'%.2e')])
    fprintf('%s = %.2e , ','eps_r',plasma.eps_r);
    fprintf('%s = %.2e \n','tandelta',plasma.tan_delta);
    
    disp('复电导率')
    disp(['sigmap = ' num2str(plasma.sigma_p,'%.2e')])
    disp('if v/w_RF >> 1, v+jω≈v,可以使用sigma_dc表达式')
    fprintf('%s = %.2e \n','veff/w_RF',plasma.nu_eff/plasma.w_RF);
    fprintf('%s = %.2e \n','sigma_dc',plasma.sigma_dc);
    disp('if wpe^2/(w_RF*sqrt(w_RF^2+v^2)) >> 1, then sigmap without jw*constants.eps0 can be used.')
    fprintf('%s = %.2e \n','for v=veff,',plasma.ratio_displacement_current);
    disp('')
    disp('特征长度')
    fprintf('%s = %.2e \n','debye length',get_debye_length(plasma.ne,plasma.Te));
    fprintf('%s = %.2e , ','effective_skin_depth',plasma.skin_depth);
    fprintf('%s = %.2e \n','wavelength',plasma.wavelength);
    fprintf('%s = %.2e \n','r_plasma',plasma.r);
    disp('特征时间')
    fprintf('%s = %.2e , ','趋肤层渡越时间τ',plasma.delta_st/plasma.ve);
    fprintf('%s = %.2e \n','RF周期T',2*pi/plasma.w_RF);
    
    disp('碰撞分析')
    fprintf('%s = %.2e , ','venp',plasma.nu_enp);
    fprintf('%s = %.2e , ','veniz',plasma.nu_eniz);
    fprintf('%s = %.2e \n','veip',plasma.nu_eip);
    %         fprintf('%s = %.2e , ','λenp',lambda_coll_enp);
    %         fprintf('%s = %.2e , ','λeniz',lambda_coll_eniz);
    %         fprintf('%s = %.2e \n','λeip',lambda_coll_eip);
    
else
    
    if flag_output_for_paper
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output for paper
        
        
        return
        
        % 部分图需要手动调整图形大小以避免legend覆盖曲线，因此不自动保存
        %             save_path='d:\School\DoctorProgram\eP-项目笔记文件夹\eP-190821-01激励器FEM模型\200221期刊论文\';
        %             saveas(gcf,[save_path name_Y '.svg'],'svg')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output for paper
    else
        if flag.data_map
            %         % benchmark with 2019Raunera fig7a
            %         vm_ref=[5.12E+06,7.71E+06,1.52E+07,4.45E+07,7.27E+07,1.37E+08];
            %         vst_ref=[7.77E+06,9.19E+06,9.58E+06,9.58E+06,9.17E+06,8.27E+06];
            %                 figure
            %         loglog(X1,plasma.nu_m,'-r','LineWidth',1.5);
            %         hold on
            %         loglog(X1,vst,'-b','LineWidth',1.5);
            %         hold on
            %         loglog(X1,vm_ref,'--dr','LineWidth',1.5);
            %         hold on
            %         loglog(X1,vst_ref,'--db','LineWidth',1.5);
            %         xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            %         ylabel('\it\nu\rm[Hz]');
            %         L1=legend('\it\nu\rm_{m}','\it\nu\rm_{st}','\it\nu\rm_{m}-Rauner','\it\nu\rm_{st}-Rauner');
            %         set(L1,'FontSize',10);
            %         set(L1,'location','southeast');
            %         %         set(L1,'box','off')
            %         title(['碰撞频率 \rmat \rm' now_str]);
            %         %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            %         grid on%显示网格
            %         pause
            
            % 等离子体等效电磁参数
            figure
            yyaxis left
            semilogx(X1,sigma_medium,'-r','LineWidth',1.5);
            ylabel('\it\sigma_{\rmeff}');
            yyaxis right
            semilogx(X1,-eps_medium/constants.eps0,'-b','LineWidth',1.5);
            ylabel('\it-\epsilon_{\rmr-eff}');
            xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            title(['等离子体等效电磁参数 \rmat \rm' now_str]);
            L1=legend('\it\sigma_{\rmeff}','\it-\epsilon_{\rmr-eff}');
            set(L1,'FontSize',10);
            %                     set(L1,'location','southeast');
            %         set(L1,'box','off')
            %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            grid on%显示网格
            %             fig_background_transparent(gcf,gca)  %临时使用
            
            figure
            semilogx(X1,sigma_medium./(-w_RF*eps_medium),'-r','LineWidth',1.5);
            ylabel('\it\sigma_{\rmeff}/\it\omega|\epsilon_{\rmeff}|');
            xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            title(['\it\sigma_{\rmeff}/\it\omega|\epsilon_{\rmeff}|' ' \rmat \rm' now_str]);
            line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
            line([X1(1),X1(end)],[10,10],'linestyle','-.','linewidth',1.5,'color','k');
            %                     L1=legend('\it\sigma_{\rmeff}','\it-\epsilon_{\rmr-eff}');
            %                     set(L1,'FontSize',10);
            %                     set(L1,'location','southeast');
            %         set(L1,'box','off')
            %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            grid on%显示网格
            pause
        else
            disp('模型参数')
            fprintf('%s = %.2e ,','f',f);
            fprintf('%s = %.2e , ',X_var{idx_X3},X3);
            fprintf('%s = %.2e , ','Tg',Tg)
            fprintf('%s = %.2e~%.2e , ',X_var{idx_X1},X1(1),X1(end))
            fprintf('%s = %.2e~%.2e\n',X_var{idx_X2},X2(1),X2(end))
            
            % TODO：除Te外，全部为SI单位，待添加进图
            %                     name_Y='等离子体等效电磁参数';
            %                     Y1=sigma_medium;
            %                     name_Y1='\it\sigma_{\rmeff}';
            %                     Y2=-w_RF*eps_medium;
            %                     name_Y2='\it-\omega\epsilon_{\rmeff}';
            %                     handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
            %                         X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            %                     % 200529 试过logy的图，不怎么好看，只是看上去近了一点
            %                     L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
            %                         [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
            %                     set(L1,'location','northwest');
            
            name_Y='等离子体等效电磁参数';
            Y1=sigma_medium;
            name_Y1='\it\sigma_{\rmeff}';
            Y2=-eps_medium/constants.eps0;
            name_Y2='\it-\epsilon_{\rmr-eff}';
            handle_fig=plot_parametric_2Ylogaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
            set(L1,'location','northwest');
            %         ne_center=5e18;
            %         ne_skin=ne_center*0.55;
            %         ne_axi_ave=ne_center*0.42;
            %         line([,X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
            
            Y1=sigma_medium./(-w_RF*eps_medium);
            name_Y1='\it\sigma_{\rmeff}/\it\omega|\epsilon_{\rmeff}|';
            handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            title([name_Y1 ' \rmat \rm' now_str]);
            line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
            line([X1(1),X1(end)],[10,10],'linestyle','-.','linewidth',1.5,'color','k');
            L1=legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
            set(L1,'location','northwest');
            
            pause
            
            %                 name_Y='有损介质参数';
            %                 Y1=-epsp_r_real;
            %                 name_Y1='\it-\epsilon_{\rmr}';
            %                 Y2=-tandelta;
            %                 name_Y2='\rm-tan\it\delta';
            %                 handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
            %                     X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            %                 L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
            %             [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
            %                 set(L1,'location','northwest');
            %                 pause
            
            %         name_Y='复电导率\it\sigma_{\rmp}';
            %         Y1=sigmap_real;
            %         name_Y1='\rmRe(\it\sigma_{\rmp}\rm)[S/m]';
            %         Y2=-sigmap_imag;
            %         name_Y2='-\rmIm(\it\sigma_{\rmp}\rm) = - \it\omega\epsilon[S/m]';
            %         handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
            %             X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            %         name_Y2='-\rmIm(\it\sigma_{\rmp}\rm)[S/m]';
            %         handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
            %             X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            %         pause
            
            Y1=skin_depth_eff;
            name_Y1='\it\delta_{\rmeff}\rm[m]';
            handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            title([name_Y1 ' \rmat \rm' now_str ', if \it\delta_{\rmeff}>\itr_{\rmchamber},电磁波穿透等离子体']);
            line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle','-.','linewidth',1.5,'color','k');
            legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end},'\itr_{\rmchamber}');
            
            Y1=wavelength_wave;
            name_Y1='\it\lambda_{\rmwave}\rm[m]';
            handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            title([name_Y1 ' \rmat \rm' now_str ', if \it\lambda_{\rmwave}<\itr_{\rmchamber},各处不同相位']);
            line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle','-.','linewidth',1.5,'color','k');
            legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end},'\itr_{\rmchamber}');
            pause
            
            name_Y='带电粒子对电磁场的响应';
            Y1=wpi_per_w;
            name_Y1='\it\omega_{\rmpi}/\it\omega_{\rmRF}';
            Y2=wpe_per_w;
            name_Y2='\it\omega_{\rmpe}/\it\omega_{\rmRF}';
            handle_fig=plot_parametric_2Ylog(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
            L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                '\it\omega_{\rmpi}/\it\omega_{\rmRF}\rm=1');
            set(L1,'location','northwest');
            pause
            
            name_Y='碰撞频率';
            Y1=plasma.nu_m;
            name_Y1='\nu_{\rmm}\rm[Hz]'; %可能是rad/s
            Y2=vst;
            name_Y2='\nu_{\rmst}\rm[Hz]';
            handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            line([X1(1),X1(end)],[w_RF,w_RF],'linestyle','-.','linewidth',1.5,'color','k');
            L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                '\it\omega_{\rmRF}');
            set(L1,'location','northwest');
            set(L1,'box','off')
            pause
            
            Y1=veff_per_w;
            name_Y1='\it\nu_{\rmeff}/\it\omega_{\rmRF}';
            handle_fig=plot_parametric_1Y(Y1,name_Y1,X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            title([name_Y1 ' \rmat \rm' now_str ', if>>1, \sigma_{\rmDC}可用']);
            line([X1(1),X1(end)],[1,1],'linestyle','-.','linewidth',1.5,'color','k');
            line([X1(1),X1(end)],[10,10],'linestyle','-.','linewidth',1.5,'color','k');
            legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
            pause
            
            % 欧姆加热不同种类碰撞频率-for 聚变负源典型参数
            figure
            loglog(X1,venp(:,1),'-r','LineWidth',1.5);
            hold on
            loglog(X1,veniz(:,1),'--r','LineWidth',1.5);
            hold on
            loglog(X1,veip(:,1),'-.r','LineWidth',1.5);
            hold on
            loglog(X1,veip(:,no_mid_X2),'-.b','LineWidth',1.5);
            hold on
            loglog(X1,veip(:,end),'-.g','LineWidth',1.5);
            xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            ylabel('碰撞频率');
            L1=legend('venp','veniz',...
                ['veip at ' legend_X2{1}],['veip at ' legend_X2{no_mid_X2}],['veip at ' legend_X2{end}]);
            set(L1,'FontSize',10);
            title(['不同种类碰撞频率 \rmat \rm' now_str]);
            axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            grid on%显示网格
            
            %                         % 欧姆加热不同种类碰撞频率-for 2014Cazzador
            %                         figure
            %                         loglog(X1,venp(:,1),'-r','LineWidth',1.5);
            %                         hold on
            %                         loglog(X1,veniz(:,1),'--r','LineWidth',1.5);
            %                         hold on
            %                         loglog(X1,veip(:,1),'-.r','LineWidth',1.5);
            %                         xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            %                         ylabel('碰撞频率');
            %                         L1=legend('venp','veniz',['veip at ' legend_X2{1}]);
            %                         set(L1,'FontSize',10);
            %                         title(['不同种类碰撞频率 \rmat \rm' now_str]);
            %                         axis([X1(1),X1(end),1e6,1e8]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            %                         grid on%显示网格
            
            % 测试stoc模型
            % stoc模型参数
            name_Y='stoc模型参数';
            Y1=alpha_st;
            name_Y1='\rm\alpha_{\rmst}';
            Y2=delta_st;
            name_Y2='\rm\delta_{\rmst}';
            handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            yyaxis left;
            line([X1(1),X1(end)],[0.03,0.03],'linestyle','-.','linewidth',1.5,'color','k');
            L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                '\it\alpha\rm=0.03',...
                [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
            %         set(L1,'location','northwest');
            
            % vst简化表达式与Cazzador fit表达式对比
            switch flag.stoc_model
                case {'Vahedi-simplify','Cazzador-simplify'}
                    name_Y='\it\nu_{\rmst}';
                    Y1=vst_fit;
                    name_Y1='Cazzador fit';
                    Y2=vst;
                    name_Y2=flag.stoc_model;
                    handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                        X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                    line([X1(1),X1(end)],[w_RF,w_RF],'linestyle','-.','linewidth',1.5,'color','k');
                    L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                        [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                        '\it\omega_{\rmRF}');
                    set(L1,'location','northwest');
                    
                    name_Y='\it\delta_{\rmst}';
                    Y1=delta_st_fit;
                    name_Y1='Cazzador-fit';
                    Y2=delta_st;
                    name_Y2=flag.stoc_model;
                    handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                        X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
                    %                     line([X1(1),X1(end)],[w_RF,w_RF],'linestyle','-.','linewidth',1.5,'color','k');
                    %                     L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                    %                         [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                    %                         '\it\omega_{\rmRF}');
                    %                     set(L1,'location','northwest');
            end
            pause
        end
    end
end
end