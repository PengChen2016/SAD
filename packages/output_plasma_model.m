function output_plasma_model(flag, plasma)
% 对plasma model的输出与可视化

if 1==plasma.size
    fprintf('%s = %.2e , ','[INFO] Results from used ICP heating model: ν_m= , ν_st= \n',plasma.nu_m,plasma.nu_st);
end

    if ~flag.sweep
        %用于单点计算
        disp('等离子体参数')
        fprintf('%s = %constants.e , ','ne',ne);
        fprintf('%s = %.2e \n','Te',Te);
        fprintf('%s = %.2e ,','f',f);
        fprintf('%s = %.2e , ','p',p);
        fprintf('%s = %.2e , ','Tg',Tg)
        fprintf('%s = %constants.e \n','ng',ng);
        disp('特征频率')
        fprintf('%s = %.2e , ','vm',vm);
        fprintf('%s = %.2e , ','vst',vst);
        fprintf('%s = %.2e \n','veff',veff);
        fprintf('%s = %.2e , ','ω_RF',w_RF);
        fprintf('%s = %.2e , ','ωpe',wpe);
        fprintf('%s = %.2e \n','ωpi',wpi);
        disp('模型适用条件') %参考利伯曼教材
        disp('if wp/w_RF > 1, 带电粒子响应RF电场')
        fprintf('%s = %.2e , ','wpe/w_RF',wpe_per_w);
        fprintf('%s = %.2e \n','wpi/w_RF',wpi_per_w);
        if ~flag.sweep && wpi_per_w>1
            disp('1995Vahedia的ICP heating model要求wpi<w_RF,但当前等离子体参数集并不满足')
        end
        disp('if vm_per_veff << 1, 加热机制以随机加热为主.')
        fprintf('%s = %.2e \n','vm_per_veff',vm_per_veff);
        disp('')
        disp('复介电常数')
        disp(['epsp_r = ' num2str(epsp_r,'%.2e')])
        fprintf('%s = %.2e , ','eps_r',epsp_r_real);
        fprintf('%s = %.2e \n','tandelta',tandelta);
        
        disp('复电导率')
        disp(['sigmap = ' num2str(sigmap,'%.2e')])
        disp('if v/w_RF >> 1, v+jω≈v,可以使用sigma_dc表达式')
        fprintf('%s = %.2e \n','veff/w_RF',veff_per_w);
        fprintf('%s = %.2e \n','sigma_dc',sigmaeff_dc);
        disp('if wpe^2/(w_RF*sqrt(w_RF^2+v^2)) >> 1, then sigmap without jw*constants.eps0 can be used.')
        fprintf('%s = %.2e \n','for v=veff,',veff_displacement_ratio);
        disp('')
        disp('特征长度')
        fprintf('%s = %.2e \n','debye length',sqrt(constants.eps0*Te*constants.e/(ne*constants.e*constants.e)));
        fprintf('%s = %.2e , ','effective_skin_depth',skin_depth_eff);
        fprintf('%s = %.2e \n','wavelength_from_epsD',wavelength_wave);
        fprintf('%s = %.2e \n','r_plasma',r_plasma);
        % TODO：待输出δstoc
        disp('特征时间')
        fprintf('%s = %.2e , ','趋肤层渡越时间τ',skin_depth_eff/Ve);
        fprintf('%s = %.2e \n','RF周期T',2*pi/w_RF);
        
        disp('碰撞分析')
        fprintf('%s = %.2e , ','venp',venp);
        fprintf('%s = %.2e , ','veniz',veniz);
        fprintf('%s = %.2e \n','veip',veip);
        fprintf('%s = %.2e , ','λenp',lambda_coll_enp);
        fprintf('%s = %.2e , ','λeniz',lambda_coll_eniz);
        fprintf('%s = %.2e \n','λeip',lambda_coll_eip);
    else
        %         % 输出到未打开的excel表
        %         filename = 'testdata.xlsx';
        % sheet = 2;
        % xlRange = 'E1';
        % xlswrite(filename,A,sheet,xlRange)
        
        if flag_output_for_paper
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% output for paper
            plot_line_width=3;
            gca_line_width=1;
            marker_size=8;
            font_size=15;
            marker_indices=1:5:length(X1);
            left_color = [0 0 0];
            right_color =left_color;
            %             right_color = [0.8500    0.3250    0.0980];
            group_color_list={'r','b','g','constants.c','m','k','y'};
            group_style_list={'-','--','-.',':'};
            
            flag_multi_Te=true;
%             flag_multi_Te=false;
            
            
            % 等离子体等效电磁参数
            name_Y='等离子体等效电磁参数';
            Y1=sigma_medium;
            name_Y1='\it{\bf\sigma}\rm';
            Y2=-eps_medium/constants.eps0;
            name_Y2='-\it{\bf\epsilon}\rm_r';
            handle_fig=figure;
            set(handle_fig,'defaultAxesColorOrder',[left_color; right_color]);
            % marker会遮掩线型，而且也不好看
            flag_logy=true;
            %             flag_logy=false;
            if flag_logy
                % 左右轴取值：为了便于联系两个曲线族与相应坐标轴
                no_group=1;
                yyaxis right;
                loglog(X1,Y1(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                ylabel([name_Y1 ' [S/m]']);
                yyaxis left;
                loglog(X1,Y2(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                ylabel(name_Y2);
                hold on
                for no_group=2:num_X2
                    yyaxis right;
                    loglog(X1,Y1(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    yyaxis left;
                    loglog(X1,Y2(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                end
            else
            end
            xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
            set(gca,'FontSize',font_size)
            set(gca, 'LineWidth',gca_line_width)
            %             title([name_Y ' \rmat \rm' now_str]);
            grid on%显示网格
            yyaxis left;
            text(0.4*X1(1),0.5e5,'(b)','FontSize',font_size)
            
            % legend按默认顺序：先左边全部，再右边全部
            % 自定义legend
            legend_style_group = zeros(1,num_X2);
            legend_text_group=cell(1,num_X2);
            for no_group=1:num_X2
                legend_style_group(no_group) = plot(NaN,NaN,['s' group_color_list{no_group}],'MarkerSize',marker_size,'MarkerEdgeColor',group_color_list{no_group},'MarkerFaceColor',group_color_list{no_group});
                legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
            end
            L1=legend(legend_style_group, legend_text_group{:});
            set(L1,'FontSize',font_size);
            set(L1,'location','southeast');
            set(L1,'box','off')
            set(L1,'AutoUpdate','off')
     
            legend_style_group = zeros(1,2);
            legend_text_group={name_Y1,name_Y2};
            for no_group=1:2
                legend_style_group(no_group) = plot(NaN,NaN,[group_style_list{no_group} 'k'],'LineWidth',plot_line_width);
            end
            axes2 = axes('position',get(gca,'position'),'visible','off');
            L2=legend(axes2,legend_style_group, legend_text_group{:});
            set(L2,'FontSize',font_size);
            set(L2,'location','south');
            set(L2,'box','off')
            
            fprintf('%s = %.2e \n',[ num2str(X1(21)) ',15eV,σ= '],sigma_medium(21,3));
            fprintf('%s = %.2e \n',[ num2str(X1(11)) ',15eV,σ= '],sigma_medium(11,3));
            
                        
            name_Y='等离子体等效电磁参数比值';
            Y1=sigma_medium./(-w_RF*eps_medium);
            name_Y1='\it{\bf\sigma}/{\bf\omega}|{\bf\epsilon}|\rm';
            handle_fig=figure;
            flag_logy=true;
            %             flag_logy=false;
            if flag_logy
                for no_group=1:num_X2
                    loglog(X1,Y1(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    hold on
                end
                ylabel(name_Y1);
            else
            end
            xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
            set(gca,'FontSize',font_size)
            set(gca, 'LineWidth',gca_line_width)
            axis([X1(1),X1(end),1,2e1])
            line([X1(1),X1(end)],[1,1],'linestyle',':','linewidth',1*plot_line_width,'color','k');
            line([X1(1),X1(end)],[10,10],'linestyle',':','linewidth',1*plot_line_width,'color','k');
            grid on%显示网格
            text(0.4*X1(1),0.6e0,'(a)','FontSize',font_size)
            %             title([name_Y1 ' \rmat \rm' now_str]);
            legend_text_group=cell(1,num_X2);
            for no_group=1:num_X2
                legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
            end
            L1=legend(legend_text_group{:});
            set(L1,'FontSize',font_size);
            set(L1,'location','northwest');
            set(L1,'box','off')
            set(L1,'AutoUpdate','off')
            
            no_X2=no_mid_X2;
            if ~flag.data_map
                fprintf(['论文绘图用：' X_var{idx_X2} '=' num2str(X2(no_X2)) ' ' unit_X_var{idx_X2} '\n'])
            end
            
            % 碰撞频率
            %             name_Y='碰撞频率';
            %             handle_fig=figure;
            %             loglog(X1,venp(:,no_X2),'-r','LineWidth',plot_line_width);
            %             %                 axis equal
            %             %                 axis([X1(1),X1(end),1e4,1e8]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            %             hold on
            %             loglog(X1,veniz(:,no_X2),'--r','LineWidth',plot_line_width);
            %             loglog(X1,veip(:,no_X2),'-.r','LineWidth',plot_line_width);
            %             line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',0.5*plot_line_width,'color','k');
            %             loglog(X1,vm(:,no_X2),'-ob','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','b',...
            %                 'LineWidth',plot_line_width);
            %             loglog(X1,vst(:,no_X2),'-sb','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','b',...
            %                 'LineWidth',plot_line_width);
            %             L1=legend('\it{\bf\nu}\rm_{en}^{p}','\it{\bf\nu}\rm_{en}^{iz}','\it{\bf\nu}\rm_{ei}^{p}','{\it\bf\omega}','\it{\bf\nu}\rm_{m}','\it{\bf\nu}\rm_{st}');
            %
            %             set(L1,'FontSize',font_size);
            %             set(L1,'location','best');
            %             set(L1,'Orientation','horizontal');
            %             xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
            %             ylabel('\it{\bf\nu} \rm[Hz]');
            %             set(gca,'FontSize',font_size)
            %             set(gca, 'LineWidth',gca_line_width)
            %
            %             %         set(L1,'box','off')
            %             %             title(['碰撞频率 \rmat \rm' now_str]);
            %             grid on%显示网格
            
            name_Y='等效碰撞频率';
            handle_fig=figure;
            if flag_multi_Te
                for no_group=1:num_X2
                    loglog(X1,vm(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    hold on
                    loglog(X1,vst(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                end
                line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',1*plot_line_width,'color','k');
                
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('\it{\bf\nu} \rm[Hz]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                %             title(['碰撞频率 \rmat \rm' now_str]);
                grid on%显示网格
                axis([X1(1),X1(end),2e6,2e8])
                text(0.4*X1(1),1e6,'(a)','FontSize',font_size)
                
                % 自定义legend
                legend_style_group = zeros(1,num_X2);
                legend_text_group=cell(1,num_X2);
                for no_group=1:num_X2
                    legend_style_group(no_group) = plot(NaN,NaN,['s' group_color_list{no_group}],'MarkerSize',marker_size,'MarkerEdgeColor',group_color_list{no_group},'MarkerFaceColor',group_color_list{no_group});
                    legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
                end
                L1=legend(legend_style_group, legend_text_group{:});
                set(L1,'FontSize',font_size);
                set(L1,'location','northwest');
                set(L1,'box','off')
                %                 set(L1,'Orientation','horizontal');
                set(L1,'AutoUpdate','off')
                text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
                
                legend_style_group = zeros(1,2);
                legend_text_group={'\it{\bf\nu}\rm_{m}','\it{\bf\nu}\rm_{st}'};
                for no_group=1:2
                    legend_style_group(no_group) = plot(NaN,NaN,[group_style_list{no_group} 'k'],'LineWidth',plot_line_width);
                end
                axes2 = axes('position',get(gca,'position'),'visible','off');
                L2=legend(axes2,legend_style_group, legend_text_group{:});
                set(L2,'FontSize',font_size);
                set(L2,'location','north');
%                 set(L2,'Orientation','horizontal');
                set(L2,'box','off')
            else
            end
            
            flag_multi_Te=false;
            
            name_Y='碰撞频率';
            handle_fig=figure;
            if flag_multi_Te
                for no_group=1:num_X2
                    loglog(X1,venp(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    hold on
                    loglog(X1,veniz(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                    loglog(X1,veip(:,no_group),[group_style_list{3} group_color_list{no_group}],'LineWidth',plot_line_width);
                end
                line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',1*plot_line_width,'color','k');
                
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('\it{\bf\nu} \rm[Hz]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                %             title(['碰撞频率 \rmat \rm' now_str]);
                grid on%显示网格
                axis([X1(1),X1(end),1e4,2e8])
                
                % 自定义legend
                legend_style_group = zeros(1,num_X2);
                legend_text_group=cell(1,num_X2);
                for no_group=1:num_X2
                    legend_style_group(no_group) = plot(NaN,NaN,['s' group_color_list{no_group}],'MarkerSize',marker_size,'MarkerEdgeColor',group_color_list{no_group},'MarkerFaceColor',group_color_list{no_group});
                    legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
                end
                L1=legend(legend_style_group, legend_text_group{:});
                set(L1,'FontSize',font_size);
                set(L1,'location','north');
                set(L1,'box','off')
                set(L1,'Orientation','horizontal');
                set(L1,'AutoUpdate','off')
                text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
                
                legend_style_group = zeros(1,3);
                legend_text_group={'\it{\bf\nu}\rm_{en}^{(p)}','\it{\bf\nu}\rm_{en}^{(iz)}','\it{\bf\nu}\rm_{ei}^{(p)}'};
                for no_group=1:3
                    legend_style_group(no_group) = plot(NaN,NaN,[group_style_list{no_group} 'k'],'LineWidth',plot_line_width);
                end
                axes2 = axes('position',get(gca,'position'),'visible','off');
                L2=legend(axes2,legend_style_group, legend_text_group{:});
                set(L2,'FontSize',font_size);
                set(L2,'location','north');
                set(L2,'Orientation','horizontal');
                set(L2,'box','off')
            else
                loglog(X1,venp(:,no_X2),'-r','LineWidth',plot_line_width);
                %                 axis equal
                axis([X1(1),X1(end),1e4,5e7]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
                hold on
                loglog(X1,veniz(:,no_X2),'--r','LineWidth',plot_line_width);
                loglog(X1,veip(:,no_X2),'-.b','LineWidth',plot_line_width);
                line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',plot_line_width,'color','k');
                L1=legend('\it{\bf\nu}\rm_{en}^{(p)}','\it{\bf\nu}\rm_{en}^{(iz)}','\it{\bf\nu}\rm_{ei}^{(p)}');
                text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
                set(L1,'FontSize',font_size);
                set(L1,'location','southeast');
                set(L1,'Orientation','horizontal');
                set(L1,'box','off');
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('\it{\bf\nu} \rm[Hz]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                
                %         set(L1,'box','off')
                %             title(['碰撞频率 \rmat \rm' now_str]);
                grid on%显示网格
                text(0.4*X1(1),0.2e4,'(b)','FontSize',font_size)
            end            
            %计算线性性
            fprintf('%s = %.2e \n','15eV,k= ',(log(veip(end,no_X2))-log(veip(1,no_X2)))/(log(X1(end))-log(X1(1))));
            fprintf('%s = %.2e \n','15eV,k= ',(log(veip(11,no_X2))-log(veip(1,no_X2)))/(log(X1(11))-log(X1(1))));
            fprintf('%s = %.2e \n','15eV,k= ',(log(veip(21,no_X2))-log(veip(11,no_X2)))/(log(X1(21))-log(X1(11))));
            fprintf('%s = %.2e \n','15eV,k= ',(log(veip(31,no_X2))-log(veip(21,no_X2)))/(log(X1(31))-log(X1(21))));
            
            % 碰撞截面随能量变化
%             [E_enp,a_enp]=textread('ELASTIC.txt','%n %n','headerlines',10);
%             [E_eniz,a_eniz]=textread('IONIZATION.txt','%n %n','headerlines',10);
%             handle_fig=figure;
%             loglog(E_enp,a_enp,'-r','LineWidth',plot_line_width)
%             hold on
%             loglog(E_eniz,a_eniz,'-.b','LineWidth',plot_line_width)
%             % line([X1(1),X1(end)],[1,1],'linestyle',':','linewidth',1*plot_line_width,'color','k');
%             
%             ylabel('Cross section [m^2]');
%             xlabel('Energy [eV]')
%             set(gca,'FontSize',font_size)
%             set(gca, 'LineWidth',gca_line_width)
%             %             title([name_Y ' \rmat \rm' now_str]);
%             grid on%显示网格
%             % text(0.4*X1(1),0.5e5,'(b)','FontSize',font_size)
%             
%             L1=legend('{\it\bf\sigma}_{en}^{(p)}','{\it\bf\sigma}_{en}^{(iz)}');
%             set(L1,'FontSize',font_size);
%             set(L1,'location','northwest');
%             set(L1,'box','off')
%             set(L1,'AutoUpdate','off')
            
            % 特征尺寸
            handle_fig=figure;
            if flag_multi_Te
                for no_group=1:num_X2
                    loglog(X1,skin_depth_eff(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
                    hold on
                    loglog(X1,delta_st(:,no_group),[group_style_list{3} group_color_list{no_group}],'LineWidth',plot_line_width);
                    loglog(X1,wavelength_wave(:,no_group),[group_style_list{2} group_color_list{no_group}],'LineWidth',plot_line_width);
                end
                line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle',':','linewidth',0.5*plot_line_width,'color','k');
                
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('Characteristic length [m]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                %             title(['特征尺寸 \rmat \rm' now_str]);
                grid on%显示网格
                axis([X1(1),X1(end),5e-3,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
                
                % 自定义legend
                legend_style_group = zeros(1,num_X2);
                legend_text_group=cell(1,num_X2);
                for no_group=1:num_X2
                    legend_style_group(no_group) = plot(NaN,NaN,['s' group_color_list{no_group}],'MarkerSize',marker_size,'MarkerEdgeColor',group_color_list{no_group},'MarkerFaceColor',group_color_list{no_group});
                    legend_text_group{no_group}=[name_X_var{idx_X2} '=' num2str(X2(no_group)) ' eV'];
                end
                L1=legend(legend_style_group, legend_text_group{:});
                set(L1,'FontSize',font_size);
                set(L1,'location','north');
                set(L1,'box','off')
                set(L1,'Orientation','horizontal');
                set(L1,'AutoUpdate','off')
                text(1.1*X1(1),1.5*r_plasma,'\itr\rm_{plasma}','FontSize',font_size)
                
                legend_style_group = zeros(1,3);
                legend_text_group={'\it{\bf\delta}_{\rmeff}','\it{\bf\lambda}','\it{\bf\delta}_{\rmst}'};
                for no_group=1:3
                    legend_style_group(no_group) = plot(NaN,NaN,[group_style_list{no_group} 'k'],'LineWidth',plot_line_width);
                end
                axes2 = axes('position',get(gca,'position'),'visible','off');
                L2=legend(axes2,legend_style_group, legend_text_group{:});
                set(L2,'FontSize',font_size);
                set(L2,'location','north');
                set(L2,'Orientation','horizontal');
                set(L2,'box','off')
            else
                loglog(X1,skin_depth_eff(:,no_X2),'-r','LineWidth',plot_line_width);
                hold on
                loglog(X1,delta_st(:,no_X2),'--r','LineWidth',plot_line_width);
                loglog(X1,delta_m(:,no_X2),'-.r','LineWidth',plot_line_width);
                loglog(X1,wavelength_wave(:,no_X2),'-ob','LineWidth',plot_line_width);
                line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle',':','linewidth',plot_line_width,'color','k');
                %                     line([X1(1),5e17],[0.05,0.05],'linestyle','-.','linewidth',1.5,'color','k');
                xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                ylabel('Characteristic length [m]');
                set(gca,'FontSize',font_size)
                set(gca, 'LineWidth',gca_line_width)
                %                 L1=legend('\it{\bf\delta}_{\rmeff}','\it{\bf\delta}_{\rmst}','\it{\bf\lambda}','\itr\rm_{plasma}');
%                 L1=legend('\it{\bf\delta}_{\rmeff}','\it{\bf\delta}_{\rmst}','\it{\bf\lambda}');
                L1=legend('\it{\bf\delta}_{\rmeff}','\it{\bf\delta}_{\rmst}','\it{\bf\delta}_{\rmm}','\it{\bf\lambda}');
                text(1.1*X1(1),1.5*r_plasma,'\itr\rm_{plasma}','FontSize',font_size)
                set(L1,'FontSize',font_size);
                set(L1,'location','northeast');
                %                 set(L1,'Orientation','horizontal');
                set(L1,'box','off')
                %             title(['特征尺寸 \rmat \rm' now_str]);
                grid on%显示网格
                axis([X1(1),X1(end),5e-3,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            end
            
            % 特征频率
            name_Y='特征频率';
            handle_fig=figure;
            loglog(X1,wpi(:,no_X2),'-r','LineWidth',plot_line_width);
            %             axis equal
            %             axis([X1(1),X1(end),1e6,1e12]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            hold on
            loglog(X1,wpe(:,no_X2),'--b','LineWidth',plot_line_width);
            xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
            ylabel('\it{\bf\omega} \rm[Hz]');
            set(gca,'FontSize',font_size)
            set(gca, 'LineWidth',gca_line_width)
            line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',plot_line_width,'color','k');
            %             L1=legend('\it{\bf\omega}_{\rmpi}','\it{\bf\omega}_{\rmpe}','{\it\bf\omega}');
            L1=legend('\it{\bf\omega}_{\rmpi}','\it{\bf\omega}_{\rmpe}');
            text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
            set(L1,'FontSize',font_size);
            set(L1,'location','northwest');
            set(L1,'box','off')
            %             title(['特征频率 \rmat \rm' now_str]);
            grid on%显示网格
            
            % vst简化表达式与Cazzador fit表达式对比
            switch flag.stoc_model
                case {'Vahedi-simplify','Cazzador-simplify'}
                    name_Y='对比st表达式';
                    handle_fig=figure;
                    set(handle_fig,'defaultAxesColorOrder',[left_color; right_color]);
                    % marker会遮掩线型，而且也不好看
                    flag_logy=true;
                    %             flag_logy=false;
                    if flag_logy
                        yyaxis left;
                        h1=loglog(X1,vst_fit(:,no_X2),'-dr','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','r',...
                            'LineWidth',plot_line_width);
                        yyaxis right;
                        h2=loglog(X1,delta_st_fit(:,no_X2),'-b','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','r',...
                            'LineWidth',plot_line_width);
                        hold on
                        yyaxis left;
                        h3=loglog(X1,vst(:,no_X2),'-sr','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','r',...
                            'LineWidth',plot_line_width);
                        line([X1(1),X1(end)],[w_RF,w_RF],'linestyle',':','linewidth',plot_line_width,'color','r');
                        axis([X1(1),X1(end),4e6,1e8])
                        ylabel('{\it\bf\nu}_{st} [Hz]');
                        yyaxis right;
                        h4=loglog(X1,delta_st(:,no_X2),'-ob','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','b',...
                            'LineWidth',plot_line_width);
                        line([X1(1),X1(end)],[r_plasma,r_plasma],'linestyle',':','linewidth',plot_line_width,'color','b');
                        axis([X1(1),X1(end),1e-3,3e-1])
                        ylabel('{\it\bf\delta}_{st} [m]');
                    else
                    end
                    xlabel([name_X_var{idx_X1} ' \rm[' unit_X_var{idx_X1} ']']);
                    set(gca,'FontSize',font_size)
                    set(gca, 'LineWidth',gca_line_width)
                    grid on%显示网格
                    L1=legend([h1,h3],'{\it\bf\nu}_{st}, numerical','{\it\bf\nu}_{st}, piecewise');
                    set(L1,'FontSize',font_size);
                    set(L1,'location','southwest');
                    set(L1,'Orientation','horizontal');
                    set(L1,'box','off')
                    
                    yyaxis left
                    text(0.7*X1(1),w_RF,'{\it\bf\omega}','FontSize',font_size)
                    yyaxis right
                     text(1.05*X1(end),1.2*r_plasma,'\itr\rm_{plasma}','FontSize',font_size)
                    
                    axes2 = axes('position',get(gca,'position'),'visible','off');
                    L2=legend(axes2,[h2,h4], '{\it\bf\delta}_{st}, numerical','{\it\bf\delta}_{st}, piecewise');
                    set(L2,'FontSize',font_size);
                    set(L2,'location','northeast');
                    set(L2,'box','off')
                    set(L2,'Orientation','horizontal');
                   
                    
                    no_X1=11;
                    (vst(no_X1)-vst_fit(no_X1))/(vst_fit(no_X1))
                    
                    
            end
            
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
                %         loglog(X1,vm,'-r','LineWidth',1.5);
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
                Y1=vm;
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