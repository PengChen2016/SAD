% output for paper            

if flag.output_for_paper
    flag.input_plasma='ELISE_base';
    
    flag.data_map=false;
    flag.using_stored_data=false;
    
    flag.stoc_model='Cazzador-fit';
            flag.stoc_model='Vahedi-simplify'; %用于对比st模型
    flag.electric_model='transformer_base';
    flag.good_conductor_approximation=false;
    
%     flag.input_plasma='CHARLIE_base';
%     flag.sweep=false; 
end


%% 输出及可视化
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
            %             loglog(X1,plasma.nu_m(:,no_X2),'-ob','MarkerIndices',marker_indices,'MarkerSize',marker_size,'MarkerEdgeColor','b',...
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
                    loglog(X1,plasma.nu_m(:,no_group),[group_style_list{1} group_color_list{no_group}],'LineWidth',plot_line_width);
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