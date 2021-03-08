% RF-ICP源功率耦合模型
% 主要用于分析负源激励器功率耦合，
% 计算等效阻抗与RF传输效率
% MATLAB R2017a

%% 初始化
close all
clear

addpath(genpath('./packages')) 

% 如无说明，则单位为国际单位制
% 全局变量
constants=get_constants();% 全局常数 结构体

now_str=datestr(now,'yyyy-mm-dd HH:MM:SS');
% now_str=datestr(now,'yyyy-mm-dd');
disp(now_str)

%% 输入
%%%% 控制位

% flag.input_plasma='2018Jainb_ELISE_typical';
% flag.input_plasma='2019Raunera_CHARLIE_sweep';
% flag.input_plasma='BATMAN_typical';
% flag.input_plasma='2021Zielke_BATMAN_sweep';
% flag.input_plasma='small_source1_LZS';
% flag.input_plasma='2018Jainb_NIO1_sweep';
% flag.input_plasma='CHARLIE_10Pa_4MHz_520W';
flag.input_plasma='2020Chen_NIS_sweep1';
% flag.input_plasma='2020Chen_NIS_sweep_p';

% flag.input_geometry='ELISE_base'; 
% flag.input_geometry='BATMAN_base'; %no data
% flag.input_geometry='HUST_small_driver_base'; %by ZP
% flag.input_geometry='CHARLIE_base'; 
% flag.input_geometry='small_source1_LZS'; %Ref from 李增山-整体模型耦合解析电模型-2020.03.30\ by CP
% flag.input_geometry='NIO1_base'; %no data
flag.input_geometry='HUST_large_driver_base'; %by CP, LJW 201029

flag.Rmetal='measured-Rmetal-woplasma';
% flag.Rmetal=calculated-Rcoil-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';

% flag.electric_model='';
flag.electric_model='transformer_base';
% flag.electric_model='analytical_base';
% flag.electric_model='transformer_2011Chabert';

% stoc表达式
% flag.stoc_model='';
% flag.stoc_model='Vahedi-simplify';
% flag.stoc_model='Cazzador-simplify';
flag.stoc_model='Cazzador-fit';


flag.sweep=false;
flag.sweep=true; % 参数扫描则为真

flag.data_map=false;
flag_using_stored_data=false;


flag_good_conductor_approximation=false;
% flag_good_conductor_approximation=true;

flag_output_plasma_model=true;
flag_output_plasma_model=false;
flag_output_electric_model=true;
% flag_output_electric_model=false;
flag_output_for_paper=true;
flag_output_for_paper=false;
if flag_output_for_paper
    flag.input_plasma='ELISE_base';
    flag.sweep=true; 
%     flag.sweep=false;
    
    flag.data_map=false;
    flag_using_stored_data=false;
    
    flag.stoc_model='Cazzador-fit';
            flag.stoc_model='Vahedi-simplify'; %用于对比st模型
    flag_electric_model='transformer_base';
    flag_good_conductor_approximation=false;
    
%     flag.input_plasma='CHARLIE_base';
%     flag.sweep=false; 
end

input=get_input_data( flag );

    %% ICP的等效电磁媒质模型 equivalent_medium_model_of_plasma
    %%%% 如果改成函数，虽然简洁不需要多处修改了，但传递数据似乎比较麻烦？不能改成函数，顶多多脚本文件
    % 主要基于2014Cazzador、1995Vahedia，借鉴2018Jain、

file_name='stored_data_test200529_2.mat';
% 已存储则注释掉下列语句
if flag.sweep && flag_using_stored_data
    warning('Store data first.')
    error('end')
end
% 在flag_using_stored_data=false时运行以下三行一次,然后将flag_using_stored_data改回true
% if flag.sweep
%     flag_using_stored_data=true;
%     save(file_name)
%     error('end')
% end

if flag.sweep && flag_using_stored_data
    load(file_name)
    warning('using stored data')
end
%%%%%%%%%%%%%%%%%%%%%%%%%% 后处理 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_output_plasma_model
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

% error('end')

%% 螺线管线圈ICP源的电模型
% 改成函数，变量传递是可以接受的。但作为一个实验室代码，可能这样更便于阅读――做好版本管理就可以了
switch flag_electric_model
    case 'transformer_base'
        % 变压器模型
        disp('## ICP源电模型：变压器模型')
        if l_plasma<r_plasma
            warn('plasma为短螺线管')
        end
        % fprintf('%s = %.2e \n','IC的长冈系数',check_Nagaoka(2*r_plasma/l_plasma));
        Lmp=constants.mu0*pi*r_plasma^2*check_Nagaoka(2*r_plasma/l_plasma)/l_plasma;  %等离子体感应电感-螺线管电感经验公式。相较于ZP，考虑了长冈系数
        % Lmp1=constants.mu0*pi*r_plasma^2*(1-8*r_plasma/(3*pi*l_plasma))/l_plasma; % 某论文
        % (1-8*r_plasma/(3*pi*l_plasma))与长冈系数结果有显著差别。长冈系数计算结果更接近2018Jain
        for X1i=1:num_X1%不同的ne
            for X2i=1:num_X2%不同的Te
                Rp(X1i,X2i)=pi*r_plasma/(skin_depth_eff(X1i,X2i)*sigmap_real(X1i,X2i)*l_plasma); %二次侧等离子体电阻，2011Chabert
                Lp(X1i,X2i)=Rp(X1i,X2i)/veff(X1i,X2i); %等离子体中电子惯性形成的电感。
                if flag_good_conductor_approximation
                    % 这似乎并不是忽略复电导率虚部
                    disp('变压器模型中使用良导体近似')
                    Rp(X1i,X2i)=2*pi*r_plasma/(skin_depth_eff(X1i,X2i)*sigmap_real(X1i,X2i)*l_plasma); %二次侧等离子体电阻，2011Chabert
                    Lp(X1i,X2i)=0;
                end
                
                M=sqrt(Lcoil*(Lmp+Lp(X1i,X2i)))*r_plasma^2/r_coil^2; %线圈和等离子体之间的互感
                %         M=(N-0)*Lmp;          %by ZP,但该表达式未考虑线圈半径
                %         fprintf('%s = %.2e , %s = %.2e , %s = %.2e\n','1995Vehadi,Rp',Rp,'Lmp',Lmp,'M',M);
                
                %         %2018Jainb中
                %         if flag_single_independent_variable
                %             disp('使用2018Jainb的Rp、Lmp、M表达式')
                %             Rp(X1i,X2i)=2*pi*r_plasma/(skin_depth_eff(X1i,X2i)*sigmaeff_dc(X1i,X2i)*l_plasma);
                %             Lmp=0.002*pi*(2*r_plasma*100)*(log(4*2*r_plasma/l_plasma)-0.5)*1e-6; %该表达式结果可能为负值，则bug;即使取绝对值，也超出适用范围
                %             if Lmp<0
                %                 warning('Lmp<0. We will use Lmp=-Lmp.')
                %                 pause
                %                 Lmp=-Lmp;
                %             end
                %             M=0.0095*N*1e-6*(2*r_plasma*100)^2/sqrt((2*r_coil*100)^2+(l_coil*100)^2);
                %             fprintf('%s = %.2e , %s = %.2e , %s = %.2e\n','2018Jainb,Rp',Rp,'Lmp',Lmp,'M',M);
                %         end
                
                PER(X1i,X2i)=Rp(X1i,X2i)*w_RF^2*M^2/(Rp(X1i,X2i)^2+w_RF^2*(Lmp+Lp(X1i,X2i))^2);  %换算到一次侧的等离子体等效电阻
                Rsys(X1i,X2i)=PER(X1i,X2i)+Rmetal(X1i,X2i); %系统等效阻抗的电阻分量
                Lplasma(X1i,X2i)=-M^2*w_RF^2*(Lmp+Lp(X1i,X2i))/(Rp(X1i,X2i)^2+w_RF^2*(Lmp+Lp(X1i,X2i))^2); %换算到一次侧的等离子体等效电感
                if ~flag.sweep && Lplasma(X1i,X2i)>0
                    warning('Lplasma应<0, please stop and check')
                    pause
                end
                Xplasma(X1i,X2i)=w_RF*Lplasma(X1i,X2i);
                Lsys(X1i,X2i)=Lcoil+Lplasma(X1i,X2i);   %系统等效阻抗的电感分量
                Xsys(X1i,X2i)=w_RF*Lsys(X1i,X2i); %系统等效阻抗的电抗分量
                P_abs(X1i,X2i)=PER(X1i,X2i)*Icoil_rms(X1i,X2i)^2;            %吸收功率 即Pcp_real
            end
        end
    case 'analytical_base'
        % 解析模型 2011Chabert
        disp('## ICP源电模型：轴对称平行平面解析模型')
        % 几何校正
        l_equ=l_chamber; %等离子体、线圈均与腔室同长
        %         N_equ=N_coil;
        N_equ=N_coil*l_equ/l_coil; %变换线圈长度
        r_p=r_chamber; %等离子体与腔室同半径
        
        if flag.sweep
            
        else
            % 增大线圈阻抗的校正方法
            %             % 理论计算线圈阻抗
            %             fprintf('几何校正前线圈阻抗：')
            %             fprintf('%s = %.2e , ','Rcoil',Rcoil);
            %             fprintf('%s = %.2e \n','Lcoil',Lcoil);
            %             Rcoil_modified=N_equ*2*pi*r_coil*sqrt(constants.mu0*w_RF*constants.sigma_Cu/2)/2/pi/r_wire_of_coil/constants.sigma_Cu;%已考虑集肤效应，未考虑邻近效应
            %             Lcoil_modified=constants.mu0*check_Nagaoka(2*r_coil/l_equ)*pi*r_coil^2*N_equ^2/l_equ;  %一次螺线管线圈电感理论计算表达式
            %             fprintf('几何校正后，理论计算线圈阻抗：')
            %             fprintf('%s = %.2e , ','Rcoil',Rcoil_modified);
            %             fprintf('%s = %.2e \n','Lcoil',Lcoil_modified);
            %             if Rcoil<Rcoil_modified
            %                 Rcoil=Rcoil_modified*ones(size(Rcoil));
            %                 fprintf('使用几何校正后')
            %             else
            %                 fprintf('使用几何校正前')
            %             end
            %             fprintf('%s = %.2e \n','Rcoil',Rcoil);
            %             if abs(Lcoil)<abs(Lcoil_modified)
            %                 Lcoil=Lcoil_modified;
            %                 fprintf('使用几何校正后')
            %             else
            %                 fprintf('使用几何校正前')
            %             end
            %             fprintf('%s = %.2e \n','Lcoil',Lcoil);
            %             Qcoil=w_RF*Lcoil/Rcoil;
        end
        
        % 介质窗相对介电常数，玻璃约2，氧化铝陶瓷约9
        epst_r=9;
        k_t_wave=w_RF*sqrt(constants.mu0*constants.eps0*epst_r);
        if k_t_wave*r_p>1
            warn('介质管中磁场不是常数，该模型不适用')
        end
        % 重构：将以上参数代码提到最前面输入部分去
        Hz0=N_equ*Icoil_rms(X1i,X2i)/l_equ;
        for X1i=1:num_X1%不同的ne
            for X2i=1:num_X2%不同的Te
                Hz_p=@(r)Hz0*besselj(0,k_p_wave(X1i,X2i)*r)/besselj(0,k_p_wave(X1i,X2i)*r_p);
                % 与2018Zhao中带FS的FEM模型的20A时磁场结果相比，小了
                Etheta_p=@(r)-1i*k_p_wave(X1i,X2i)*Hz0*besselj(1,k_p_wave(X1i,X2i)*r)/...
                    (w_RF*constants.eps0*epsp_r_real(X1i,X2i)+besselj(0,k_p_wave(X1i,X2i)*r_p));
                % 与2018Zhao中带FS的FEM模型的20A时电场结果相比，大了
                Hz_t=@(r)Hz0; %介质管中磁场均匀
                
                %检验表达式用-断点调试
                % 复坡印廷定理计算复功率
                Etheta_t_at_rc=Etheta_p(r_p)*r_p/r_coil-1i*w_RF*constants.mu0*Hz0*(r_coil^2-r_p^2)/(2*r_coil); %rc处电场
                Pcp1=-pi*r_coil*l_equ*Etheta_t_at_rc*Hz_t(r_coil);
                % ERROR 该表达式与简化表达式结果不一致
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %简化表达式计算复功率
                Pcp_part1=k_p_wave(X1i,X2i)*r_p*besselj(1,k_p_wave(X1i,X2i)*r_p)/...
                    (w_RF*constants.eps0*epsp_r(X1i,X2i)*besselj(0,k_p_wave(X1i,X2i)*r_p));
                Pcp_part2=Pcp_part1+w_RF*constants.mu0*(r_coil^2-r_p^2)/2;
                Pcp=1i*pi*N_equ^2*Icoil_rms(X1i,X2i)^2*Pcp_part2/l_equ;
                
                Pcp_real=real(Pcp); %即Pabs
                Pabs(X1i,X2i)=Pcp_real;
                %                 %检验表达式用-断点调试
                %                 Pcp_real1_part1=1i*k_p_wave(X1i,X2i)*r_p*besselj(1,k_p_wave(X1i,X2i)*r_p)/(w_RF*constants.eps0*epsp_r(X1i,X2i)*besselj(0,k_p_wave(X1i,X2i)*r_p));
                %                 Pcp_real1=pi*N_equ^2*I_coil^2*real(Pcp_real1_part1)/l_equ; %与Pcp_real一致，说明上下表达式自洽
                PER(X1i,X2i)=2*Pcp_real/Icoil_rms(X1i,X2i)^2; %换算到一次侧的等离子体等效电阻
                I_p=l_equ*Hz0*(1/besselj(0,k_p_wave(X1i,X2i)*r_p)-1); %等离子体中电流
                Rp(X1i,X2i)=2*Pcp_real/abs(I_p)^2; %二次侧等离子体电阻
                %                 %检验表达式用-断点调试
                % 与LZS实现结果一致
                %                 Rp1_part1=1i*k_p_wave(X1i,X2i)*r_p*besselj(1,k_p_wave(X1i,X2i)*r_p)/(w_RF*constants.eps0*epsp_r(X1i,X2i)*besselj(0,k_p_wave(X1i,X2i)*r_p));
                % 可能少了一个N^2
                %                 Rp1=2*pi*real(Rp1_part1)/l_equ/(1/besselj(0,k_p_wave(X1i,X2i)*r_p)-1)^2; %与Rp一致，说明上下表达式自洽
                Rsys(X1i,X2i)=PER(X1i,X2i)+Rmetal(X1i,X2i); %系统等效阻抗的电阻分量
                Xsys(X1i,X2i)=2*imag(Pcp)/Icoil_rms(X1i,X2i)^2; %系统等效阻抗的电抗分量
                Lsys(X1i,X2i)=Xsys(X1i,X2i)/w_RF; %系统等效阻抗的电感分量
                Lplasma(X1i,X2i)=Lsys(X1i,X2i)-Lcoil;  % 一次侧的等离子体电感（定义的）
                if ~flag.sweep && Lplasma(X1i,X2i)>0
                    warning('Lplasma应<0, please stop and check')
                    pause
                end
                Xplasma(X1i,X2i)=w_RF*Lplasma(X1i,X2i); %一次侧的等离子体导致电抗
            end
        end
    case 'transformer_2011Chabert'
        warning('施工中. Please stop and check.')
    otherwise
        warning('Unexpected electric model. Please stop and check.')
        pause
end
PTE=PER./Rsys;   %射频功率传输效率 PTE
PCF=Rsys./Xsys;  %激励器射频功率耦合因数

%%%%%%%%%%%%%%%%%%%%%%%%%%% 后处理 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_output_electric_model
    if ~flag.sweep
        disp('等离子体等效阻抗')
        fprintf('%s = %.2e Ω, ','PER',PER);
        fprintf('%s = %.2e , ','Lplasma',Lplasma);
        fprintf('%s = %.2e \n','Xplasma',Xplasma);
        disp('系统等效阻抗Zsys')
        fprintf('%s = %.2e , ','Rsys',Rsys);
        fprintf('%s = %.2e , ','Lsys',Lsys);
        fprintf('%s = %.2e \n','Xsys',Xsys);
        fprintf('%s = %.2e \n','Rmetal',Rmetal);
        fprintf('%s = %.2e \n','射频功率传输效率 PTE',PTE);
        fprintf('%s = %.2e \n','射频功率耦合因数 PCF',PCF);
    else
        if flag.data_map
            % 等离子体等效阻抗
            figure
            semilogx(X1,PER,'-r','LineWidth',1.5);
            hold on
            semilogx(X1,-Xplasma,'-b','LineWidth',1.5);
            ylabel('\itZ_{\rmeff}');
            xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            title(['等离子体等效阻抗 \rmat \rm' now_str]);
            L1=legend('\itPER','\it-X_{\rmplasma}');
            set(L1,'FontSize',10);
            %                     set(L1,'location','southeast');
            %         set(L1,'box','off')
            %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            grid on%显示网格
            %             fig_background_transparent(gcf,gca)  %临时使用
            
            % 系统等效阻抗
            figure
            semilogx(X1,Rsys,'-r','LineWidth',1.5);
            hold on
            semilogx(X1,Xsys,'-b','LineWidth',1.5);
            ylabel('\itZ_{\rmeff}');
            xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            title(['系统等效阻抗 \rmat \rm' now_str]);
            L1=legend('\itR_{\rmsys}','\itX_{\rmsys}');
            set(L1,'FontSize',10);
            %                     set(L1,'location','southeast');
            %         set(L1,'box','off')
            %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            grid on%显示网格
            %             fig_background_transparent(gcf,gca)  %临时使用
            
            % PTE与PCF
            figure
            yyaxis left
            semilogx(X1,PTE,'-r','LineWidth',1.5);
            ylabel('\itPTE');
            yyaxis right
            semilogx(X1,PCF,'-b','LineWidth',1.5);
            ylabel('\itPCF');
            xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
            title(['PTE&PCF \rmat \rm' now_str]);
            L1=legend('\itPTE','\itPCF');
            set(L1,'FontSize',10);
            %                     set(L1,'location','southeast');
            %         set(L1,'box','off')
            %         axis([X1(1),X1(end),1e4,1e10]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
            grid on%显示网格
            %             fig_background_transparent(gcf,gca)  %临时使用
            
        else
            name_Y='等离子体等效阻抗\itZ_{\rmplasma}';
            Y1=PER;
            name_Y1='\itPER\rm, namely \itR_{\rmplasma}\rm[\Omega]';
            Y2=Xplasma;
            name_Y2='\itX_{\rmplasma}\rm[\Omega]';
            handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            line([X1(1),X1(end)],[0,0],'linestyle','-.','linewidth',1.5,'color','k');
            L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
                'Z=0');
            set(L1,'location','southwest');
            
            Y2=-Lplasma*1e6;
            name_Y2='-\itL_{\rmplasma}\rm[\muH]';
            handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            pause
            
            name_Y='系统等效阻抗\itZ_{\rmsys}即实验可测的端口阻抗';
            Y1=Rsys;
            name_Y1='\itR_{\rmsys}\rm[\Omega]';
            Y2=Xsys;
            name_Y2='\itX_{\rmsys}\rm[\Omega]';
            handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            Y2=Lsys*1e6;
            name_Y2='\itL_{\rmsys}\rm[\muH]';
            handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            pause
            

            
            
            name_Y='系统功率传输效率PTE和耦合因数PCF';
            Y1=PTE;
            name_Y1='\itPTE';
            Y2=PCF;
            name_Y2='\itPCF';
            handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
                X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
            L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
                [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
            set(L1,'location','northwest');
            pause
            
            %     figure(2)
            %     [constants.c,h]=contour(ne,Te,PTE,'LevelList',[0.7:0.05:1]);
            %     clabel(constants.c,h,'FontSize',11);
            %     title('Power tranfer efficinecy, at \itT_{\rmgas}\rm=400 K，\itp_{gas}\rm=0.3 Pa');
            %     xlabel('\itn_{\rme} \rm[m^{-3}]');
            %     ylabel('\itT_{\rme} \rm[eV]');
            %     set(gca,'XScale','log');
            %     grid on
            %
            %     figure(3)
            %     plot(Te,PER(11,:),'-g','LineWidth',1.5);
            %     hold on
            %     plot(Te,PER(21,:),'-b','LineWidth',1.5);
            %     hold on
            %     plot(Te,PER(31,:),'-k','LineWidth',1.5);
            %     hold on
            %     grid on
            %     L1=legend('\itn_{constants.e}\rm=1e16 m^{-3}','\itn_{constants.e}\rm=1e17 m^{-3}','\itn_{constants.e}\rm=1e18 m^{-3}');
            %     set(L1,'FontSize',10);
            %     title('Plasma equivalent resitance,at at \itT_{\rmgas}\rm=400 K，\itp_{\rmgas}\rm=0.3 Pa');
            %     xlabel('\itT_{\rme} \rm[eV]');
            %     ylabel('Plasma equivalent resitance,\itR_{\rmeq} \rm[\Omega]');
        end
    end
end

%% aid function
% 使用函数，主要是方便测试和复用。
% 弃用函数，在main中直接写表达式。因为在main中只使用一处，且基本不会复用，而写函数传参麻烦。

%美观目的的微调在origin中进行
function handle_fig=plot_parametric_1Y(Y1,name_Y1,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    semilogx(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,end),'-g','LineWidth',1.5);
else
    plot(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,end),'-g','LineWidth',1.5);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']); %\it斜体 \rm 正体 ^上标 _下标 {}组合
ylabel(name_Y1);
L1=legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
set(L1,'FontSize',10);
title([name_Y1 ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1与Y2在单纵轴下随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    semilogx(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,end),'-g','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,end),'--g','LineWidth',1.5);
else
    plot(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,end),'-g','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,end),'--g','LineWidth',1.5);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
ylabel([name_Y1 ' & ' name_Y2]);
L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
set(L1,'FontSize',10);
title([name_Y ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1与Y2双纵轴随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    yyaxis left;
    semilogx(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogx(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    semilogx(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    semilogx(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
else
    yyaxis left;
    plot(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    plot(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    plot(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    plot(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
set(L1,'FontSize',10);
title([name_Y ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_1Ylog(Y1,name_Y1,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    loglog(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,end),'-g','LineWidth',1.5);
else
    semilogy(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,end),'-g','LineWidth',1.5);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']); %\it斜体 \rm 正体 ^上标 _下标 {}组合
ylabel(name_Y1);
L1=legend(legend_X2{1},legend_X2{no_mid_X2},legend_X2{end});
set(L1,'FontSize',10);
title([name_Y1 ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_2Ylog(name_Y,Y1,name_Y1,Y2,name_Y2,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1与Y2在单对数纵轴下随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    loglog(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,end),'-g','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,end),'--g','LineWidth',1.5);
else
    semilogy(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,end),'-g','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,end),'--g','LineWidth',1.5);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
ylabel([name_Y1 ' & ' name_Y2]);
L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
set(L1,'FontSize',10);
title([name_Y ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function handle_fig=plot_parametric_2Ylogaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
    X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str)
% 首末及正中这三个代表性X2下，Y1与Y2双纵轴随X1的变化
handle_fig=figure;
if strcmp(X_var(idx_X1),'ne')||strcmp(X_var(idx_X1),'p')
    yyaxis left;
    loglog(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    loglog(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    loglog(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    loglog(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
else
    yyaxis left;
    semilogy(X1,Y1(:,1),'-r','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,no_mid_X2),'-b','LineWidth',1.5);
    hold on
    semilogy(X1,Y1(:,end),'-g','LineWidth',1.5);
    ylabel(name_Y1);
    yyaxis right;
    semilogy(X1,Y2(:,1),'--r','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,no_mid_X2),'--b','LineWidth',1.5);
    hold on
    semilogy(X1,Y2(:,end),'--g','LineWidth',1.5);
    ylabel(name_Y2);
end
xlabel([name_X_var{idx_X1} '\rm[' unit_X_var{idx_X1} ']']);
L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
    [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
set(L1,'FontSize',10);
title([name_Y ' \rmat \rm' now_str]);
%         axis([X1(1),X1(end),0,1]) %绘图显示范围，即[xmin,xmax,ymin,ymax]
grid on%显示网格
end

function fig_background_transparent(gcf,gca)
% 使图片背景透明，用于复制到文档
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
L=legend('show');
set(L,'box','off')
end

% skin_depth(5.8e7,50,1)=0.0093 经测试，skin_depth可用
function skin_d=skin_depth(sigma,f,mu_r)
% calculate skin depth of metal
constants.mu0=4*pi*1e-7;         %真空磁导率
skin_d=sqrt(1/(pi*f*mu_r*constants.mu0*sigma));
end

% lambda(1,1e6,1)=300 经测试，lambda可用
function lambda_d=lambda(eps_r,f,mu_r)
% calculate wave length in low-loss medium or lossless medium
constants.mu0=4*pi*1e-7;         %真空磁导率
constants.eps0=8.85e-12;          %真空介电常数
medium_v=1/sqrt(constants.mu0*mu_r*constants.eps0*abs(eps_r));
lambda_d=medium_v/f;
end