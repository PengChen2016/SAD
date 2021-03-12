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
% flag.input_plasma='given_directly'; %不能用于get_input_data()

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

flag.medium_approximation='';
% flag.medium_approximation='consider_real(sigma)_only';
% flag.medium_approximation='sigma_dc';




flag.sweep=false;
flag.sweep=true; % 参数扫描则为真

flag.data_map=false;
flag_using_stored_data=false;


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
plasma=plasma_model(flag, input.plasma);


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