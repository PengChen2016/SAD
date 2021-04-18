function output_electric_model( flag, source )
% 对electric model的输出与可视化
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
%             name_Y='等离子体等效阻抗\itZ_{\rmplasma}';
%             Y1=PER;
%             name_Y1='\itPER\rm, namely \itR_{\rmplasma}\rm[\Omega]';
%             Y2=Xplasma;
%             name_Y2='\itX_{\rmplasma}\rm[\Omega]';
%             handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
%                 X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%             line([X1(1),X1(end)],[0,0],'linestyle','-.','linewidth',1.5,'color','k');
%             L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
%                 [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}],...
%                 'Z=0');
%             set(L1,'location','southwest');
%             
%             Y2=-Lplasma*1e6;
%             name_Y2='-\itL_{\rmplasma}\rm[\muH]';
%             handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
%                 X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%             pause
%             
%             name_Y='系统等效阻抗\itZ_{\rmsys}即实验可测的端口阻抗';
%             Y1=Rsys;
%             name_Y1='\itR_{\rmsys}\rm[\Omega]';
%             Y2=Xsys;
%             name_Y2='\itX_{\rmsys}\rm[\Omega]';
%             handle_fig=plot_parametric_2Y(name_Y,Y1,name_Y1,Y2,name_Y2,...
%                 X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%             Y2=Lsys*1e6;
%             name_Y2='\itL_{\rmsys}\rm[\muH]';
%             handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
%                 X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%             pause
%             
%             
%             
%             
%             name_Y='系统功率传输效率PTE和耦合因数PCF';
%             Y1=PTE;
%             name_Y1='\itPTE';
%             Y2=PCF;
%             name_Y2='\itPCF';
%             handle_fig=plot_parametric_2Yaxi(name_Y,Y1,name_Y1,Y2,name_Y2,...
%                 X1,idx_X1,X_var,name_X_var,unit_X_var,legend_X2,no_mid_X2,now_str);
%             L1=legend([name_Y1 ' at ' legend_X2{1}],[name_Y1 ' at ' legend_X2{no_mid_X2}],[name_Y1 ' at ' legend_X2{end}],...
%                 [name_Y2 ' at ' legend_X2{1}],[name_Y2 ' at ' legend_X2{no_mid_X2}],[name_Y2 ' at ' legend_X2{end}]);
%             set(L1,'location','northwest');
%             pause
            
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
end
