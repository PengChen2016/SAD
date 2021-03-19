function [ Xsec ] = get_Xsec( flag_plot )
% 读赵鹏截面文件, 得到弹性与电离反应截面。
[E,a]=textread('ELASTIC.txt','%n %n','headerlines',10);
Xsec.enp=[E,a];
[E,a]=textread('IONIZATION.txt','%n %n','headerlines',10);
Xsec.eniz=[E,a];

if flag_plot
    plot_line_width=3;
    gca_line_width=1;
    font_size=15;
    
    handle_fig=figure;
    loglog(Xsec.enp(:,1),Xsec.enp(:,2),'-r','LineWidth',plot_line_width)
    hold on
    loglog(Xsec.eniz(:,1),Xsec.eniz(:,2),'--k','LineWidth',plot_line_width)
    
    ylabel('Cross section (m^2)');
    xlabel('Energy (eV)')
    set(gca,'FontSize',font_size)
    set(gca, 'LineWidth',gca_line_width)
    %             title([name_Y ' \rmat \rm' now_str]);
    grid on%显示网格
%     text(0.4*X1(1),0.5e5,'(b)','FontSize',font_size)
    
    L1=legend('{\it\bf\sigma}_{en}^p','{\it\bf\sigma}_{en}^{iz}');
    set(L1,'FontSize',font_size);
    set(L1,'location','southwest');
    set(L1,'box','off')
    set(L1,'AutoUpdate','off')
    hold off
end
end

