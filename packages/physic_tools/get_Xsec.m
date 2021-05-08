function [ Xsec ] = get_Xsec( type, flag_plot )
%% get Xsec from source
switch type
    case 'e-H2-1990Tawara'
        % PengZhao get data from 1990Tawara - Cross Sections and Related Data for
        % Electron Collisions with Hydrogen Molecules and Molecular Ions
        [E,a]=textread('e-H2-ELASTIC-1990Tawara.txt','%n %n','headerlines',10);
        Xsec.enp=[E,a];
        [E,a]=textread('e-H2-IONIZATION-1990Tawara.txt','%n %n','headerlines',10);
        Xsec.eniz=[E,a];
    case 'e-H2-Phelps'
        % PengChen2016 get data from Phelps (details are recorded in the data
        % files)
        [E,a]=textread('e-H2-ELASTIC-Phelps.txt','%n %n','headerlines',12);
        Xsec.enp=[E,a];
        [E,a]=textread('e-H2-IONIZATION-Phelps.txt','%n %n','headerlines',10);
        Xsec.eniz=[E,a];
    case 'e-Ar-Biagi'
        % PengChen2016 get data from LXCat-Biagi-v7.1
        [E,a]=textread('e-Ar-ELASTIC-Biagi-v7.1.txt','%n %n','headerlines',10);
        Xsec.enp=[E,a];
        [E,a]=textread('e-Ar-IONIZATION-Biagi-v7.1.txt','%n %n','headerlines',9);
        Xsec.eniz=[E,a];
end

%% plot
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
    title(type)
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

