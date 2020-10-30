function [ k ] = k_eniz( Te )
%输入电子温度Te(eV),得到其电子中性粒子电离碰撞反应系数。
%   此处显示详细说明
e=1.6e-19;           %电子电量
me=9.1e-31;           %电子质量
mi=1.67e-27;          %氢离子质量
epsilon0=8.85e-12;    %真空介电常数
c=3e8;                %真空光速[m/s]
kB=1.381e-23;         %玻尔兹曼常数
u0=4*pi*1e-7;         %真空磁导率

%p=1;
%Tg=400;

%ng=p/(kB*Tg);  

[E,a]=textread('IONIZATION.txt','%n %n','headerlines',10); %读截面文件
[m,n]=size([E,a]);
k=trapz(E,a.*E.*exp(-E/Te))/Te^1.5*sqrt(8*e/pi/me); %截面对Maxwellian分布积分

% TODO：参考Main.m中代码，做截面-能量关系可视化
% plot_line_width=3;
% gca_line_width=1;
% marker_size=8;
% font_size=15;
% marker_indices=1:5:length(X1);
% 
% handle_fig=figure;
% loglog(E,a,'-r','LineWidth',plot_line_width)
% hold on
% % line([X1(1),X1(end)],[1,1],'linestyle',':','linewidth',1*plot_line_width,'color','k');
% 
% ylabel('Cross section (m2)');
% xlabel('Energy (eV)')
% set(gca,'FontSize',font_size)
% set(gca, 'LineWidth',gca_line_width)
% %             title([name_Y ' \rmat \rm' now_str]);
% grid on%显示网格
% text(0.4*X1(1),0.5e5,'(b)','FontSize',font_size)
% 
% L1=legend('{\it\bf\sigma}_{en}^p');
% set(L1,'FontSize',font_size);
% set(L1,'location','northwest');
% set(L1,'box','off')
% set(L1,'AutoUpdate','off')
        
        
end

