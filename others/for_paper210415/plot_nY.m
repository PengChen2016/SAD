function plot_nY(varargin)
% 批量绘图
% (2n+2)个输入，依次为Yi, name_Yi, name_Y, scale_type


assert(nargin>4 && rem(nargin,2) == 0)
num_Y=floor(nargin/2)-1;


scale_type={'loglog','semilogx','plot'};




p=[0.3, 0.5, 1, 3, 5, 10]';

figure;
legend_text={};

loglog(p,experiment.PER(:,1),'-xk');
legend_text{end+1}='subtractive method, 1MHz';

hold on
xticks(p)
xlabel('{\itp} [Pa]');
axis([0.3,10,-inf,inf])
grid on

loglog(p,fem.dielectric_PER(:,1),'-or','MarkerSize',8);
legend_text{end+1}='FEM model, 1MHz';
loglog(p,source_raza.PER(:,1),'-sb');
legend_text{end+1}='transformer model, 1MHz';
ylabel('PER');
yticks([0.1,0.2,0.5,1,2,5,10])
axis([0.3,10,0.1,10])

L1=legend(legend_text);
set(L1,'Location','best');
set(L1,'AutoUpdate','off');

end