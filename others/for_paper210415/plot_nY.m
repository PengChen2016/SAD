function handle_fig=plot_nY(varargin)
% 批量绘图
% (2n+2)个输入，依次为Yi, name_Yi, name_Y, scale_type
%% basic
line_style={'-','-.','--'};
% line_marker={'','o','x','+'};
line_marker={'o','x','+'};
line_color={'r','g','b','c','m','y','k'};
line_d={};
% 全部轮换
for i=line_style
    for j=line_marker
        for k=line_color
            line_d(end+1)=strcat(i,j,k);
        end
    end
end


%% input operation
assert(nargin>4 && rem(nargin,2) == 0 && nargin<=length(line_d))
num_Y=floor(nargin/2)-1;

%% plot
if strcmp(varargin{end},'all')
    handle_fig(1)=plot_nY(varargin{1:end-1},'plot');
    handle_fig(2)=plot_nY(varargin{1:end-1},'loglog');
    handle_fig(3)=plot_nY(varargin{1:end-1},'semilogx');
else
    p=[0.3, 0.5, 1, 3, 5, 10]';
    handle_fig{1}=figure;
    legend_text={};
    switch varargin{end}
        case 'plot'
            plot(p,varargin{1},line_d{1});
        case 'loglog'
            loglog(p,varargin{1},line_d{1});
        case 'semilogx'
            semilogx(p,varargin{1},line_d{1});
    end
    legend_text{end+1}=varargin{2};
    
    hold on
    xticks(p)
    xlabel('{\itp} [Pa]');
    axis([0.3,10,-inf,inf])
    grid on
    
    for i=2:num_Y
        switch varargin{end}
            case 'plot'
                plot(p,varargin{2*i-1},line_d{i});
            case 'loglog'
                loglog(p,varargin{2*i-1},line_d{i});
            case 'semilogx'
                semilogx(p,varargin{2*i-1},line_d{i});
        end
        legend_text{end+1}=varargin{2*i};
    end

    ylabel(varargin{end-1});
    
    L1=legend(legend_text);
    set(L1,'Location','best');
    set(L1,'AutoUpdate','off');
end

end