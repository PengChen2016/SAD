function [ k ] = get_Nagaoka( a)
% 对无限长螺线管电感公式做几何修正，以计算有限长螺线管的电感
% a：螺线管线圈2R/l 直径/长度，0<a<10
% k：长冈系数

%% 导入数据
% 长冈系数的求解：见https://ja.wikipedia.org/wiki/長岡係数
% 长冈系数的数据：见1909Nagaoka - The inductance coefficients of solenoids
Nagaoka_data=textread('Nagaoka.txt');

if ~isnan(a)
    %% 插值
    k = interp1(Nagaoka_data(:,1),Nagaoka_data(:,2),a,'linear',0);
    idx=find(a<Nagaoka_data(1,1) | a<Nagaoka_data(1,end),1);
    if ~isempty(idx)
        warning('[WARN] Out of the range of Nagaoka.txt, make k=0')
        disp(idx')
    end
else
    k=nan;
    %% 显示
    handle_fig=figure;
    plot(Nagaoka_data(:,1),Nagaoka_data(:,2))
    ylabel('k');
    xlabel('a=2R/l')
    title('Nagaoka coefficient');
    grid on%显示网格
end
end

