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
flag.using_stored_data=false;
% flag.using_stored_data=true;
%%%%%%%% plasma model
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

% stoc表达式
% flag.stoc_model='';
% flag.stoc_model='Vahedi-simplify';
% flag.stoc_model='Cazzador-simplify';
flag.stoc_model='Cazzador-fit';

flag.medium_approximation='';
% flag.medium_approximation='consider_real(sigma)_only';
% flag.medium_approximation='sigma_dc';

% flag.output_plasma_model=true;
flag.output_plasma_model=false;

%%%%%%%% electric model
% flag.electric_model='';
flag.electric_model='transformer-base';
% flag.electric_model='transformer-2011Chabert';
% flag.electric_model='transformer-2015Bandyopadhyay';
% flag.electric_model='transformer-2018Jainb';
% flag.electric_model='analytical_base';

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

% flag.output_electric_model=true;
flag.output_electric_model=false;

if ~flag.using_stored_data
    %% 计算
    input=get_input_data( flag );
    input.plasma=plasma_model(flag, input.plasma);
    
    if ~isempty(flag.electric_model)
        source=electric_model( flag, input );
    end
else
    %% 载入已计算结果
    % 参考cubeproject
    file_name='stored_data_test200529_2.mat';
    % 已存储则注释掉下列语句
    if flag.using_stored_data
        warning('Store data first.')
        error('end')
    end
    % 在flag.using_stored_data=false时运行以下三行一次,然后将flag.using_stored_data改回true
    % if flag.sweep
    %     flag.using_stored_data=true;
    %     save(file_name)
    %     error('end')
    % end
    
%     if flag.using_stored_data
%         load(file_name)
%         warning('using stored data')
%     end
end

%% 输出与可视化
% output_plasma_model(flag, input.plasma)
% output_electric_model( flag, source )
