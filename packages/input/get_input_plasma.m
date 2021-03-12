function [ input ] = get_input_plasma( flag_input_plasma )
% 返回预置 等离子体输入数据 结构体
% sweep数据都是coupled的，参数一一对应

% input.f %驱动频率，单位Hz
% input.ne %电子密度[m^-3]
% input.Te  %电子温度[eV]
% input.p %气压[Pa]
% input.Tg %气体温度[K]
% input.w_RF %驱动角频率，单位Hz
% input.ng %中性气体分子密度[m^-3]
% input.Pin % 总输入功率 W

input.flag=flag_input_plasma;

constants=get_constants();% 全局常数 结构体
Tg_CHARLIE=@(p) 128.7*log10(p)+568.5; % 根据2017Rauner_CHARLIE数据拟合

input.size=1; %用于指示size
switch flag_input_plasma
    %% 单点计算
    case '2018Jainb_ELISE_typical'
        % 假设值
        input.f=1e6;
        input.ne=5e18; %中心值
        input.Te=9; % 2018Jainb
        input.p=0.3;
        input.Tg=1200;
        %2018Jain的ELISE case使用5e18，15eV，1200-1500K
        % input.Te=15;  % 2018Jain
        input.Pin=nan;
    case 'BATMAN_typical'
        input.f=1e6;
        %2010Mcneely中为中心1.5e18，>10eV
        input.ne=1.5e18;  %中心值
        input.Te=10;
        %2006Fantz报道激励器内1000K，2018Fantz报道整体630K
        input.p=0.3;
        input.Tg=630;
        input.Pin=nan;
    case '2021Zielke_BATMAN_sweep'
        % TODO
        % 尚未处理
        error('no data')
    case 'small_source1_LZS'
        % 用于与LZS的理论电磁模型做对比
        input.f=1e6;
        input.ne=1.6e17;
        input.Te=6.5;
        input.p=0.6;
        input.Tg=600;
        input.Pin=nan;
    case '2018Jainb_NIO1_sweep'
        % 尚未处理
        error('no data')
        input.f=2.1e6;
        input.p=1;
        input.Tg=400;
        input.Pin=nan;
        %% 多点计算
        % length(input.size)>1，存在多值变量
        % 根据非耦合的多值变量，扩展为多维数组
        % 使用repmat扩展，请注意先reshape
        
        %%%%%%%%% 情况一：所有多值变量非耦合
        % length(input.size)>3，各多值变量扩展为多维数组
        % 一般用于扫描参数空间
        % 主要使用length(size)信息
        % 单值变量不做扩展
    case '2020Chen_NIS_sweep1'
        input.p=0.3;
        input.f=1e6;
        input.Pin=nan;
        dne=(16:0.1:19)';  % 指数
        input.ne=10.^dne;
        input.Te=5:5:25;
        input.Tg=630;
        
        input.size=[length(input.p),length(input.f),length(input.Pin),...
            length(input.ne),length(input.Te)];
        % 扩展
        input.ne=repmat(input.ne, 1,input.size(5));
        input.Te=repmat(input.Te, input.size(4),1);
    case '2020Chen_NIS_sweep_p'
        input.p=0.3:0.3:10; % 2019Raunera使用0.3~10Pa
        % input.p=0.1:0.1:1; % 2018Jainb使用0.1~1Pa
        input.f=1e6;
        input.Pin=nan;
        dne=16:0.1:19;
        input.ne=10.^dne;
        input.Te=1:4:21;
        input.Tg=630;
        
        input.size=[length(input.p),length(input.f),length(input.Pin),...
            length(input.ne),length(input.Te)];
%         size(size>1);        
        % 扩展
        input.p=reshape(input.p,[],1,1);
        input.ne=reshape(input.ne,1,[],1);
        input.Te=reshape(input.Te,1,1,[]);
        input.p=repmat(input.p, 1,input.size(4),input.size(5));
        input.ne=repmat(input.ne, input.size(1),1,input.size(5));
        input.Te=repmat(input.Te, input.size(1),input.size(4),1);
        
        %%%%%%%%% 情况二：存在耦合的多值变量
        % length(input.size)==3，全部变量扩展为三维数组
        % 等离子体参数与工况参数三元组一一对应
        % 工况参数p、f、Pin依次为第1/2/3维
        % 其他变量按与工况参数映射关系，存在三维数组中
    case '2019Raunera_CHARLIE_sweep'
        % 原始数据
        p=[0.3, 0.5, 1, 3, 5, 10]'; % 维度一
        f=[1e6, 4e6]; % 维度二
        Pin=520; % 维度三
        num_p=length(p);
        num_f=length(f);
        num_Pin=length(Pin);
        % 相同size的多维数组作为结构体元素
        input.p=repmat(p,1,num_f,num_Pin);
        input.f=repmat(f,num_p,1,num_Pin);
        input.Pin=repmat(Pin,num_p,num_f,1);
        input.ne=[...
            7.0E+16, 1.2E+17;...
            1.3E+17, 1.6E+17;...
            1.7E+17, 2.3E+17;...
            2.4E+17, 2.7E+17;...
            2.4E+17, 3.3E+17;...
            2.2E+17, 3.7E+17]; % 轴向平均径向中心
        Te=[5.5, 4.79, 4.06, 3, 2.62, 2.1]'; % 轴向平均径向中心
        input.Te=[Te, Te];
        input.Tg=Tg_CHARLIE(input.p);
        % Tg=[500.364, 529.95, 574.556, 623.714, 652.39, 704.051]';
        
        input.size=[num_p, num_f, num_Pin];
        %         % 结构体数组：冗余存储，但便于作为对象被统一引用
        %         % 虽然组建较麻烦，不适用于矩阵操作，但概念简单
        %         input=struct('p',0,'f',0,'Pin',0,'ne',0,'Te',0,'Tg',0,'w_RF',0,'ng',0);
        %         input=repmat(input,num_p,num_f,num_Pin);
        %         for i_p=1:num_p
        %             for i_f=1:num_f
        %                 for i_Pin=1:num_Pin
        %                     input(i_p,i_f,i_Pin).p=p(i_p);
        %                     input(i_p,i_f,i_Pin).f=f(i_f);
        %                     input(i_p,i_f,i_Pin).Pin=Pin(i_Pin);
        %                     input(i_p,i_f,i_Pin).ne=ne(i_p,i_f);
        %                     input(i_p,i_f,i_Pin).Te=Te(i_p,i_f);
        %                     input(i_p,i_f,i_Pin).Tg=Tg(i_p);
        %                 end
        %             end
        %         end
    otherwise
        error('no data')
end

input.w_RF=2*pi*input.f;
% 理想气体状态方程，ng用于碰撞频率计算
input.ng=input.p./(constants.kB*input.Tg);
end


