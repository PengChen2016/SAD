function [ plasma ] = get_input_plasma( flag_input_plasma )
% 返回预置 等离子体输入数据 结构体
% Please use get_input_data with flag.electric_model='' rather than use
% this function directly.

% input.f %驱动频率，单位Hz
% input.ne %电子密度[m^-3]
% input.Te  %电子温度[eV]
% input.p %气压[Pa]
% input.Tg %气体温度[K]
% input.w_RF %驱动角频率，单位Hz
% input.ng %中性气体分子密度[m^-3]
% input.Pin % 总输入功率 W

plasma.flag=flag_input_plasma;
constants=get_constants();% 全局常数 结构体

plasma.size=1; %用于指示size
switch flag_input_plasma
    %% 单点计算
    case '2018Jainb_ELISE_typical'
        % Used in get_example_flag(0). Please do not modified.
        plasma.f=1e6;
        plasma.ne=5e18; %center
        plasma.Te=9; % 2018Jainb
        plasma.p=0.3;
        plasma.Tg=1200;
        plasma.Pin=nan;
    case 'BATMAN_typical'
        plasma.f=1e6;
        %2010Mcneely中为中心1.5e18，>10eV
        plasma.ne=1.5e18;  %中心值
        plasma.Te=10;
        %2006Fantz报道激励器内1000K，2018Fantz报道整体630K
        plasma.p=0.3;
        plasma.Tg=630;
        plasma.Pin=nan;
    case '2021Zielke_BATMAN_sweep'
        % TODO
        % 尚未处理
        error('no data')
    case 'small_source1_LZS'
        % 用于与LZS的理论电磁模型做对比
        plasma.f=1e6;
        plasma.ne=1.6e17;
        plasma.Te=6.5;
        plasma.p=0.6;
        plasma.Tg=600;
        plasma.Pin=nan;
    case '2018Jainb_NIO1_sweep'
        % 尚未处理
        error('no data')
        plasma.f=2.1e6;
        plasma.p=1;
        plasma.Tg=400;
        plasma.Pin=nan;
        %% 多点计算
        % length(input.size)>1，存在多值变量
        % 根据非耦合的多值变量，扩展为多维数组
        % 使用repmat扩展，请注意先reshape
        
        %%%%%%%%% 情况一：所有多值变量非耦合
        % length(input.size)>3，各多值变量扩展为多维数组
        % 扩展后断言检验
        % 一般用于扫描参数空间
        % 主要使用length(size)信息
        % 单值变量不做扩展
    case '2020Chen_NIS_sweep1'
        % Used in get_example_flag(2). Please do not modified.
        plasma.p=0.3;
        plasma.f=1e6;
        plasma.Pin=nan;
        plasma.ne=logspace(16,19,40);
        plasma.Te=5:5:25;
        plasma.Tg=630;
        
        plasma.size=[length(plasma.p),length(plasma.f),length(plasma.Pin),...
            length(plasma.ne),length(plasma.Te)];
        % 扩展
        % 各变量，依次占据一个维度
        plasma.ne=reshape(plasma.ne,[],1);
        plasma.Te=reshape(plasma.Te,1,[]);
        % 各变量，在自身维度上为1，在其他维度上按相应变量size扩展
        plasma.ne=repmat(plasma.ne, 1,plasma.size(5));
        plasma.Te=repmat(plasma.Te, plasma.size(4),1);
        % 扩展后断言检验
        assert(isequal(plasma.ne(:,1),logspace(16,19,40)'))
        assert(isequal(plasma.Te(1,:),5:5:25))
        
    case '2020Chen_NIS_sweep_p'
        plasma.p=0.3:0.3:10; % 2019Raunera使用0.3~10Pa
        plasma.f=1e6;
        plasma.Pin=nan;
        plasma.ne=logspace(16,19,31);
        plasma.Te=1:4:21;
        plasma.Tg=630;
        
        plasma.size=[length(plasma.p),length(plasma.f),length(plasma.Pin),...
            length(plasma.ne),length(plasma.Te)];
%         size(size>1);        
        % 扩展
        % 各变量，依次占据一个维度
        plasma.p=reshape(plasma.p,[],1,1);
        plasma.ne=reshape(plasma.ne,1,[],1);
        plasma.Te=reshape(plasma.Te,1,1,[]);
        % 各变量，在自身维度上为1，在其他维度上按相应变量size扩展
        plasma.p=repmat(plasma.p, 1,plasma.size(4),plasma.size(5));
        plasma.ne=repmat(plasma.ne, plasma.size(1),1,plasma.size(5));
        plasma.Te=repmat(plasma.Te, plasma.size(1),plasma.size(4),1);
        
        % 扩展后断言检验
        assert(isequal(plasma.p(:,1,1),(0.3:0.3:10)'))
        assert(isequal(plasma.ne(1,:,1),logspace(16,19,31)))
        assert(isequal(plasma.Te(1,1,:),reshape(1:4:21,1,1,[])))
    case '2018Jainb_ELISE_sweep_f'
        % fig40 of 2018Jainb - Studies and experimental activities to
        % qualify the behaviour of RF power circuits for Negative Ion
        % Sources of Neutral Beam Injectors for ITER and fusion experiments
         plasma.p=0.3;
        plasma.f=logspace(4,9,141);
        plasma.Pin=nan;
        dne=16:1:19;
        plasma.ne=5*10.^dne;
        plasma.Te=9;
        plasma.Tg=1200;
        
        plasma.size=[length(plasma.p),length(plasma.f),length(plasma.Pin),...
            length(plasma.ne),length(plasma.Te)];
        % 扩展
        % 各变量，依次占据一个维度
        plasma.f=reshape(plasma.f,[],1);
        plasma.ne=reshape(plasma.ne,1,[]);
        % 各变量，在自身维度上为1，在其他维度上按相应变量size扩展
        plasma.f=repmat(plasma.f, 1,plasma.size(4));
        plasma.ne=repmat(plasma.ne, plasma.size(2),1);
        
        % 扩展后断言检验
        assert(isequal(plasma.f(:,1),logspace(4,9,141)'))
        assert(isequal(plasma.ne(1,:),5*logspace(16,19,4)))
    case '2018Jainb_ELISE_sweep_all'
        % fig40/41 of 2018Jainb - Studies and experimental activities to
        % qualify the behaviour of RF power circuits for Negative Ion
        % Sources of Neutral Beam Injectors for ITER and fusion experiments
        plasma.p=0.05:0.05:1;
        plasma.f=logspace(4,9,141);
        plasma.Pin=nan;
        plasma.ne=5*logspace(16,19,4);
        plasma.Te=[5,9,10,15,20,25];
        plasma.Tg=1200;
        
        plasma.size=[length(plasma.p),length(plasma.f),length(plasma.Pin),...
            length(plasma.ne),length(plasma.Te)];
        % 扩展
        % 各变量，依次占据一个维度
        plasma.p=reshape(plasma.p,[],1,1,1);
        plasma.f=reshape(plasma.f,1,[],1,1);
        plasma.ne=reshape(plasma.ne,1,1,[],1);
        plasma.Te=reshape(plasma.Te,1,1,1,[]);
        % 各变量，在自身维度上为1，在其他维度上按相应变量size扩展
        plasma.p=repmat(plasma.p, 1,plasma.size(2),plasma.size(4),plasma.size(5));
        plasma.f=repmat(plasma.f, plasma.size(1),1,plasma.size(4),plasma.size(5));
        plasma.ne=repmat(plasma.ne, plasma.size(1),plasma.size(2),1,plasma.size(5));
        plasma.Te=repmat(plasma.Te, plasma.size(1),plasma.size(2),plasma.size(4),1);
        
        % 扩展后断言检验
        assert(isequal(plasma.p(:,1,1,1),(0.05:0.05:1)'))
        assert(isequal(plasma.f(1,:,1,1),logspace(4,9,141)))
        assert(isequal(plasma.ne(1,1,:,1),reshape(5*logspace(16,19,4),1,1,[],1)))
        assert(isequal(plasma.Te(1,1,1,:),reshape([5,9,10,15,20,25],1,1,1,[])))
        
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
        plasma.p=repmat(p,1,num_f,num_Pin);
        plasma.f=repmat(f,num_p,1,num_Pin);
        plasma.Pin=repmat(Pin,num_p,num_f,1);
        plasma.ne=[...
            7.0E+16, 1.2E+17;...
            1.3E+17, 1.6E+17;...
            1.7E+17, 2.3E+17;...
            2.4E+17, 2.7E+17;...
            2.4E+17, 3.3E+17;...
            2.2E+17, 3.7E+17]; % 轴向平均径向中心
        Te=[5.5, 4.79, 4.06, 3, 2.62, 2.1]'; % 轴向平均径向中心
        plasma.Te=[Te, Te];
        Tg_CHARLIE=@(p) 128.7*log10(p)+568.5; % 根据2017Rauner_CHARLIE数据拟合
        plasma.Tg=Tg_CHARLIE(plasma.p);
        % Tg=[500.364, 529.95, 574.556, 623.714, 652.39, 704.051]';
        
        plasma.size=[num_p, num_f, num_Pin];
    otherwise
        error('no data')
end

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

%% 派生参数
plasma.w_RF=2*pi*plasma.f;
% 理想气体状态方程，ng用于碰撞频率计算
plasma.ng=plasma.p./(constants.kB*plasma.Tg);

%% 等离子体不均匀分布处理
plasma.r=inf; %均匀无限大/半无限大

%% 外场信息
% 磁场

end


