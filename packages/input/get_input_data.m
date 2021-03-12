function [ input ] = get_input_data( flag )
% 返回预置 模型输入数据

%% input of  equivalent_medium_model_of_plasma
fprintf('[INFO] Use input plasma dataset: %s \n',flag.input_plasma);
% 导入预置数据
switch flag.input_plasma
    case 'CHARLIE_10Pa_4MHz_520W'
        % 一个欧姆加热占主的特殊情况
        input.plasma.flag=flag.input_plasma;
        % 从一个 multi_coupled 数据中取出部分
        temp=get_input_plasma( '2019Raunera_CHARLIE_sweep' );
        idx=(temp.p==10) & (temp.f==4e6) & (temp.Pin==520);
        input.plasma.idx=idx; % 用于get_external
        input.plasma.p=temp.p(idx);
        input.plasma.f=temp.f(idx);
        input.plasma.Pin=temp.Pin(idx);
        input.plasma.ne=temp.ne(idx);
        input.plasma.Te=temp.Te(idx);
        input.plasma.Tg=temp.Tg(idx);
        input.plasma.w_RF=temp.w_RF(idx);
        input.plasma.ng=temp.ng(idx);
        input.plasma.size=[length(input.plasma.p),...
            length(input.plasma.f),length(input.plasma.Pin)];
    otherwise
        input.plasma=get_input_plasma( flag.input_plasma );
end

%% 等离子体不均匀分布处理
input.plasma.r=inf; %均匀无限大/半无限大

%% 外场信息
% 磁场

%% 等离子体特征参数计算
constants=get_constants();
input.plasma.wpe=get_omega_pe(input.plasma.ne); %电子等离子体频率
input.plasma.wpi=get_omega_pi(input.plasma.ne,1,1); %离子等离子体频率
input.plasma.ve=sqrt(8*input.plasma.Te*constants.e/(pi*constants.me));    %电子平均速率，计算自由程用
%             plasma.wce=constants.e;     %电子拉莫尔运动频率

%% input of electric model
if ~isempty(flag.electric_model)
    % 几何数据
    input.geometry=get_input_geometry( flag.input_geometry ); % 导入预置数据
    
    % 假设等离子体几何尺寸，与放电腔同中心的圆柱
    %简化：以放电腔半径为等离子体半径
    input.geometry.r_plasma_eff=input.geometry.r_chamber;
    input.geometry.l_plasma=input.geometry.l_chamber;
    % input.geometry.l_plasma=2*input.geometry.l_coil; %几何校正
    % 当Lmp≥Rp/veff，可见Rp对PER影响小，则lplasma对PER影响小
    input.plasma.r=input.geometry.r_plasma_eff;
    
    % 外电路数据
    switch input.plasma.flag
        case 'CHARLIE_10Pa_4MHz_520W'
            temp=get_input_external(flag, input.geometry, input.plasma);
            input.external=temp;
            % Rcoil_th, Lcoil_th, Lcoil_ex等只有单个数据
            % 从多维数组中选择
            idx=input.plasma.idx;
            input.external.Rmetal_ex=temp.Rmetal_ex(idx);
            input.external.Rcoil_ex=temp.Rcoil_ex(idx);
            input.external.Icoil_rms=temp.Icoil_rms(idx);
            if strcmp(flag.Rmetal,'measured-Rmetal-woplasma')
                input.external.Rmetal=temp.Rmetal(idx);
            end
            input.external.Qcoil=temp.Qcoil(idx);
        otherwise
            input.external=get_input_external( flag, input.geometry, input.plasma );
    end
end

end