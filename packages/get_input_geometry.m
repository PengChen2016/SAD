function [ input ] = get_input_geometry( flag_input_geometry )
% 返回预置 ICP源几何数据 结构体

% 放电腔为圆柱，线圈为螺线管
% input.r_chamber %放电腔内半径[m]
% input.l_chamber %放电腔长度[m]
% input.r_coil %线圈半径[m]
% input.l_coil %线圈长度m
% input.N_coil %线圈匝数
% input.r_wire %线圈导线半径[m]

fprintf('Use input ICPs geometry dataset: %s \n',flag_input_geometry );
input.flag=flag_input_geometry ;

%% 导入几何参数
switch flag_input_geometry
    case 'ELISE_base'
        % Jain
        input.r_chamber=0.14;
        input.l_chamber=0.14;
        input.r_coil=0.158;
        input.l_coil=0.08; % 一次线圈为短螺线管
        input.N_coil=6;
        input.r_wire=0.004; %200625改正，之前为0.003
    case 'HUST_large_driver_base'
        input.r_chamber=0.142;
        input.l_chamber=0.147;
        input.r_coil=0.155;
        input.l_coil=0.08935;
        input.N_coil=9;
        input.r_wire=0.003;
    case 'BATMAN_base'
        warning('no data')
        % TODO
    case 'HUST_small_driver_base'
        % ZP % 不确定
        input.r_chamber=0.051;
        input.l_chamber=0.137;
        input.r_coil=0.063;
        input.l_coil=0.06; 
        input.N_coil=6;
        input.r_wire=0.003;
    case 'HUST_small_driver_ZL'
        % ZL
        input.r_chamber=0.051;
        input.l_chamber=0.137;
        input.r_coil=0.063;
        input.N_coil=7;
        input.l_coil=12*input.N_coil;
        input.r_wire=0.003;
    case 'CHARLIE_base'
        % Rauner
        input.r_chamber=0.05-0.005;
        input.l_chamber=0.4;
        input.r_coil=0.055;
        input.l_coil=0.1;
        input.r_wire=0.003;
        input.N_coil=5;
    case 'small_source1_LZS'
        % 解析电磁模型：无限长直
        input.r_chamber=0.06;
        input.l_chamber=0.16;
        input.r_coil=input.r_chamber;
        input.l_coil=input.l_chamber;
        input.r_wire=0.003;
        input.N_coil=5;
    case 'NIO1_base'
        warning('no data')
        %TODO: 待修改。目前是copy ELISE数据
end

end