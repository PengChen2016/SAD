function [ source ] = electric_model( flag, input )
% electric model of ICP source (cylinder coil)
% 如果将input_geo和input_external赋给source，则改为
% function [ source ] = electric_model( flag, plasma, source )
    
%% preparation
fprintf('[INFO] Use electric model: %s\n',flag.electric_model);

source.w_RF=input.plasma.w_RF;
% 若将Icoil_rms移出具体model，则以下条件语句移到本文件最后
if ~isfield(input.external,'Icoil_rms') || isempty(input.external.Icoil_rms)...
        || isempty(find(~isnan(input.external.Icoil_rms),1))
    input.external.Icoil_rms=1;
    disp('[WARN] No given Icoil. Use Icoil_rms=1A to calculate power.')
end

%% choose electric model
if strfind(flag.electric_model,'transformer')
    source=transformer_model( flag, input );
else
    switch flag.electric_model
        case 'analytical_base'
            source=analytical_EM_model( flag, input );
        case 'multi-filament'
            error('To be realized.')
        otherwise
            error('Unexpected electric model. Please stop and check.')
    end
end

%% derived parameters
source.size=input.plasma.size;

% 如果 Pin或Icoil_rms相关的运算 与具体电模型无关，则提出来放在这里
% source.Pin=input.plasma.Pin; % Pin需保留在input_plasma内，是常用表征参数
% 与Pin或Icoil_rms相关的运算
% 一般至少有Pin
% 有Pin，则由Rsys计算Isys
% 有Icoil_rms，则由Rsys计算Psys
% 同时有Pin和Icoil_rms，则计算Rsys_given
% source.Rsys_experiment=input.plasma.Pin./input.external.Icoil_rms

source.PTE=source.PER./source.Rsys;   %射频功率传输效率 PTE
% source.PCF=source.Rsys./source.Xsys;  %激励器射频功率耦合因数

%% output
if isfield(flag,'output_electric_model') && flag.output_electric_model
    output_electric_model( flag, source )
end

end