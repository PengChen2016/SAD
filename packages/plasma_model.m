function [ plasma ] = plasma_model( flag, plasma)
% plasma model for electric model of ICP source

%% ICP heating model
plasma=ICP_heating_model( flag, plasma);

%% equivalent EM medium model of plasma
plasma=equivalent_EM_medium_model( flag, plasma);

if isinf(plasma.r)
    disp('[INFO] Radius of plasma is assumed to be infinite.')
else
    idx=find(plasma.wavelength<3*plasma.r);
    if ~isempty(idx)
        disp('[WARN] λ＜≈R的元素的索引')
        disp(idx')
        if isfield(flag,'electric_model') && ~isempty(strfind(flag.electric_model,'transformer'))
            warning('集中参数电路模型有bug')
            %         pause flag.electric_model='transformer-base';
        end
    end
    
    idx=find(plasma.r<3*plasma.skin_depth);
    if ~isempty(idx)
        disp('[WARN] δ>≈R的元素的索引')
        disp(idx')
        if isfield(flag,'stoc_model') && ~isempty(flag.stoc_model)
            warning('随机加热模型计算 有bug')
        end
        if isfield(flag,'electric_model') && ~isempty(strfind(flag.electric_model,'transformer'))
            warning('变压器模型中，计算Rp-Lp的圆柱媒质涡流模型适用性存疑.')
        end
    end
end
%% output
if flag.output_plasma_model
    output_plasma_model(flag,plasma)
end

end

