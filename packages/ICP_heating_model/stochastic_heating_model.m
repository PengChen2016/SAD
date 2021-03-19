function [ plasma ] = stoc_model( type, plasma )
% 随机加热模型
% 在plasma_model.m中使用
constants=get_constants();

fprintf('[INFO] Use stoc model: %s\n',type);

%% 计算
% 基于1995Vahedi的随机加热模型
% 表征电子穿过趋肤层耗时与RF周期比值的参数
alpha_st_fun=@(delta_st,Ve,w_RF) 4*delta_st.^2.*w_RF.^2/pi./(Ve.^2);

switch type
    case 'Cazzador-fit'
        % nu_st Cazzador fit
        f_Cazzador=2.1e6; %Cazzador拟合式使用的驱动频率，单位Hz
        w_RF_Cazzador=2*pi*f_Cazzador;
        nu_st_fun_Cazzador=@(x)10.^(-244.1+48.1*log10(x)-3.467*log10(x).^2+...
            0.1113*log10(x).^3-0.001336*log10(x).^4); % 2.1MHz下的nu_st-x关系，x=neTe
        nu_st_fun_new_f=@(ne, Te, w_RF) nu_st_fun_Cazzador(...
            ne.*Te.*w_RF_Cazzador^2./w_RF.^2).*w_RF/w_RF_Cazzador;
        
        % 拟合表达式使用范围
        X_temp=plasma.ne.*plasma.Te*w_RF_Cazzador^2./plasma.w_RF.^2;
        [row,col] = find(X_temp<1e14 | X_temp>1e21);
        num=length(row);
        if num>0
            % warning(['X=' num2str(X_temp) ' out of the range of Cazzador-fit'])
            X_out_range=X_temp(row,col);
            warning('Some X are out of the range of Cazzador-fit (1e14<X<1e21):')
            for i=1:num
                fprintf('[WARN] X(%d, %d)=%e\n',row(i),col(i),X_out_range(i))
            end
        end
        
        plasma.nu_st=nu_st_fun_new_f(plasma.ne, plasma.Te, plasma.w_RF);
        plasma.delta_st=get_plasma_skin_depth('as-medium-simplified',...
            plasma.f,plasma.nu_st,plasma.wpe,plasma.r);
        plasma.alpha_st=alpha_st_fun(plasma.delta_st,plasma.ve, plasma.w_RF);
        
        pause(0.1)
    case 'Vahedi-simplify' %case 'Cazzador-fit'时，欲测试stoc则注释掉，顺序执行两种nu_st计算
        % 基于1995Vahedia model，略有修正
        size_mat=size(plasma.ve);
        plasma.alpha_st=zeros(size_mat);
        plasma.delta_st=zeros(size_mat);
        plasma.nu_st=zeros(size_mat);
        
        delta_st_Vahedi1=(constants.c^2*plasma.ve./(plasma.wpe.^2)/pi./plasma.w_RF).^(1/3);
        alpha_st_Vahedi1=alpha_st_fun(delta_st_Vahedi1,plasma.ve,plasma.w_RF);
        %%%% if alpha_st_Vahedi1<0.03
        idx1=alpha_st_Vahedi1<0.03;
        plasma.alpha_st(idx1)=alpha_st_Vahedi1(idx1);
        plasma.delta_st(idx1)=delta_st_Vahedi1(idx1);
        plasma.nu_st(idx1)=plasma.ve(idx1)./(2*pi*delta_st_Vahedi1(idx1)); %190922PC，根据1995Vahedia
        %%%% else
        idx2=~idx1;
        plasma.delta_st(idx2)=constants.c./plasma.wpe(idx2);
        if isscalar(plasma.w_RF)
            plasma.alpha_st(idx2)=alpha_st_fun(plasma.delta_st(idx2),plasma.ve(idx2),plasma.w_RF);
        else
            plasma.alpha_st(idx2)=alpha_st_fun(plasma.delta_st(idx2),plasma.ve(idx2),plasma.w_RF(idx2));
        end
        %%%%%%%% if alpha_st_Vahedi2<0.03
        idx3=idx2 & (plasma.delta_st<0.03);
        temp=plasma.delta_st(idx3);
        if ~isempty(temp)
            warning('α分段存在问题')
            plasma.nu_st(idx3)=plasma.ve(idx3)./temp/4;
        end
        %%%% elseif alpha_st_Vahedi2<10
        idx4=idx2 & (plasma.delta_st>0.03) & (plasma.delta_st<10);
        plasma.nu_st(idx4)=plasma.ve(idx4)./plasma.delta_st(idx4)/4;
        %%%% else
        idx5=idx2 & (plasma.delta_st>10);
        plasma.nu_st(idx5)=pi/4./(plasma.w_RF(idx5).^2).*(plasma.ve(idx5)./plasma.delta_st(idx5)).^3;
        %%% end
    case 'Cazzador-simplify'
        error('Not implemented.')
        % 基于2018Jain model
        %         试探逻辑与'Vahedi-simplify'类似，待修改
    otherwise
        error('No such case.')
end

end