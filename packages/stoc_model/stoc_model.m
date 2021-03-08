function [ plasma ] = stoc_model( type, plasma )
% 随机加热模型
% 在plasma_model.m中使用
constants=get_constants();

fprintf('Use stoc model: %s\n',type);

%% 计算
% 基于1995Vahedi的随机加热模型
% 表征随机加热显著程度的因子
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
        X_temp=plasma.ne.*plasma.Te*w_RF_Cazzador^2./w_RF.^2;
        X_out_range=Xtemp(X_temp<1e14 | X_temp>1e21);
        num=length(X_out_range(:));
        if num>0
            % warning(['X=' num2str(X_temp) ' out of the range of Cazzador-fit'])
            warning('Some X are out of the range of Cazzador-fit:')
            disp(X_out_range)
            pause
        end
        
        plasma.nu_st=nu_st_fun_new_f(plasma.ne, plasma.Te, plasma.w_RF);
        plasma.delta_st=get_plasma_skin_depth('as-medium-simplified',...
    plasma.f,plasma.nu_st,plasma.wpe,plasma.r);
        plasma.alpha_st=alpha_st_fun(plasma.delta_st,plasma.Ve, plasma.w_RF);
    case 'Vahedi-simplify' %case 'Cazzador-fit'时，欲测试stoc则注释掉，顺序执行两种nu_st计算
        % 基于1995Vahedia model，略有修正
        delta_st_Vahedi1=(constants.c^2*plasma.ve/(plasma.wpe^2)/pi/w_RF)^(1/3); %趋肤深度
        delta_st_Vahedi2=constants.c/plasma.wpe;
        alpha_st_Vahedi1=4*delta_st_Vahedi1^2*w_RF^2/pi/(plasma.ve^2); %表征电子穿过趋肤层耗时与RF周期比值的参数
        alpha_st_Vahedi2=4*delta_st_Vahedi2^2*w_RF^2/pi/(plasma.ve^2);
        if alpha_st_Vahedi1<0.03
            plasma.nu_st=plasma.ve/(2*pi*delta_st_Vahedi1); %190922PC，根据1995Vahedia
            % 测试stoc model
            plasma.alpha_st=alpha_st_Vahedi1;
            plasma.delta_st=delta_st_Vahedi1;
            %             vst2(X1i,X2i)=(4/pi)^0.2*(w_RF^0.4)*(plasma.ve/delta_st_Vahedi1)^0.6; %待询问2019Zuo
            %             vst3(X1i,X2i)=-pi*plasma.ve/4/delta_st_Vahedi1/(log(alpha_st_Vahedi1)+1.58); %待询问2019Zuo
            %             fprintf('%s = %.2e \n','nu_st',plasma.nu_st);
            %             fprintf('%s = %.2e \n','vst2',vst2(X1i,X2i));
            %             fprintf('%s = %.2e \n','vst3',vst3(X1i,X2i));
        else
            if alpha_st_Vahedi2<0.03
                warning('α分段存在问题')
                %                 pause %测试时使用
                % 200422 ne=1e16,Te=15/20两次,出现bug
                % 仍然使用α>0.03时的δst算出来的α
                plasma.nu_st=plasma.ve/delta_st_Vahedi2/4;
            else
                if alpha_st_Vahedi2<10
                    plasma.nu_st=plasma.ve/delta_st_Vahedi2/4;
                else
                    plasma.nu_st=pi/4/(w_RF^2)*(plasma.ve/delta_st_Vahedi2)^3;
                end
            end
            plasma.alpha_st=alpha_st_Vahedi2;
            plasma.delta_st=delta_st_Vahedi2;
        end
    case 'Cazzador-simplify'
        % 基于2018Jain model，略有修正。copy LZS代码
        %         试探逻辑与'Vahedi-simplify'类似，待修改
        %
        %         delta_st_Vahedi1=(constants.c^2*plasma.ve/(plasma.wpe^2)/pi/w_RF)^(1/3); %趋肤深度
        %         delta_st_Vahedi2=constants.c/plasma.wpe;
        %         alpha_st_Vahedi1=4*delta_st_Vahedi1^2*w_RF^2/pi/(plasma.ve^2); %表征电子穿过趋肤层耗时与RF周期比值的参数
        %         alpha_st_Vahedi2=4*delta_st_Vahedi2^2*w_RF^2/pi/(plasma.ve^2);
        %         if alpha_st_Vahedi1<0.03
        %             plasma.nu_st=plasma.ve/(2*pi*delta_st_Vahedi1); %190922PC，根据1995Vahedia
        %             plasma.delta_st=delta_st_Vahedi1;
        %             %             vst2(X1i,X2i)=(4/pi)^0.2*(w_RF^0.4)*(plasma.ve/delta_st_Vahedi1)^0.6; %待询问2019Zuo
        %             %             vst3(X1i,X2i)=-pi*plasma.ve/4/delta_st_Vahedi1/(log(alpha_st_Vahedi1)+1.58); %待询问2019Zuo
        %             %             fprintf('%s = %.2e \n','nu_st',plasma.nu_st);
        %             %             fprintf('%s = %.2e \n','vst2',vst2(X1i,X2i));
        %             %             fprintf('%s = %.2e \n','vst3',vst3(X1i,X2i));
        %         else
        %             if alpha_st_Vahedi2<0.03
        %                 warning('α分段存在问题')
        %                 %                 pause %测试时使用
        %                 % 200422 ne=1e16,Te=15/20两次,出现bug
        %                 % 仍然使用α>0.03时的δst算出来的α
        %                 plasma.nu_st=plasma.ve/delta_st_Vahedi2/4;
        %             else
        %                 if alpha_st_Vahedi2<10
        %                     plasma.nu_st=plasma.ve/delta_st_Vahedi2/4;
        %                 else
        %                     plasma.nu_st=pi/4/(w_RF^2)*(plasma.ve/delta_st_Vahedi2)^3;
        %                 end
        %                 plasma.delta_st=delta_st_Vahedi2;
        %             end
        %         end
end
%case {'Vahedi-simplify','Cazzador-simplify'}时测试stoc模型用

switch flag_vst_expression
    case {'Vahedi-simplify','Cazzador-simplify'}
        nu_st_fit(X1i,X2i)=nu_st_fun_new_f(plasma.ne*plasma.Te);
        delta_st_fit(X1i,X2i)=delta_simplified_fun(nu_st_fit(X1i,X2i),plasma.wpe);
        alpha_st_fit(X1i,X2i)=alpha_st_fun(delta_st_fit(X1i,X2i),plasma.ve);
end


end