function [ plasma ] = stochastic_heating_model( flag, plasma )
% 随机加热模型
% 在plasma_model.m中使用
constants=get_constants();

fprintf('[INFO] Use stoc model: %s\n',flag.stoc_model);
% 表征电子穿过趋肤层耗时与RF周期比值的参数
alpha_st_fun=@(delta_st,Ve,w_RF) 4*delta_st.^2.*w_RF.^2/pi./(Ve.^2);

switch flag.stoc_model
    case 'Cazzador-fit'
        %% Cazzador-fit
        % 2014Cazzador - Analytical and numerical models and first
        % operations on the negative ion source NIO1
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
    case 'Vahedi-simplify'
        %% Vahedi-simplify
        % 1995Vahedia - Analytic model of power deposition in inductively
        % coupled plasma sources. 
        % pre-processing
        if 1==plasma.size
            size_mat=[1,1];
        else
            size_mat=plasma.size(plasma.size>1);
        end
        plasma.alpha_st=zeros(size_mat);
        plasma.delta_st=zeros(size_mat);
        plasma.nu_st=zeros(size_mat);
        % For vectorization parallel
        if isscalar(plasma.Te)
            temp_ve=repmat(plasma.ve,size_mat);
        else
            temp_ve=plasma.ve;
        end
        if isscalar(plasma.ne)
            temp_wpe=repmat(plasma.wpe,size_mat);
        else
            temp_wpe=plasma.wpe;
        end
        if isscalar(plasma.w_RF)
            temp_w_RF=repmat(plasma.w_RF,size_mat);
        else
            temp_w_RF=plasma.w_RF;
        end
        
        delta_st_Vahedi1=(constants.c^2*temp_ve./(temp_wpe.^2)/pi./temp_w_RF).^(1/3);
        alpha_st_Vahedi1=alpha_st_fun(delta_st_Vahedi1,temp_ve,temp_w_RF);
        % α(δst for α<=0.03) <= 0.03
        idx1=alpha_st_Vahedi1<=0.03;
        plasma.alpha_st(idx1)=alpha_st_Vahedi1(idx1);
        plasma.delta_st(idx1)=delta_st_Vahedi1(idx1);
        plasma.nu_st(idx1)=temp_ve(idx1)./(2*pi*delta_st_Vahedi1(idx1));
        % α(δst for α<=0.03) > 0.03
        idx2=~idx1;
        plasma.delta_st(idx2)=constants.c./temp_wpe(idx2); %delta_st_Vahedi2
        plasma.alpha_st(idx2)=alpha_st_fun(plasma.delta_st(idx2),...
            temp_ve(idx2),temp_w_RF(idx2)); %alpha_st_Vahedi2
        % α(δst for α>0.03) <= 0.03 : Wrong expression used
        idx3=idx2 & (plasma.alpha_st<=0.03);
        temp=plasma.delta_st(idx3); %delta_st_Vahedi2
        if ~isempty(temp)
            warning('α(δst for α>0.03) <= 0.03 : Wrong expression used.')
            disp('[WARN] Treated as 0.03 < α(δst for α>0.03) <= 10 instead.')
            plasma.nu_st(idx3)=temp_ve(idx3)./temp/4; 
        end
        % 0.03 < α(δst for α>0.03) <= 10
        idx4=idx2 & (plasma.alpha_st>0.03) & (plasma.alpha_st<=10);
        plasma.nu_st(idx4)=temp_ve(idx4)./plasma.delta_st(idx4)/4;
        % α(δst for α>0.03) > 10
        idx5=idx2 & (plasma.alpha_st>10);
        plasma.nu_st(idx5)=pi/4./(temp_w_RF(idx5).^2).*(temp_ve(idx5)./plasma.delta_st(idx5)).^3;
    case 'Cazzador-simplify'
        %% Cazzador-simplify
        % 2014Cazzador - Analytical and numerical models and first
        % operations on the negative ion source NIO1       
        % Also used in 
        % Try and error at the gap.
        error('Not implemented.') % 20210407 
        
    case '2018Jainb-simplify'
        % 2018Jainb - Studies and experimental activities to qualify the
        % behaviour of RF power circuits for Negative Ion Sources of
        % Neutral Beam Injectors for ITER and fusion experiments
        % 2018Jainb modify Vahedi/Cazzador-simplify.
        % 1. use ve=thermal velocity rather than mean thermal velocity
        % 2. use nu_m to calculate delta_st, then alpha_st, then nu_st
        % pre-processing
        if 1==plasma.size
            size_mat=[1,1];
        else
            size_mat=plasma.size(plasma.size>1);
        end
        plasma.alpha_st=zeros(size_mat);
        plasma.delta_st=zeros(size_mat);
        plasma.nu_st=zeros(size_mat);
        % For vectorization parallel
        if isscalar(plasma.Te)
            temp_ve=repmat(plasma.ve,size_mat);
        else
            temp_ve=plasma.ve;
        end
        if isscalar(plasma.ne)
            temp_wpe=repmat(plasma.wpe,size_mat);
        else
            temp_wpe=plasma.wpe;
        end
        if isscalar(plasma.w_RF)
            temp_f=repmat(plasma.f,size_mat);
            temp_w_RF=repmat(plasma.w_RF,size_mat);
        else
            temp_f=plasma.f;
            temp_w_RF=plasma.w_RF;
        end
        
        if ~isfield(plasma,'nu_m') || isempty(plasma.nu_m)
            warning('No nu_m value. Run Ohmic heating model first.')
            plasma=Ohmic_heating_model(plasma, flag.type_Xsec);
        end
        
        plasma.delta_st=get_plasma_skin_depth('as-medium-simplified-finite-radius',...
            temp_f,plasma.nu_m,temp_wpe,plasma.r);
        plasma.alpha_st=alpha_st_fun(plasma.delta_st,temp_ve,temp_w_RF);
        % case 4: Vahedi/Cazzador-simplify 2 phase. Used.
        x=plasma.ne.*plasma.Te;
        % α <= 1
        idx1=plasma.alpha_st<=1;
        plasma.nu_st(idx1)=temp_ve(idx1)./(2*pi*plasma.delta_st(idx1));
        % α > 1
        idx2=plasma.alpha_st>1;
        % 2018Jain Maybe equal Vahedi phase 3.
        plasma.nu_st(idx2)=pi/4./(temp_w_RF(idx2).^2).*...
            (8*constants.mu0*constants.e^3*x(idx2)/pi/constants.me^2).^(3/2);
    otherwise
        error('No such case.')
end

end