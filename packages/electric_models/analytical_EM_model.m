function [ source ] = analytical_EM_model( flag, input )
% 解析模型 2011Chabert
 disp('## ICP源电模型：轴对称平行平面解析模型')

        
        % 介质窗相对介电常数，玻璃约2，氧化铝陶瓷约9
        epst_r=9;
        k_t_wave=w_RF*sqrt(constants.mu0*constants.eps0*epst_r);
        if k_t_wave*r_p>1
            warn('介质管中磁场不是常数，该模型不适用')
        end
        % 重构：将以上参数代码提到最前面输入部分去
        Hz0=N_equ*Icoil_rms(X1i,X2i)/l_equ;
        for X1i=1:num_X1%不同的ne
            for X2i=1:num_X2%不同的Te
                Hz_p=@(r)Hz0*besselj(0,k_p_wave(X1i,X2i)*r)/besselj(0,k_p_wave(X1i,X2i)*r_p);
                % 与2018Zhao中带FS的FEM模型的20A时磁场结果相比，小了
                Etheta_p=@(r)-1i*k_p_wave(X1i,X2i)*Hz0*besselj(1,k_p_wave(X1i,X2i)*r)/...
                    (w_RF*constants.eps0*epsp_r_real(X1i,X2i)+besselj(0,k_p_wave(X1i,X2i)*r_p));
                % 与2018Zhao中带FS的FEM模型的20A时电场结果相比，大了
                Hz_t=@(r)Hz0; %介质管中磁场均匀
                
                %检验表达式用-断点调试
                % 复坡印廷定理计算复功率
                Etheta_t_at_rc=Etheta_p(r_p)*r_p/r_coil-1i*w_RF*constants.mu0*Hz0*(r_coil^2-r_p^2)/(2*r_coil); %rc处电场
                Pcp1=-pi*r_coil*l_equ*Etheta_t_at_rc*Hz_t(r_coil);
                % ERROR 该表达式与简化表达式结果不一致
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %简化表达式计算复功率
                Pcp_part1=k_p_wave(X1i,X2i)*r_p*besselj(1,k_p_wave(X1i,X2i)*r_p)/...
                    (w_RF*constants.eps0*epsp_r(X1i,X2i)*besselj(0,k_p_wave(X1i,X2i)*r_p));
                Pcp_part2=Pcp_part1+w_RF*constants.mu0*(r_coil^2-r_p^2)/2;
                Pcp=1i*pi*N_equ^2*Icoil_rms(X1i,X2i)^2*Pcp_part2/l_equ;
                
                Pcp_real=real(Pcp); %即Pabs
                Pabs(X1i,X2i)=Pcp_real;
                %                 %检验表达式用-断点调试
                %                 Pcp_real1_part1=1i*k_p_wave(X1i,X2i)*r_p*besselj(1,k_p_wave(X1i,X2i)*r_p)/(w_RF*constants.eps0*epsp_r(X1i,X2i)*besselj(0,k_p_wave(X1i,X2i)*r_p));
                %                 Pcp_real1=pi*N_equ^2*I_coil^2*real(Pcp_real1_part1)/l_equ; %与Pcp_real一致，说明上下表达式自洽
                PER(X1i,X2i)=2*Pcp_real/Icoil_rms(X1i,X2i)^2; %换算到一次侧的等离子体等效电阻
                I_p=l_equ*Hz0*(1/besselj(0,k_p_wave(X1i,X2i)*r_p)-1); %等离子体中电流
                Rp(X1i,X2i)=2*Pcp_real/abs(I_p)^2; %二次侧等离子体电阻
                %                 %检验表达式用-断点调试
                % 与LZS实现结果一致
                %                 Rp1_part1=1i*k_p_wave(X1i,X2i)*r_p*besselj(1,k_p_wave(X1i,X2i)*r_p)/(w_RF*constants.eps0*epsp_r(X1i,X2i)*besselj(0,k_p_wave(X1i,X2i)*r_p));
                % 可能少了一个N^2
                %                 Rp1=2*pi*real(Rp1_part1)/l_equ/(1/besselj(0,k_p_wave(X1i,X2i)*r_p)-1)^2; %与Rp一致，说明上下表达式自洽
                Rsys(X1i,X2i)=PER(X1i,X2i)+Rmetal(X1i,X2i); %系统等效阻抗的电阻分量
                Xsys(X1i,X2i)=2*imag(Pcp)/Icoil_rms(X1i,X2i)^2; %系统等效阻抗的电抗分量
                Lsys(X1i,X2i)=Xsys(X1i,X2i)/w_RF; %系统等效阻抗的电感分量
                Lplasma(X1i,X2i)=Lsys(X1i,X2i)-Lcoil;  % 一次侧的等离子体电感（定义的）
                if ~flag.sweep && Lplasma(X1i,X2i)>0
                    warning('Lplasma应<0, please stop and check')
                    pause
                end
                Xplasma(X1i,X2i)=w_RF*Lplasma(X1i,X2i); %一次侧的等离子体导致电抗
            end
        end
end