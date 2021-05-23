function [input, source]=get_magnetized(flag,input,source)
% 磁化处理
% 目前仅适用于CHARLIE

%% 初始化
constants=get_constants();
p0=input.plasma;
size_mat=size(p0.k_wave);
r1=0;
r2=input.geometry.r_plasma_eff;
r=r1:(r2-r1)/100:r2;
len_r=length(r);
Bm=@(H,r) constants.mu0*abs(H(r));

%% get B
if isnumeric(source)
    p0.Bz_mean=source;
else
    if ~isa(source,'struct') || ~isfield(source,'emf')
        source=analytical_EM_model( flag, input_m );
    end
    f_integrand=cell(size_mat);
    B0=zeros(1,len_r);
    for i_r=1:len_r
        temp=Bm(source.emf.Hz_plasma,r(i_r));
        B0(i_r)=temp(end,1);
        % 使用柱坐标系下平均磁场
        for i=1:prod(size_mat)
            f_integrand{i}(i_r)=temp(i)*r(i_r);
        end
    end
    p0.Bz_mean=zeros(size_mat);
    for i=1:prod(size_mat)
        p0.Bz_mean(i)=2*trapz(r,f_integrand{i})/(r2*r2-r1*r1);
    end
end
s0=source;

p0.w_ce=constants.e*p0.Bz_mean/constants.me;
p0.w_ci=constants.e*p0.Bz_mean/constants.mH;

%% 考虑RF磁场的等效电导率
factor=sqrt(1+p0.w_ce.^2./p0.nu_eff.^2);
sigma_dc=p0.ne*constants.e^2./(constants.me*p0.nu_eff);
% 1996Tuszewski
p0.sigma_c=sigma_dc./factor...
    -1i*sigma_dc.*(p0.w_RF./p0.nu_eff)./factor.^3;
p0.sigma=real(p0.sigma_c);
p0.eps_prime=imag(p0.sigma_c)./p0.w_RF;
p0.eps_r=p0.eps_prime/constants.eps0;
p0.eps_c=p0.eps_prime-1i*p0.sigma./p0.w_RF;
p0 = wave_analysis( p0, '' );
p0.magnetized='done';
input.plasma=p0;
source=analytical_EM_model( flag, input );

%% outplot
if flag.output_plasma_model
    % PER
%     if isa(s0,'struct')
%     %     p=[0.3, 0.5, 1, 3, 5, 10]';
%     %
%     %     h_fig=plot_nY(p,'{\itp} [Pa]',...
%     %         source.PER(:,1),'after, 1MHz',...
%     %         s0.PER(:,1),'before, 1MHz',...
%     %         source.PER(:,2),'after, 4MHz',...
%     %         s0.PER(:,2),'before, 4MHz',...
%     %         'PER [\Omega]','plot');
%     elseif isnumeric(s0)
%         
%     end
    
    % 磁场
    B1=zeros(1,len_r);
    for i_r=1:len_r
        temp=Bm(source.emf.Hz_plasma,r(i_r));
        B1(i_r)=temp(end,1);
    end
    
    if isa(s0,'struct')
        h_fig=plot_nY(r,'r [m]',...
            B1,'after',...
            B0,'before',...
            'Bm [T]','plot');
        xticks('auto')
        line([r1,r2],[p0.Bz_mean(end,1), p0.Bz_mean(end,1)])
        legend('after','before','used, before')
    elseif isnumeric(s0)
        h_fig=plot_nY(r,'r [m]',...
            B1,'after',...
            'Bm [T]','plot');
        xticks('auto')
        line([r1,r2],[p0.Bz_mean(end,1), p0.Bz_mean(end,1)])
        legend('after','used, before')
    end
end
end