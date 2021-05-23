function new=test_magnetized(old_input,old_source)
constants=get_constants();
p0=old_input.plasma;
% Bm
Bm=@(H,r) constants.mu0*abs(H(r));
% 使用柱坐标系下平均磁场
r1=0;
r2=45.5e-3;
r=r1:(r2-r1)/100:r2;
len_r=length(r);
f_integrand={};
for i_r=1:len_r
    temp=Bm(old_source.emf.Hz_plasma,r(i_r));
    for j=1:6
        for k=1:2
            f_integrand{j,k}(i_r)=temp(j,k)*r(i_r);
        end
    end
end
for j=1:6
    for k=1:2
        Bz_mean(j,k)=2*trapz(r,f_integrand{j,k})/(r2*r2-r1*r1);
    end
end
w_ce=constants.e*Bz_mean/constants.me;
w_ci=constants.e*Bz_mean/constants.mH;

% Bz_half=Bm(old_source.emf.Hz_plasma,0.5*45.5e-3);
% % 使用边缘磁场
% Bz_edge=Bm(old_source.emf.Hz_plasma,45.5e-3);
% w_ce=constants.e*Bz_edge/constants.me;
% w_ci=constants.e*Bz_edge/constants.mH;


p=[0.3, 0.5, 1, 3, 5, 10]';
% 特征频率
h_fig=plot_nY(p,'{\itp} [Pa]',...
    w_ce(:,1),'wce, 1MHz',...
    w_ce(:,2),'wce, 4MHz',...
    w_ci(:,1),'wci, 1MHz',...
    w_ci(:,2),'wci, 4MHz',...
    p0.nu_eff(:,1),'\nu_{eff}, 1MHz',...
    p0.nu_eff(:,2),'\nu_{eff}, 4MHz',...
    p0.w_RF(:,1),'\omega, 1MHz',...
    p0.w_RF(:,2),'\omega, 4MHz',...
    '\omega [s^{-1}]','loglog');

% 考虑RF磁场的等效电导率
factor=sqrt(1+w_ce.^2./p0.nu_eff.^2);

get_eps_c=@(sigma, eps_r, omega) eps_r*constants.eps0-1i*sigma./omega;

p1=p0;
% 1998Tuszewski  错误
sigma1=p0.sigma_dc.*factor; 
eps_r1=1;
p1.eps_c=get_eps_c(sigma1,eps_r1,p1.w_RF);
p1 = wave_analysis( p1, '' );
i1=old_input;
i1.plasma=p1;
s1=analytical_EM_model( flag, i1 );

p2=p0;
% 根据对1998Tuszewski猜测做调整，影响不大，但理应更优
sigma_c=p0.sigma_c.*factor;
sigma2=real(sigma_c);
eps_r2=imag(sigma_c)./p0.w_RF/constants.eps0;
p2.eps_c=get_eps_c(sigma2,eps_r2,p2.w_RF);
p2 = wave_analysis( p2, '' );
i2=old_input;
i2.plasma=p2;
s2=analytical_EM_model( flag, i2 );

p3=p0;
% 1996Tuszewski
sigma_c=p0.sigma_dc./factor...
    -1i*p0.sigma_dc.*(p0.w_RF./p0.nu_eff)./factor.^3;
sigma3=real(sigma_c);
eps_r3=imag(sigma_c)./p0.w_RF/constants.eps0;
p3.eps_c=get_eps_c(sigma3,eps_r3,p3.w_RF);
p3 = wave_analysis( p3, '' );
i3=old_input;
i3.plasma=p3;
s3=analytical_EM_model( flag, i3 );

h_fig=plot_nY(p,'{\itp} [Pa]',...
    p0.sigma(:,1),'old sigma, 1MHz',...
    p0.sigma(:,2),'old sigma, 4MHz',...
    sigma1(:,1),'new sigma1, 1MHz',...
    sigma1(:,2),'new sigma1, 4MHz',...
    sigma3(:,1),'new sigma3, 1MHz',...
    sigma3(:,2),'new sigma3, 4MHz',...
    '\sigma','loglog');

h_fig=plot_nY(p,'{\itp} [Pa]',...
    s1.PER(:,1),'analytical-new source1',...
    s2.PER(:,1),'analytical-new source2',...
    s3.PER(:,1),'analytical-new source3',...
    old_source.PER(:,1),'analytical-old source',...
    'PER [\Omega]','plot');

len_r=length(r);
B0=zeros(1,len_r);
B1=B0;
B2=B0;
B3=B0;
for i_r=1:len_r
    temp=Bm(old_source.emf.Hz_plasma,r(i_r));
    B0(i_r)=temp(end,1);
    temp=Bm(s1.emf.Hz_plasma,r(i_r));
    B1(i_r)=temp(end,1);
    temp=Bm(s2.emf.Hz_plasma,r(i_r));
    B2(i_r)=temp(end,1);
    temp=Bm(s3.emf.Hz_plasma,r(i_r));
    B3(i_r)=temp(end,1);
end

h_fig=plot_nY(r,'r [m]',...
    B1,'analytical-new source1',...
    B2,'analytical-new source2',...
    B3,'analytical-new source3',...
    B0,'analytical-old source',...
    'Bm [T]','plot');
xticks('auto') 
line([r1,r2],[Bz_mean(end,1) Bz_mean(end,1)])

new.i1=i1;
new.i2=i2;
new.s1=s1;
new.s2=s2;
new.factor=factor;
end

% input0=input;
% input0.plasma=plasma_raza;
% new=test_magnetized(input0,source_a);
% new_za=new;
% 
% load('CHARLIE_razcoil_sweep210428.mat')
% input0=input;
% s0=analytical_EM_model( flag, input0);
% new=test_magnetized(input0,s0);
% 
% h_fig=plot_nY(p,'{\itp} [Pa]',...
%     experiment.PER(:,1),'subtractive method',...
%     fem.dielectric_PER(:,1),'FEM model',...
%     source_a.PER(:,1),'analytical EM-za',...
%     s0.PER(:,1),'analytical EM-zcoil',...
%     new_za.s1.PER(:,1),'analytical EM-za-magnetized',...
%     new.s1.PER(:,1),'analytical EM-zcoil-magnetized',...
%     'PER [\Omega]','plot');
% 
% % 迭代
% % za-source 1做迭代
% new1=test_magnetized(new_za.i1,new_za.s1);
% new1=test_magnetized(new1.i1,new1.s1);
% new1=test_magnetized(new1.i1,new1.s1);
% % za-source 2做迭代
% new2=test_magnetized(new_za.i2,new_za.s2);
% new2=test_magnetized(new2.i2,new2.s2);
% new2=test_magnetized(new2.i2,new2.s2);
% 
% h_fig=plot_nY(p,'{\itp} [Pa]',...
%     experiment.PER(:,1),'subtractive method',...
%     fem.dielectric_PER(:,1),'FEM model',...
%     source_a.PER(:,1),'analytical EM-za',...
%     s0.PER(:,1),'analytical EM-zcoil',...
%     new_za.s1.PER(:,1),'analytical EM-za-magnetized-1',...
%     new_za.s2.PER(:,1),'analytical EM-za-magnetized-2',...
%     new1.s1.PER(:,1),'analytical EM-za-magnetized-iter1',...
%     new2.s2.PER(:,1),'analytical EM-za-magnetized-iter2',...
%     'PER [\Omega]','plot');
% title('1MHz')
% h_fig=plot_nY(p,'{\itp} [Pa]',...
%     experiment.PER(:,2),'subtractive method',...
%     fem.dielectric_PER(:,2),'FEM model',...
%     source_a.PER(:,2),'analytical EM-za',...
%     s0.PER(:,2),'analytical EM-zcoil',...
%     new_za.s1.PER(:,2),'analytical EM-za-magnetized-1',...
%     new_za.s2.PER(:,2),'analytical EM-za-magnetized-2',...
%     new1.s1.PER(:,2),'analytical EM-za-magnetized-iter1',...
%     new2.s2.PER(:,2),'analytical EM-za-magnetized-iter2',...
%     'PER [\Omega]','plot');
% title('4MHz')