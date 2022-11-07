% -*- coding: utf-8 -*-
% ----------------------------------------------
%{
 * brief #Abstract  
 * Created 13:10:30 2022/11/07
 * author 
   Zengshan Li, HUST: origin 新的可用的整体模型代码.rar myOde2 2020/09/03
 * note #Detail
    自变量不含时，但输出结果为时间微分量
 * #TODO
%}
% ----------------------------------------------
function [Res]=global_model_H2_LZS2020(t,X,Pabs,Press,Tg,Q0,R,L)
elec=1.602E-19;
Res=zeros(7,1);
nHi=X(1);
nH2i=X(2);
nH3i=X(3);
nHNi=X(4);
ne=X(1)+X(2)+X(3)-X(4);
Te=X(5)./elec./ne/1.5;
nH=X(6);
nHH=X(7);

[nHiDt1,nH2iDt1,nH3iDt1,nHNiDt1,DTe1,nHDt1,nHHDt1]...
    =Dn(Pabs,Press,R,L,Tg,Q0,nHi,nH2i,nH3i,nHNi,ne,Te,nH,nHH);


Res(1)=nHiDt1; % 时间微分
Res(2)=nH2iDt1;
Res(3)=nH3iDt1;
Res(4)=nHNiDt1;
Res(5)=DTe1;
Res(6)=nHDt1;
Res(7)=nHHDt1;



end

function [nHiDt,nH2iDt,nH3iDt,nHNiDt,DTe,nHDt,nHHDt]...
    =Dn(Pabs,Press,R,L,Tg,Q0,nHi,nH2i,nH3i,nHNi,ne,Te,nH,nHH)

nHiDt=0;
nH2iDt=0;
nH3iDt=0;
nHNiDt=0;
DTe=0;
nHDt=0;
nHHDt=0;

e=1.602E-19;

%反应1.1，H弹性，无粒子净通量
name='R1d1.mat';
me=9.10956E-31;
mH=1.674E-27;
Eth=3*me/mH*Te;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nH*Rate*e*Eth;

%1.2.1
name='R1d2d1.mat';
Eth=10.2;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nH*Rate*e*Eth;

%1.2.2
name='R1d2d2.mat';
Eth=12.09;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nH*Rate*e*Eth;

%1.2.3
name='R1d2d3.mat';
Eth=12.75;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nH*Rate*e*Eth;

%1.3
name='R1d3.mat';
Eth=13.6;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nH*Rate*e*Eth;
nHDt=nHDt-ne*nH*Rate;
nHiDt=nHiDt+ne*nH*Rate;

%反应2.1
name='R2d1.mat';
me=9.10956E-31;
mH=1.674E-27;
Eth=3*me/mH/2*Te;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.2.1
name='R2d2d1.mat';
Eth=0.516;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.2.2
name='R2d2d2.mat';
Eth=1.003;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.1
name='R2d3d1.mat';
Eth=12.42;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.2
name='R2d3d2.mat';
Eth=14.62;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.3
name='R2d3d3.mat';
Eth=13.84;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.4
name='R2d3d4.mat';
Eth=11.37;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.5
name='R2d3d5.mat';
Eth=12.41;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.6
name='R2d3d6.mat';
Eth=14.74;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.7
name='R2d3d7.mat';
Eth=14.13;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.8
name='R2d3d8.mat';
Eth=11.89;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.9
name='R2d3d9.mat';
Eth=13.98;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.10
name='R2d3d10.mat';
Eth=13.36;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.11
name='R2d3d11.mat';
Eth=13.98;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%反应2.3.12
name='R2d3d12.mat';
Eth=11.9;
[Rate]=eleColl(name,Te);
DTe=DTe-ne*nHH*Rate*e*Eth;

%2.4
name='R2d4.mat';
Eth=10;
[Rate]=eleColl(name,Te);
Sulv=ne*nHH*Rate;
DTe=DTe-Sulv*e*Eth;
nHDt=nHDt+2*Sulv;
nHHDt=nHHDt-Sulv;

%2.5
name='R2d5.mat';
Eth=15.4;
[Rate]=eleColl(name,Te);
Sulv=ne*nHH*Rate;
DTe=DTe-Sulv*e*Eth;
nH2iDt=nH2iDt+Sulv;
nHHDt=nHHDt-Sulv;

%2.6.1
name='R2d6d1.mat';
Eth=18.15;
[Rate]=eleColl(name,Te);
Sulv=ne*nHH*Rate;
DTe=DTe-Sulv*e*Eth;
nHiDt=nHiDt+Sulv;
nHDt=nHDt+Sulv;
nHHDt=nHHDt-Sulv;

%2.6.2
name='R2d6d2.mat';
Eth=30.6;
[Rate]=eleColl(name,Te);
Sulv=ne*nHH*Rate;
DTe=DTe-Sulv*e*Eth;
nHiDt=nHiDt+Sulv;
nHDt=nHDt+Sulv;
nHHDt=nHHDt-Sulv;

%2.7
name='R2d7.mat';
Eth=3.72;
[Rate]=eleColl(name,Te);
Sulv=ne*nHH*Rate;
DTe=DTe-Sulv*e*Eth;
nHNiDt=nHNiDt+Sulv;
nHDt=nHDt+Sulv;
nHHDt=nHHDt-Sulv;

%3.1
name='R3d1.mat';
Eth=0;
[Rate]=eleColl(name,Te);
Sulv=ne*nHi*Rate;
DTe=DTe-Sulv*e*Eth;
nHNiDt=nHNiDt-Sulv;
nHDt=nHDt+Sulv;

%4.1
name='R4d1.mat';
Eth=2.4;
[Rate]=eleColl(name,Te);
Sulv=ne*nH2i*Rate;
DTe=DTe-Sulv*e*Eth;
nHiDt=nHiDt+Sulv;
nHDt=nHDt+Sulv;
nH2iDt=nH2iDt-Sulv;

%4.2
name='R4d2.mat';
Eth=0;
[Rate]=eleColl(name,Te);
Sulv=ne*nH2i*Rate;
DTe=DTe-Sulv*e*Eth;
nHDt=nHDt+2*Sulv;
nH2iDt=nH2iDt-Sulv;

%5.2
name='R5d2.mat';
Eth=14;
[Rate]=eleColl(name,Te);
Sulv=ne*nH3i*Rate;
DTe=DTe-Sulv*e*Eth;
nHiDt=nHiDt+Sulv;
nHDt=nHDt+2*Sulv;
nH3iDt=nH3iDt-Sulv;

%5.3
name='R5d3.mat';
Eth=0;
[Rate]=eleColl(name,Te);
Sulv=ne*nH3i*Rate;
DTe=DTe-Sulv*e*Eth;
nHDt=nHDt+3*Sulv;
nH3iDt=nH3iDt-Sulv;

%5.4
name='R5d4.mat';
Eth=0;
[Rate]=eleColl(name,Te);
Sulv=ne*nH3i*Rate;
DTe=DTe-Sulv*e*Eth;
nHHDt=nHHDt+Sulv;
nHDt=nHDt+Sulv;
nH3iDt=nH3iDt-Sulv;

%6.1
name='R6d1.mat';
Eth=0.754;
[Rate]=eleColl(name,Te);
Sulv=ne*nHNi*Rate;
DTe=DTe-Sulv*e*Eth;
nHDt=nHDt+Sulv;
nHNiDt=nHNiDt-Sulv;

%%%%% 下面是重粒子碰撞

%7.1
[Rate]=R7d1(Tg);
Sulv=nH*nHNi*Rate;
nHDt=nHDt-Sulv;
nHNiDt=nHNiDt-Sulv;
nHHDt=nHHDt+Sulv;

%7.4.1
[Rate]=R7d4d1(Tg);
Sulv=nHi*nHNi*Rate;
nHiDt=nHiDt-Sulv;
nHNiDt=nHNiDt-Sulv;
nHDt=nHDt+2*Sulv;

%7.4.2
[Rate]=R7d4d2(Tg);
Sulv=nHi*nHNi*Rate;
nHiDt=nHiDt-Sulv;
nHNiDt=nHNiDt-Sulv;
nHDt=nHDt+2*Sulv;

%7.5
[Rate]=R7d5(Tg);
Sulv=nH2i*nHNi*Rate;
nH2iDt=nH2iDt-Sulv;
nHNiDt=nHNiDt-Sulv;
nHDt=nHDt+3*Sulv;

%7.6
[Rate]=R7d6(Tg);
Sulv=nH2i*nHNi*Rate;
nH2iDt=nH2iDt-Sulv;
nHNiDt=nHNiDt-Sulv;
nHDt=nHDt+Sulv;
nHHDt=nHHDt+Sulv;

%7.9
[Rate]=R7d9(Tg);
Sulv=nH3i*nHNi*Rate;
nH3iDt=nH3iDt-Sulv;
nHNiDt=nHNiDt-Sulv;
nHDt=nHDt+2*Sulv;
nHHDt=nHHDt+Sulv;

%7.10
[Rate]=R7d10(Tg);
Sulv=nH2i*nH*Rate;
nH2iDt=nH2iDt-Sulv;
nHDt=nHDt-Sulv;
nHiDt=nHiDt+Sulv;
nHHDt=nHHDt+Sulv;

%7.11
[Rate]=R7d11(Tg);
Sulv=nH2i*nHH*Rate;
nH2iDt=nH2iDt-Sulv;
nHHDt=nHHDt-Sulv;
nH3iDt=nH3iDt+Sulv;
nHDt=nHDt+Sulv;

%%%%%% 下面是器壁反应

uBHi=sqrt(e*Te/mH);
uBH2i=sqrt(e*Te/mH/2);
uBH3i=sqrt(e*Te/mH/3);
Ethe=2*Te;
Ethi=0.5*Te;
Vs=-log(4*(nHi*uBHi+nH2i*uBH2i+nH3i*uBH3i)/ne...
    *sqrt(pi*me/8/e/Te))*Te;

%8.1 H+和器壁
HiWallRate=WallHi(R,L,nH,nHH,uBHi);
Sulv=HiWallRate*nHi;
DTe=DTe-Sulv*e*(Ethe+Ethi+Vs);
nHiDt=nHiDt-Sulv;
nHDt=nHDt+Sulv;

%8.2 H2+和器壁
H2iWallRate=WallH2i(R,L,nH,nHH,uBH2i);
Sulv=H2iWallRate*nH2i;
DTe=DTe-Sulv*e*(Ethe+Ethi+Vs);
nH2iDt=nH2iDt-Sulv;
nHHDt=nHHDt+Sulv;

%8.3 H3+和器壁
H3iWallRate=WallH3i(R,L,nH,nHH,uBH3i);
Sulv=H3iWallRate*nH3i;
DTe=DTe-Sulv*e*(Ethe+Ethi+Vs);
nH3iDt=nH3iDt-Sulv;
nHHDt=nHHDt+Sulv;
nHDt=nHDt+Sulv;

%8.4 H器壁复合
HWallRate=WallH(Tg,R,L,nH,nHH);
nHDt=nHDt-nH*HWallRate;
nHHDt=nHHDt+0.5*nH*HWallRate;

%进气过程
Volume=pi*R*R*L;
nHHDt=nHHDt+4.48E17*Q0/Volume;

%抽气过程
%Pout的单位为Torr

kB=1.3806505E-23;
nTot=Press/Tg/kB;
Pout0=nTot*1.27E-5/4.48E17;
Pout=Pout0;
% Pout=(nHH+nH+nHi+nH2i+nH3i+nHNi)/nTot*Pout0;

nHHDt=nHHDt-nHH*1.27E-5*Q0/Volume/Pout;
nHDt=nHDt-nH*1.27E-5*Q0/Volume/Pout;
nHiDt=nHiDt-nHi*1.27E-5*Q0/Volume/Pout;
nH2iDt=nH2iDt-nH2i*1.27E-5*Q0/Volume/Pout;
nH3iDt=nH3iDt-nH3i*1.27E-5*Q0/Volume/Pout;
% nHHDt=0;

%功率平衡
DTe=Pabs/Volume+DTe;



end
function Rate=WallH(Tg,R,L,nH,nHH)
    kB=1.3806505E-23;       %J/K
    mH=1.674E-27;
    diamH=2*0.053E-9;
    diamHH=4*0.053E-9;
    diamH_HH=(diamH+diamHH)/2;
    DH_H=3/32/(nH)/diamH/diamH*sqrt(8*kB*Tg/pi*(1/mH+1/mH));
    DH_HH=3/32/(nHH)/diamH_HH/diamH_HH*sqrt(8*kB*Tg/pi*(1/mH+0.5/mH));

    diffLength=1/sqrt((pi/L)*(pi/L)+(2.405/R)*(2.405/R));
    V=pi*R*R*L;
    A=2*pi*R*R+2*pi*R*L;
    gammRec=0.1;
    v=sqrt(8*kB*Tg/pi/mH);
    Dn=1/(1/(diffLength*v/3)+1/DH_H+1/DH_HH);
    Rate=1/(diffLength*diffLength/Dn+...
        2*V*(2-gammRec)/A/v/gammRec);   %单位1/s
end

function Rate=WallH3i(R,L,nH,nHH,uBH3i)
    sigmH3i_H=28E-20;
    sigmH3i_HH=33E-20;
    Lamdai=1./(sigmH3i_HH*nHH+nH*sigmH3i_H);        %计算平均自由程
    hL=0.86/sqrt(3+0.5*L/Lamdai);
    hR=0.8/sqrt(4+R/Lamdai);
    Aeff=2*pi*R*L*hR+2*pi*R*R*hL;
    Volume=pi*R*R*L;
    Rate=Aeff/Volume*uBH3i;
end

function Rate=WallH2i(R,L,nH,nHH,uBH2i)
    sigmH2i_H=23E-20;
    sigmH2i_HH=28E-20;
    Lamdai=1./(sigmH2i_HH*nHH+nH*sigmH2i_H);        %计算平均自由程
    hL=0.86/sqrt(3+0.5*L/Lamdai);
    hR=0.8/sqrt(4+R/Lamdai);
    Aeff=2*pi*R*L*hR+2*pi*R*R*hL;
    Volume=pi*R*R*L;
    Rate=Aeff/Volume*uBH2i;
end

function Rate=WallHi(R,L,nH,nHH,uBHi)
    sigmHi_H=18E-20;
    sigmHi_HH=23E-20;
    Lamdai=1./(sigmHi_HH*nHH+nH*sigmHi_H);        %计算平均自由程
    hL=0.86/sqrt(3+0.5*L/Lamdai);
    hR=0.8/sqrt(4+R/Lamdai);
    Aeff=2*pi*R*L*hR+2*pi*R*R*hL;
    Volume=pi*R*R*L;
    Rate=Aeff/Volume*uBHi;
end

function Rate=R7d11(Tg)
    Rate=2.1E-15;
end

function Rate=R7d10(Tg)
    Rate=6.4E-16;
end

function Rate=R7d9(Tg)
    Rate=2E-7*(300/Tg)*1E-6;
end

function Rate=R7d6(Tg)
    Rate=2E-7*(300./Tg)^(0.5)*1E-6;
end

function Rate=R7d5(Tg)
    Rate=8.29E-13*(300/Tg)^0.5;
end

function Rate=R7d4d2(Tg)
    Rate=1.77E-7*(300./Tg)^(0.5)*1E-6;
end

function Rate=R7d4d1(Tg)
    Rate=9.1E-11.*(300./Tg).^(0.83).*1E-6;
end

function Rate=R7d1(Tg)
    Rate=1.3E-9*1E-6;
end

function [Rate]=eleColl(name,Te)
   
    s=load(name);  
    sc=struct2cell(s);  
    Sulv=cell2mat(sc); 

    if Te<Sulv(1,1)
        Rate=Sulv(1,2);
    elseif Te>Sulv(end,1)
        Rate=Sulv(end,2);
    else
        x0=log10(Sulv(:,1));
        y0=log10(Sulv(:,2));
        Rate=10.^interp1(x0,y0,log10(Te));
    end

end









