function [Z] = get_impedance( type, Z, f )
% ------ mode 1：阻抗各参数之间换算
% 输入Z为结构体，不支持矢量化
% type: 输入结构体应具有的域
% '','Z'：Z.cm
% 'RL'：Z.R, Z.L
% 'RX'：类似'R_L'
% 'ap'：abs和phase，其中phase单位为degree
% ------ mode 2：简单电路计算
% 输入Z为double，支持矢量化
% 'C2X': 由C计算X
% 'X2C': 由X计算C
% 's': 串联。Zlist
% 'p': 并联。Zlist
if nargin==3 
    w=2*pi*f;
end
mode=1;
switch type
    %% mode 1：阻抗各参数之间换算
    case {0,'','Z'}
        assert(isa(Z.cm,'double'))
    case {1,'R_L','RL'}
        assert( nargin==3 )
        assert(isreal(Z.R))
        assert(isreal(Z.L))
        Z.cm=Z.R+1i*w*Z.L;
    case {2,'R_X','RX'}
        assert(isreal(Z.R))
        assert(isreal(Z.X))
        Z.cm=Z.R+1i*Z.X;
    case {3,'a_p','abs_phase','ap'}
        assert(isreal(Z.abs))
        assert(isreal(Z.phase))
        phase_rad=Z.phase*pi/180;
        Z.cm=Z.abs*cos(phase_rad)+1i*Z.abs*sin(phase_rad);
    case {4,'abs_R','R_abs','Ra','aR'}
        assert(isreal(Z.abs))
        assert(isreal(Z.R))
        X=sqrt(Z.abs^2-Z.R^2);
        Z.cm=Z.R+1i*X;
    case {5,'R_C','RC'}
        assert(isreal(Z.R))
        assert(isreal(Z.C))
        Z.cm=Z.R-1i/(w*Z.C);
    otherwise
        mode=2;
end

switch mode
    case 1
        Z.R=real(Z.cm);
        Z.X=imag(Z.cm);
        Z.abs=abs(Z.cm);
        Z.phase=phase(Z.cm)*180/pi;
        Z.Q=Z.X/Z.R; %品质因子
        Z.pF=cos(phase(Z.cm)); %power factor 功率因子
        if nargin==3 
            if Z.X>0
                Z.L=Z.X/w;
            else
                Z.C=-1/(w*Z.X);
            end
        else
            warning('not input f, so no L/C data')
        end
    case 2
        %% mode 2：简单电路计算
        switch type
            case {'C2X','c2x','C','c'}
                % 由C计算X
                % 输入为电容，实数
                assert( nargin==3 )
                assert(isreal(Z))
                Z=-1./(w*Z);
            case {'X2C','x2c'}
                % 由X计算C
                assert( nargin==3 )
                if isreal(Z)
                    Z=-1./(w*Z);
                else
                    Z=-1./(w*imag(Z));
                end
            case {'S','s','series'}
                % 串联
                Z=sum(Z);
            case {'P','p','parallel'}
                % 并联
                Z=1/sum(Z.^-1);
            otherwise
                    error('no such type')
        end
end
end