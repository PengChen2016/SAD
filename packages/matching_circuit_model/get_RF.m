function [s] = get_RF( type, s, f, Ptrans)
% 单端口-传输线/源与负载的界面 RF各参数之间换算
% 输入s为结构体，不支持矢量化
% ----- mode 1: 复数

% ------ mode 2: 仅模值

s.mode=1;
switch type
    %% mode 1: 复数
    % 先全部转为Zs
    case {0,'','Z','Zcm'}
        assert(isa(s.Z.cm,'double'))
    case {1,'gamma','Gamma','Γ'}
        assert(isa(s.gamma.cm,'double'))
        assert(s.gamma.cm~=0)
        s.Z.cm=s.Z0*(1+s.gamma.cm)/(1-s.gamma.cm);
    case {2,'gamma_ap'}
        assert(s.gamma.abs>=0 && s.gamma.abs<=1)
        phase_rad=s.gamma.phase*pi/180;
        s.gamma.cm=s.gamma.abs*cos(phase_rad)+1i*s.gamma.abs*sin(phase_rad);
        s.Z.cm=s.Z0*(1+s.gamma.cm)/(1-s.gamma.cm);
    case {3,'PfrIrms'}
        % 定向耦合器+电源电流 2021Jain
        assert(s.Pfor>0)
        assert(s.Prel>0)
        gamma_abs=sqrt(s.Prel/s.for);
        s.P=s.Pfor-s.Prel;
        R=P/(s.I.rms)^2;
        Z_abs=sqrt((s.Z0^2*(1-gamma_abs^2)-2*s.Z0*R*(1+gamma_abs^2))/(gamma_abs^2-1));
        X=sqrt(Z_abs^2-R^2);
        s.Z.cm=R+1i*X;
    case {4,'UIphase','UI','UIcm'}
        assert(isa(s.U.cm,'double'))
        assert(isa(s.I.cm,'double'))
        s.Z.cm=s.U.cm/s.I.cm;
    otherwise
        s.mode=2;
end

switch s.mode
    case 1
%         fprintf('defult：时谐系统\n')
        assert(nargin>2)
        assert(s.Z0>0)
        % 换算
        s.Z=get_impedance('Z',s.Z,f);
        s.gamma.cm=(s.Z.cm-s.Z0)/(s.Z.cm+s.Z0);
        s.gamma.abs=abs(s.gamma.cm);
        s.gamma.phase=phase(s.gamma.cm)*180/pi;
        s.S11dB=20*log10(s.gamma.abs);
        s.Prel_ratio=s.gamma.abs^2;
        s.VSWR=(1+s.gamma.abs)/(1-s.gamma.abs);
        
        % 换算功率压流
        switch type
            case {0,'','Z','Zcm',1,'gamma','Gamma','Γ',2,'gamma_ap'}
                % 先转为P
                if exist('Ptrans','var')
                    s.P=Ptrans;
                    s.Pfor=Ptrans/(1-s.Prel_ratio);
                    s.Prel=s.Pfor-Ptrans;
                else
                    warning('no input for Ps.')
                    fprintf('defult：Pfor=1\n')
                    s.Pfor=1;
                    s.Prel=s.Pfor*s.Prel_ratio;
                    s.P=s.Pfor-s.Prel;
                end
                % 已知P、Z
                s=get_UIP( 'PZ', s);
            case {3,'PfrIrms'}
                % 已知P、Z
                assert(nargin<4)
                s.I.abs=sqrt(2)*s.I.rms;
                s=get_UIP( 'IZ', s);
            case {4,'UIphase','UI','UIcm'}
                % 已知U、I
                assert(nargin<4)
                s=get_UIP( 'UIcm', s);
                s.Pfor=s.P/(1-s.Prel_ratio);
                s.Prel=s.Pfor-s.P;
                return
        end
    case 2
        %% mode 2: 仅模值
        assert(nargin>1)
        switch type
            % 先全部转为gamma.abs
            case {'gammaabs','gamma_abs'}
                assert(s.gamma.abs>0)
            case 'Pfr'
                s.gamma.abs=sqrt(s.Prel/s.for);
            case 'S11dB'
                s.gamma.abs=10^(s.S11dB/20);
            case 'VSWR'
                s.gamma.abs=(s.VSWR-1)/(s.VSWR+1);
            otherwise
                error('no such type')
        end
        warning('Do not know: phase(Γ), Z.abs/Z.cm')
        s.gamma.phase=NaN;
        s.Z.phase=NaN;
        s.Z.abs=NaN;
        s.Z.cm=NaN;
        s.Z.R=NaN;
        s.Z.X=NaN;
        s.S11dB=20*log10(s.gamma.abs);
        s.Prel_ratio=s.gamma.abs^2;
        s.VSWR=(1+s.gamma.abs)/(1-s.gamma.abs);
        
        % 换算功率
        switch type
            case 'Pfr'
                s.P=s.Pfor-s.Prel;
            otherwise
                if exist('Ptrans','var')
                    s.P=Ptrans;
                    s.Pfor=Ptrans/(1-s.Prel_ratio);
                    s.Prel=s.Pfor-Ptrans;
                else
                    warning('no input for Ps.')
                    fprintf('defult：Pfor=1\n')
                    s.Pfor=1;
                    s.Prel=s.Pfor*s.Prel_ratio;
                    s.P=s.Pfor-s.Prel;
                end
        end
end


%% 检查值有效性
assert(s.gamma.abs>=0 && s.gamma.abs<=1)

end