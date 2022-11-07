function [s] = get_UIP( type, s, in_phase)
% 计算压流功率
% 输入结构体s，其成员值默认为cm或abs。
assert(nargin>1)
switch type
    case {'UI','UIcm','UIphase'}
        if nargin==2
%             fprintf('已输入：UIcm\n')
            s.U.abs=abs(s.U.cm);
            s.U.phase=phase(s.U.cm)*180/pi;
            s.I.abs=abs(s.I.cm);
            s.I.phase=phase(s.I.cm)*180/pi;
%             fprintf('phase(I)|t=0 = %.2f °\n',s.I.phase)
        elseif nargin==3
            assert(s.U.abs>0);
            assert(s.I.abs>0);
%             fprintf('已输入：UIabs and Zphase\n')
%             fprintf('defult：phase(I)|t=0 = 0\n')
            s.I.phase=0;
            s.I.cm=s.I.abs;
            s.U.phase=s.I.phase+in_phase;
            phase_rad=s.U.phase*pi/180;
            s.U.cm=s.U.abs*cos(phase_rad)+1i*s.U.abs*sin(phase_rad);
        else
            error('Input error')
        end        
    case {'PZ','IZ'}
        % 以下先转为Iabs，然后考虑Zcm与初相得到Icm和Ucm
        switch type
            case 'PZ'
                assert(s.P>0)
                s.I.abs=sqrt(2*s.P/real(s.Z.cm));
                %     case {'PI','PIcm'}  % 无法得到Z信息
            case 'IZ'
                assert(s.I.abs>0)
        end
        if nargin==2
%             fprintf('defult：phase(I)|t=0 = 0\n')
            s.I.phase=0;
            s.I.cm=s.I.abs;
        elseif nargin==3
            s.I.phase=in_phase;
%             fprintf('phase(I)|t=0 = %.2f °\n',s.I.phase)
            phase_rad=in_phase*pi/180;
            s.I.cm=s.I.abs*cos(phase_rad)+1i*s.I.abs*sin(phase_rad);
        end
        s.U.cm=s.I.cm*s.Z.cm;
        s.U.abs=abs(s.U.cm);
        s.U.phase=phase(s.U.cm)*180/pi;
end

s.I.rms=s.I.abs/sqrt(2);
s.U.rms=s.U.abs/sqrt(2);
PQ=s.U.cm*conj(s.I.cm)/2;
s.P=real(PQ);
s.Q=imag(PQ);
% s.P=s.U.rms*s.I.rms*cos(phase(s.U.cm/s.I.cm));
% s.Q=s.U.rms*s.I.rms*sin(phase(s.U.cm/s.I.cm));
% s.P=s.U.rms*s.I.rms*s.Z.pF;
% s.Q=sign(s.Z.phase)*s.Z.Q*s.P;
end