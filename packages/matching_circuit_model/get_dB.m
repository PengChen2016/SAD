function [ out_vec ] = get_dB( type,in_vec )
% dB换算
switch type
    %% 正向
    case {1,'','P','power'}
        out_vec=10*log10(in_vec);
    case {'dBm'}
        out_vec=10*log10(in_vec*1e3);
    case {2,'U','I','UI'}
        out_vec=20*log10(in_vec);
        %% 逆运算
    case {'i','P_inverse','power_inverse','P_i','Pi'}
        out_vec=10^(in_vec/10);
    case {'dBm_inverse','dBmi'}
        out_vec=1e-3*10^(in_vec/10);
    case {'U_inverse','I_inverse','UI_inverse','Ui','Ii','UIi'}
        out_vec=10^(in_vec/20);
    otherwise
        error('no such type')
end
end