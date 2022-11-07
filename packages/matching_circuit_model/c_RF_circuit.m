% -*- coding: utf-8 -*-
% ----------------------------------------------
%{
 * brief #Abstract 用于简单RF电路模型的计算类
 * Created 18:12:41 2022/04/29
 * author PengChen, HUST(peng_chen2016@hust.edu.cn)
 *
 * note #Detail
1. Matlab是一种伪面向对象，不支持嵌套类。主要是不能调用成员变量的成员函数
2. 理想变压器模型对应代码参考2022李嘉旺，互感耦合电路模型对应代码参考2021周磊
 * #TODO
目前只适用于Γ型，以后扩展至其他形状
%}
% ----------------------------------------------
classdef c_RF_circuit
    properties
        f
        w
        % ----------- power source/matching unit port
        s % struct，见构造函数
        % ----------- matching unit
        m % struct，见构造函数
        % -----------load: driver branch
        D % struct，见构造函数
        Dlist % struct，见构造函数
    end
    
    methods
        %% 构造函数
        function obj = c_RF_circuit(f)
            % 20220430 test ok
            if nargin<1
                obj.f=1e6;
                fprintf('defult：f=%d\n',obj.f)
            elseif nargin>2
                error('too much input')
            else
                obj.f=f;
            end
            obj.w=2*pi*obj.f;
            % ----------- power source/matching unit port: s port
            % ------ basic electric parameters
            s0.name='power source/matching unit port';
            s0.P=0;
            s0.Q=0;
            s0.Pfor=0;
            s0.Prel=0;
            s0.U.cm=0;
            s0.U.abs=0;
            s0.U.phase=0;
            s0.U.rms=0;
            s0.I=s0.U;
            % ------ impedance
            s0.Z0=50; % 功率源/传输线特征阻抗
            s0.Z.cm=0; % 端口负载阻抗
            s0.Z.R=0;
            s0.Z.X=0;
            s0.Z.L=0;
            s0.Z.abs=0;
            s0.Z.phase=0; % in degree
            % ------ RF parameters
            s0.gamma.cm=0;
            s0.gamma.abs=0;
            s0.gamma.phase=0;
            s0.VSWR=0;
            s0.Prel_ratio=0;
            
            obj.s=s0;
            % ----------- matching unit
            % ------ C
            % 理想电容，支路上电阻以后文方法等效考虑
            m0.name='matching unit';
            m0.s.Z.C=0;
            m0.s.Z.X=0;
            % ------ basic electric parameters at branchs
            m0.s.U=s0.U;
            m0.s.I=s0.I;
            
            m0.p=m0.s;
            % ------ transformer
            m0.t.name='RF transformer';
            m0.t.model='';
            % for ideal transformer
            m0.t.n=1;
            % for mutual inductance coupling circuit
            m0.t.M=0;
            m0.t.k=0;
            m0.t.L1=0;
            m0.t.L2=0;
            % ------ 损耗
            m0.R=0; % 折算到ID_rms的匹配单元损耗
            m0.P=0;% 匹配单元损耗
            m0.R1=0; % 当存在变压器时，一次侧等效电阻
            m0.R2=0; % 当存在变压器时，二次侧等效电阻
            
            %             m0.slist=struct([]);
            %             m0.plist=struct([]);
            
            obj.m=m0;
            % -----------load: driver branch
            D0.Z=s0.Z;
            D0.num=0;
            % ------ basic electric parameters at driver port
            D0.U=s0.U;
            D0.I=s0.I;
            D0.P=0;
            
            obj.D=D0;
            
            obj.Dlist=struct([]);
        end
        
        %% initialization
        function obj = init_source_port( obj ,type, s, Ptrans, Z0)
            % s应具有完整的RF参数，如s.Z.cm
            if nargin==5
                obj.s.Z0=Z0;
            else
                obj.s.Z0=50;
            end
            
            switch nargin
                case 2
                    % 使用内部s
                    if obj.s.P>0
                        % 使用内部P
                        obj.s=get_RF( type, obj.s, obj.f, obj.s.P);
                    else
                        % 无P（使用默认值）
                        obj.s=get_RF( type, obj.s, obj.f);
                    end
                case 3
                    % 使用外部s
                    s.Z0=obj.s.Z0;
                    if s.P>0
                        % 使用内部P
                        obj.s=get_RF( type, s, obj.f, s.P);
                    else
                        % 无P（使用默认值）
                        obj.s=get_RF( type, s, obj.f);
                    end
                case {4,5}
                    % 使用外部s和外部P
                    s.Z0=obj.s.Z0;
                    obj.s=get_RF( type, s, obj.f, Ptrans);
                otherwise
                    error('input error')
            end
        end
        
        function obj = init_matching_C( obj, type, Cs, Cp, Cs_UI_withstand, Cp_UI_withstand)
            %
            assert(nargin>3)
            switch type
                case {'','Gamma','gamma','Γ','L'}
                    if nargin<5
                        obj.m.s=obj.init_C(Cs);
                        obj.m.p=obj.init_C(Cp);
                    else
                        obj.m.s=obj.init_C(Cs,Cs_UI_withstand);
                        obj.m.p=obj.init_C(Cp,Cp_UI_withstand);
                    end
                    obj.m.type='gamma';
                otherwise
                    error('TODO')
            end
        end
        
        function [out_struct] = init_C( obj, C, UI_withstand)
            % out_struct 类似obj.m.s/p的结构
            assert(nargin>1)
            num_C=length(C);
            out_struct=struct;
            if num_C==1
                out_struct.Z.cm=1i*get_impedance( 'C2X', C, obj.f );
                out_struct.Z=get_impedance( 'Z', out_struct.Z, obj.f );
                if exist('UI_withstand','var')
                    out_struct.U_withstand=UI_withstand(1);
                    out_struct.I_withstand=UI_withstand(2);
                end
            else
                % 输入多个并联电容，处理得到容值与耐压耐流
                Xlist=get_impedance( 'C2X', C, obj.f );
                out_struct.Z.cm=1i*get_impedance( 'p', Xlist, obj.f );
                out_struct.Z=get_impedance( 'Z', out_struct.Z, obj.f );
                if exist('UI_withstand','var')
                    out_struct.U_withstand=min(UI_withstand(:,1));
                    if size(Xlist,2)>1
                        Xlist=Xlist';
                    end
                    out_struct.I_withstand=min(UI_withstand(:,2).*Xlist/out_struct.Z.X);
                end
            end
        end
        
        function obj = init_matching_t( obj, model, NorM, L1, L2)
            %
            assert(nargin>1)
            switch model
                case ''
                    % 无RF变。
                    assert(nargin==2)
                    % 清空数据
                    % for ideal transformer
                    obj.m.t.n=1;
                    % for mutual inductance coupling circuit
                    obj.m.t.M=0;
                    obj.m.t.k=0;
                    obj.m.t.L1=0;
                    obj.m.t.L2=0;
                    obj.m.t.model='';
                case {'ideal transformer','ideal','i'}
                    % 理想变压器模型
                    assert(nargin==3)
                    obj.m.t.n=NorM;
                    obj.m.t.k=1;
                    obj.m.t.model='ideal';
                case {'mutual inductance coupling circuit','mutual','m'}
                    % 互感耦合电路模型
                    assert(nargin==5)
                    if NorM>1e-2
                        warning('input M>1e-2 H. Too large.')
                    end
                    obj.m.t.M=NorM;
                    obj.m.t.L1=L1;
                    obj.m.t.L2=L2;
                    obj.m.t.k=NorM/sqrt(L1*L2);
                    obj.m.t.n=sqrt(L1/L2);
                    obj.m.t.model='mutual';
                case {'mutual2','m2'}
                    error('请以mutual类型初始化，需要get_C时再临时手动设置 obj.m.t.model=''mutual2'' ')
                otherwise
                    obj = obj.init_transformer('');
                    fprintf('defult：No RF transformer\n')
            end
        end
        
        function obj = init_matching_R(obj, R1, R2)
            assert(nargin>0)
            obj.m.R1=R1;
            if exist('R2','var')
                obj.m.R2=R2;
            else
                obj.m.R2=0;
            end
        end
        
        function obj = init_ZD( obj, type, num_drivers, ZorZlist, ZD_leadline, f, configuration_drivers )
            % ZorZlist：单驱时为Z，结构体；多驱时为Zlist，结构体数组
            % TODO: ZD_leadline/Zcoi_lleadline
            % 输入检查
            assert(nargin>1)
            if exist('num_drivers','var')
                obj.D.num=num_drivers;
            elseif obj.D.num<1
                obj.D.num=1;
                fprintf('default: 单驱\n')
            end
            if exist('ZD_leadline','var')
                if isa(ZD_leadline,'struct')
                    assert(ZD_leadline.cm>=0)
                    obj.D.leadline.Z=ZD_leadline;
                elseif isa(ZD_leadline,'double') && ZD_leadline>=0
                    obj.D.leadline.Z.cm=ZD_leadline;
                else
                    obj.D.leadline.Z.cm=0;
                end
            else
                obj.D.leadline.Z.cm=0;
            end
            if exist('f','var')
                obj.f=f;
            elseif obj.f<1
                obj.f=1e6;
                fprintf('default: f=1MHz\n')
            end
            
            % 处理
            switch obj.D.num
                case 1
                    % single driver
                    if ~exist('ZorZlist','var')
                        % 使用内部值
                        ZorZlist=obj.D.Z;
                    end
                    obj.D.Z=get_impedance(type,ZorZlist,obj.f);
                    obj.Dlist=struct([]);
                otherwise
                    % multi driver
                    if exist('configuration_drivers','var')
                        obj.D.configuration=configuration_drivers;
                    else
                        obj.D.configuration='s';
                        fprintf('default：多驱串联\n')
                    end
                    % get Dlist.Z
                    if exist('ZorZlist','var')
                        % 使用外部输入值
                        assert(length(ZorZlist)==obj.D.num)
                        obj.Dlist=struct([]);
                        for i=1:obj.D.num
                            obj.Dlist(i).Z=get_impedance(type,ZorZlist(i),obj.f);
                        end
                    else
                        % 使用内部值
                        for i=1:obj.D.num
                            obj.Dlist(i).Z=get_impedance(type,obj.Dlist(i).Z,obj.f);
                        end
                    end
                    % 已知obj.Dlist.Z, 计算obj.D.Z
                    Zcm=0;
                    switch obj.D.configuration
                        case 's'
                            for i=1:obj.D.num
                                Zcm=Zcm+obj.Dlist(i).Z.cm;
                            end
                        case 'p'
                            for i=1:obj.D.num
                                Zcm=Zcm+1/obj.Dlist(i).Z.cm;
                            end
                            Zcm=1/Zcm;
                    end
                    obj.D.Z.cm=Zcm+obj.D.leadline.Z.cm;
                    obj.D.Z=get_impedance('Z',obj.D.Z,obj.f);
            end
        end
        
        function obj = get_Zs( obj )
            % 已知ZD与匹配电路模型，计算Zs
            assert(obj.D.Z.cm>0) % ZD有实部>0
            Z_DCs=obj.D.Z.cm+obj.m.s.Z.cm;
            switch obj.m.t.model
                case ''
                    Z1=Z_DCs;
                case 'ideal'
                    Z1=(obj.m.t.n^2)*Z_DCs;
                case 'mutual'
                    ZL1=1i*obj.w*obj.m.t.L1;
                    ZL2=1i*obj.w*obj.m.t.L2;
                    ZM=1i*obj.w*obj.m.t.M;
                    if ~(obj.m.R1>=0 && obj.m.R2>=0)
                        warning('R1 or  R2 < 0. Loss in the matching unit is ignored.')
                        obj.m.R1=0;
                        obj.m.R2=0;
                    end
                    if obj.m.R1==0 && obj.m.R2==0
                        % Lossless matching unit.
                        Z1=ZL1-ZM^2/(Z_DCs+ZL2);
                    else
                        % Lossy matching unit.
                        Z1=obj.m.R1+ZL1-ZM^2/(Z_DCs+ZL2+obj.m.R2);
                    end
                otherwise
                    error('No such type')
            end
            obj.s.Z.cm=get_impedance('parallel',[Z1,obj.m.p.Z.cm]);
            obj.s=get_RF( 'Z', obj.s, obj.f );
        end
        
        function obj = get_ZD( obj )
            % 已知Zs与匹配电路模型，计算ZD
            Z_sCp=-get_impedance('parallel',[obj.m.p.Z.cm,-obj.s.Z.cm]);
            switch obj.m.t.model
                case ''
                    Z2=Z_sCp;
                case 'ideal'
                    Z2=Z_sCp/(obj.m.t.n^2);
                case 'mutual'
                    ZL1=1i*obj.w*obj.m.t.L1;
                    ZL2=1i*obj.w*obj.m.t.L2;
                    ZM=1i*obj.w*obj.m.t.M;
                    if ~(obj.m.R1>=0 && obj.m.R2>=0)
                        warning('R1 or  R2 < 0. Loss in the matching unit is ignored.')
                        obj.m.R1=0;
                        obj.m.R2=0;
                    end
                    if obj.m.R1==0 && obj.m.R2==0
                        % Lossless matching unit.
                        Z2=ZM^2/(ZL1-Z_sCp)-ZL2;
                    else
                        % Lossy matching unit.
                        Z2=ZM^2/(obj.m.R1+ZL1-Z_sCp)-ZL2-obj.m.R2;
                    end
                otherwise
                    error('No such type')
            end
            obj.D.Z.cm=Z2 - obj.m.s.Z.cm;
            obj.D.Z=get_impedance( 'Z', obj.D.Z, obj.f );
        end
        
        function obj = get_C( obj )
            % 已知Z0与Zd，计算匹配时Cs与Cp
            % 注意：不适用于复数Zs，仅适用于Zs=Z0∈R的情况
            Z0=obj.s.Z0;
            RD=obj.D.Z.R;
            XD=obj.D.Z.X;
            switch obj.m.t.model
                case ''
                    Xs=sqrt(RD*(Z0-RD))-XD;
                    Xp=-Z0*sqrt(RD/(Z0-RD));
                case 'ideal'
                    k=obj.m.t.n^2;
                    Xs=sqrt(k*RD*(Z0-k*RD))/k-XD;
                    Xp=-Z0*sqrt(k*RD/(Z0-k*RD));
                case {'mutual','mutual2'}
                    XL1=obj.w*obj.m.t.L1;
                    XL2=obj.w*obj.m.t.L2;
                    XM=obj.w*obj.m.t.M;
                    if ~(obj.m.R1>=0 && obj.m.R2>=0)
                        warning('R1 or  R2 < 0. Loss in the matching unit is ignored.')
                        obj.m.R1=0;
                        obj.m.R2=0;
                    end
                    R1=obj.m.R1;
                    R2=obj.m.R2;
                    if R1==0 && R2==0
                        % Lossless matching unit.
                        % the default solution of matching C
                        temp=sqrt(RD*(Z0*XM^2/XL1^2-RD));
                        if strcmp(obj.m.t.model,'mutual2')
                            % the second solution of matching C
                            % result in Cp<0
                            temp=-temp;
                        end
                        Xs = temp-(XD+XL2-XM^2/XL1);
                        Xp =-1/(1/XL1+temp/(RD*Z0));
                    else
                        % Lossy matching unit.
                        var1 = R1^2 + XL1^2 - R1*Z0;
                        var2 = (R2 + RD)*var1 + R1*XM^2;
                        % the default solution of matching C
                        temp = sqrt(-var2*(var2 - Z0*XM^2)/var1^2);
                        if strcmp(obj.m.t.model,'mutual2')
                            % the second solution of matching C
                            % result in Cp<0
                            temp=-temp;
                        end
                        Xs = temp-(XD+XL2-XM^2*XL1/var1);
                        Xp =-Z0*((R2+RD)*XL1*Z0-var1*temp)/...
                            (var2-Z0*XM^2+Z0*(R2+RD)*(Z0-R1));
                    end
                otherwise
                    error('No such type')
            end
            obj.m.s.Z.cm=1i*Xs;
            obj.m.s.Z=get_impedance( 'Z', obj.m.s.Z, obj.f );
            obj.m.p.Z.cm=1i*Xp;
            obj.m.p.Z=get_impedance( 'Z', obj.m.p.Z, obj.f );
        end
        
        function obj = get_all( obj, type )
            % 根据输入，计算Z与UIP
            switch type
                case {'PD','D2s_PD'}
                    obj = obj.get_all_D2s( 'PZ' );
                case 'D2s_Ps'
                    Ps=obj.s.P;
                    obj.D.P=Ps;
                    obj = obj.get_all_D2s( 'PZ' );
                    obj.D.P=Ps*obj.D.P/obj.s.P;
                    obj = obj.get_all_D2s( 'PZ' );
                    assert(abs((obj.s.P-Ps)/Ps)<1e-4)
                case 'IDrms'
                    obj.D.I.abs=sqrt(2)*obj.D.I.rms;
                    obj = obj.get_all_D2s( 'IZ' );
                case 'IDabs'
                    obj = obj.get_all_D2s( 'IZ' );
                case {'Ps','s2D_Ps'}
                    obj = obj.get_all_s2D( 'Z' );
                otherwise
                    error('No such type')
            end
        end
        
        function obj = get_all_D2s( obj, type )
            % 已知PD或ID、ZD与匹配电路模型，计算Zs与所有UIP
            obj.D=get_UIP( type, obj.D );
            obj.m.s.I.cm=obj.D.I.cm;
            obj.m.s.U.cm=obj.m.s.I.cm*obj.m.s.Z.cm;
            obj.m.s=get_UIP( 'UIcm', obj.m.s );
            Z_DCs=obj.D.Z.cm+obj.m.s.Z.cm;
            switch obj.m.t.model
                case ''
                    Z1=Z_DCs;
                    i1=obj.D.I.cm;
                case 'ideal'
                    Z1=(obj.m.t.n^2)*Z_DCs;
                    i1=obj.D.I.cm/obj.m.t.n;
                case 'mutual'
                    ZL1=1i*obj.w*obj.m.t.L1;
                    ZL2=1i*obj.w*obj.m.t.L2;
                    ZM=1i*obj.w*obj.m.t.M;
                    if ~(obj.m.R1>=0 && obj.m.R2>=0)
                        warning('R1 or  R2 < 0. Use lossless matching unit.')
                        obj.m.R1=0;
                        obj.m.R2=0;
                    end
                    if obj.m.R1==0 && obj.m.R2==0
                        % Lossless matching unit.
                        Z2=Z_DCs+ZL2;
                        Z1=ZL1-ZM^2/Z2;
                    else
                        % Lossy matching unit.
                        Z2=Z_DCs+ZL2+obj.m.R2;
                        Z1=obj.m.R1+ZL1-ZM^2/Z2;
                    end
                    i1=obj.D.I.cm*Z2/ZM;
                otherwise
                    error('No such type')
            end
            obj.m.p.U.cm=i1*Z1;
            obj.m.p.I.cm=obj.m.p.U.cm/obj.m.p.Z.cm;
            obj.m.p=get_UIP('UIcm',obj.m.p);
            obj.s.Z.cm=get_impedance('parallel',[Z1,obj.m.p.Z.cm]);
            obj.s.I.cm=i1+obj.m.p.I.cm;
            obj.s.U.cm=obj.s.I.cm*obj.s.Z.cm;
            obj.s=get_RF( 'UIcm', obj.s, obj.f );
            obj.m.P=obj.s.P-obj.D.P;
            obj.m.R=obj.m.P/obj.D.I.rms^2;
        end
        
        function obj = get_all_s2D( obj, type )
            % 已知Ps、Zs与匹配电路模型，计算ZD与所有UIP
            obj.s=get_RF( type, obj.s, obj.f, obj.s.P );
            obj.m.p.U.cm=obj.s.U.cm;
            obj.m.p.I.cm=obj.m.p.U.cm/obj.m.p.Z.cm;
            obj.m.p=get_UIP('UIcm',obj.m.p);
            i1=obj.s.I.cm-obj.m.p.I.cm;
            Z_sCp=-get_impedance('parallel',[obj.m.p.Z.cm,-obj.s.Z.cm]);
            switch obj.m.t.model
                case ''
                    Z2=Z_sCp;
                    obj.m.s.I.cm=i1;
                case 'ideal'
                    Z2=Z_sCp/(obj.m.t.n^2);
                    obj.m.s.I.cm=i1*obj.m.t.n;
                case 'mutual'
                    ZL1=1i*obj.w*obj.m.t.L1;
                    ZL2=1i*obj.w*obj.m.t.L2;
                    ZM=1i*obj.w*obj.m.t.M;
                    if ~(obj.m.R1>=0 && obj.m.R2>=0)
                        warning('R1 or  R2 < 0. Loss in the matching unit is ignored.')
                        obj.m.R1=0;
                        obj.m.R2=0;
                    end
                    if obj.m.R1==0 && obj.m.R2==0
                        % Lossless matching unit.
                        Z2=ZM^2/(ZL1-Z_sCp)-ZL2;
                        obj.m.s.I.cm=(i1*ZL1-obj.m.p.U.cm)/ZM;
                    else
                        % Lossy matching unit.
                        Z2=ZM^2/(obj.m.R1+ZL1-Z_sCp)-ZL2-obj.m.R2;
                        obj.m.s.I.cm=(i1*(obj.m.R1+ZL1)-obj.m.p.U.cm)/ZM;
                    end
                otherwise
                    error('No such type')
            end
            obj.m.s.U.cm=obj.m.s.I.cm*obj.m.s.Z.cm;
            obj.m.s=get_UIP( 'UIcm', obj.m.s );
            obj.D.Z.cm=Z2 - obj.m.s.Z.cm;
            obj.D.Z=get_impedance( 'Z', obj.D.Z, obj.f );
            obj.D.I.cm=obj.m.s.I.cm;
            obj.D.U.cm=obj.D.I.cm*obj.D.Z.cm;
            obj.D=get_UIP( 'UIcm', obj.D );
            obj.m.P=obj.s.P-obj.D.P;
            obj.m.R=obj.m.P/obj.D.I.rms^2;
        end
        
        function obj = get_PUI_multiD( obj, type )
            % 根据PD/ID计算 单/多驱 的PUI
            % type: 同get_UIP，如PZ/IZ
            obj.D=get_UIP( type, obj.D );
            if isfield(obj.D,'configuration') && obj.D.num>1
                % 多驱处理
                UpD=obj.D.U.cm-obj.D.I.cm*obj.D.leadline.Z.cm;
                for i=1:obj.D.num
                    switch obj.D.configuration
                        case 's'
                            % 串联多驱
                            obj.Dlist(i).I.cm=obj.D.I.cm;
                            obj.Dlist(i).U.cm=obj.D.I.cm*obj.Dlist(i).Z.cm;
                        case 'p'
                            % 并联多驱
                            obj.Dlist(i).U.cm=UpD;
                            obj.Dlist(i).I.cm=UpD/obj.Dlist(i).Z.cm;
                    end
                    Di=get_UIP( 'UIcm', obj.Dlist(i) );
                    % TODO: 自动遍历field
                    obj.Dlist(i).Z=Di.Z;
                    obj.Dlist(i).U=Di.U;
                    obj.Dlist(i).I=Di.I;
                    obj.Dlist(i).P=Di.P;
                    obj.Dlist(i).Q=Di.Q;
                end
                
                
            end
        end
        
    end
    
end