function [ input ] = get_input_data( flag )
% 返回预置 模型输入数据
% If flag.electric_model='', only input_data of plasma model are got.

% TODO: 分离plasma与其他部分，以便于修改plasma参数

constants=get_constants();

%% input of  equivalent_medium_model_of_plasma
if isa(flag,'char')
    temp=flag;
    clear flag
    flag.input_plasma=temp;
end

if ~isfield(flag,'type_Xsec') || isempty(flag.type_Xsec)
    switch flag.input_plasma
        case '2011Chabert'
            flag.type_Xsec='e-Ar-Biagi';
            %             error('No given type_Xsec. ');
        otherwise
            flag.type_Xsec='e-H2-Phelps';
    end
    warning(['Use the default type_Xsec: ' flag.type_Xsec])
else
    fprintf('[INFO] Use type_Xsec: %s \n',flag.type_Xsec);
end

% 导入预置数据
fprintf('[INFO] Use input plasma dataset: %s \n',flag.input_plasma);
switch flag.input_plasma
    case {'CHARLIE_10Pa_4MHz_520W','CHARLIE_10Pa1MHz520W'}
        % 从一个 multi_coupled 数据中取出部分
        temp=get_input_plasma( '2019Raunera_CHARLIE_sweep' );
        switch flag.input_plasma
            case 'CHARLIE_10Pa_4MHz_520W'
                idx=(temp.p==10) & (temp.f==4e6) & (temp.Pin==520);
            case 'CHARLIE_10Pa1MHz520W'
                idx=(temp.p==10) & (temp.f==1e6) & (temp.Pin==520);
        end
        input.plasma.idx=idx; % 用于get_external
        input.plasma.p=temp.p(idx);
        input.plasma.f=temp.f(idx);
        input.plasma.Pin=temp.Pin(idx);
        input.plasma.ne=temp.ne(idx);
        input.plasma.Te=temp.Te(idx);
        input.plasma.Tg=temp.Tg(idx);
        input.plasma.w_RF=temp.w_RF(idx);
        input.plasma.ng=temp.ng(idx);
        input.plasma.r=inf;
        input.plasma.size=[length(input.plasma.p),...
            length(input.plasma.f),length(input.plasma.Pin)];
    case {'BUG_0.3Pa_1MHz_55kW'}
        % 从一个 multi_coupled 数据中取出部分
        temp=get_input_plasma( '2021Zielke_BUG_sweep' );
        switch flag.input_plasma
            case 'BUG_0.3Pa_1MHz_55kW'
                idx=(temp.p==0.3) & (temp.f==1e6) & (temp.Pin==55e3);
        end
        input.plasma.idx=idx; % 用于get_external
        input.plasma.p=temp.p(idx);
        input.plasma.f=temp.f(idx);
        input.plasma.Pin=temp.Pin(idx);
        input.plasma.ne=temp.ne(idx);
        input.plasma.Te=temp.Te(idx);
        input.plasma.Tg=temp.Tg;
        input.plasma.w_RF=temp.w_RF(idx);
        input.plasma.ng=temp.ng(idx);
        input.plasma.r=inf;
        input.plasma.size=[length(input.plasma.p),...
            length(input.plasma.f),length(input.plasma.Pin)];
    case 'CHARLIE_raza_sweep'
        input.plasma=get_input_plasma( '2019Raunera_CHARLIE_sweep' );
        ratio_origin2goal.ne_r=nonuniform_dist.get_ne_r_CHARLIE([0,45.5])/nonuniform_dist.get_ne_r_CHARLIE(10);
        input.plasma.ne= input.plasma.ne*ratio_origin2goal.ne_r;
    case 'CHARLIE_10Pa1MHz520W_nonuniform'
        ratio_origin2goal=nan(5,5);
        norm_ne_r10=nonuniform_dist.get_ne_r_CHARLIE(10); % origin data from experiments
        norm_ne_r0=nonuniform_dist.get_ne_r_CHARLIE(0); % peak value at the center
        for i=1:5
            dist_rp(i)=nonuniform_dist.get_nonuniform_dist_CHARLIE(['rp' num2str(i)]);
            for j=1:i
                ratio_origin2goal(i,j)=dist_rp(i).ne_r(j)/norm_ne_r10;
            end
        end
        ratio_origin2goal(1,5)=norm_ne_r0/norm_ne_r10;
        flag_temp.input_plasma='CHARLIE_10Pa1MHz520W';
        input=get_input_data( flag_temp );
        input.plasma.ne=input.plasma.ne*ratio_origin2goal;
        input.plasma.size=[1,1,1,15];
        input.plasma=rmfield(input.plasma,'idx');
        
        assert(input.plasma.ne(3,3)-6.21e16<1e14)
%     case 'CHARLIE_10Pa1MHz520W_nonuniform_r5z3'
%         ratio_origin2goal=nan(5,5);
%         norm_ne_r10=nonuniform_dist.get_ne_r_CHARLIE(10); % origin data from experiments
%         norm_ne_r0=nonuniform_dist.get_ne_r_CHARLIE(0); % peak value at the center
%         for i=1:5
%             dist_rp(i)=nonuniform_dist.get_nonuniform_dist_CHARLIE(['rp' num2str(i)]);
%             for j=1:i
%                 ratio_origin2goal(i,j)=dist_rp(i).ne_r(j)/norm_ne_r10;
%             end
%         end
%         ratio_origin2goal(1,5)=norm_ne_r0/norm_ne_r10;
%         flag_temp.input_plasma='CHARLIE_10Pa1MHz520W';
%         input=get_input_data( flag_temp );
%         input.plasma.ne=input.plasma.ne*ratio_origin2goal;
%         input.plasma.size=[1,1,1,15];
%         input.plasma=rmfield(input.plasma,'idx');
%         
%         assert(input.plasma.ne(3,3)-6.21e16<1e14)
    case 'CHARLIE_10Pa1MHz520W_sweepne'
        flag_temp.input_plasma='CHARLIE_10Pa1MHz520W';
        input=get_input_data( flag_temp );
        input.plasma.ne=logspace(16,18,21)';
        input.plasma.size=[1,1,1,21];
        input.plasma=rmfield(input.plasma,'idx');
    case 'CHARLIE_razcoil_sweep'
        input.plasma=get_input_plasma( '2019Raunera_CHARLIE_sweep' );
        ratio_origin2goal.ne_z=nonuniform_dist.get_ne_z_CHARLIE([-50,50])/nonuniform_dist.get_ne_z_CHARLIE([-200,200]);
        input.plasma.ne= input.plasma.ne*ratio_origin2goal.ne_z;
    case 'BUG_nonuniform_from_Tandem' % for the double helical coil
        % for 'BUG_z_nonuniform210716'
        %         z_part=([0 20 20+80/2 20+80 131]/131)';
        %         z_part=[z_part(1:end-1),z_part(2:end)];
        % for 'BUG_z_nonuniform220304'
        % 20220303 不使用区域离散值，而使用大量散点值
        z_part=(0:0.01:1)';       
        dist_zp=nonuniform_dist.get_ne_z_Tandem(z_part);
        norm_ne_z65=nonuniform_dist.get_ne_z_Tandem(65/131); % norm ne at z=65mm
        ratio_origin2goal=dist_zp/norm_ne_z65;

        flag_temp.input_plasma='BUG_0.3Pa_1MHz_55kW';
        flag_temp.type_Xsec='e-H2-Phelps';
        input=get_input_data( flag_temp ); % origin data from experiments % ne at z=65mm
        input.plasma.ne=input.plasma.ne*ratio_origin2goal; 
        
        input.plasma.size=[1,1,1,101];
        input.plasma=rmfield(input.plasma,'idx');
        
        assert(input.plasma.ne(1)<1e17)
    otherwise
        input.plasma=get_input_plasma( flag.input_plasma );
end
input.plasma.flag=flag.input_plasma;

%% 等离子体特征参数计算
input.plasma=get_characterstic_parameters( flag, input.plasma );

if isfield(flag,'electric_model') && ~isempty(flag.electric_model)
    %% input of electric model
    % 几何数据
    input.geometry=get_input_geometry( flag.input_geometry ); % 导入预置数据
    
    % 假设等离子体几何尺寸，与放电腔同中心的圆柱
    %简化：以放电腔半径为等离子体半径
    input.geometry.r_plasma_eff=input.geometry.r_chamber;
    input.geometry.l_plasma=input.geometry.l_chamber;
    % input.geometry.l_plasma=2*input.geometry.l_coil; %几何校正
    % 当Lmp≥Rp/veff，可见Rp对PER影响小，则lplasma对PER影响小
    input.plasma.r=input.geometry.r_plasma_eff;
    
    switch input.plasma.flag
        case 'CHARLIE_razcoil_sweep'
            input.geometry.l_plasma=input.geometry.l_coil;
    end
    
    % 外电路数据
    switch input.plasma.flag
        case {'CHARLIE_10Pa_4MHz_520W', 'CHARLIE_10Pa1MHz520W'}
            temp=get_input_external(flag, input.geometry, input.plasma.w_RF);
            input.external=temp;
            % Rcoil_th, Lcoil_th, Lcoil_ex等只有单个数据
            % 从多维数组中选择
            idx=input.plasma.idx;
            input.external.Rmetal_ex=temp.Rmetal_ex(idx);
            input.external.Rcoil_ex=temp.Rcoil_ex(idx);
            input.external.Icoil_rms=temp.Icoil_rms(idx);
            if strcmp(flag.Rmetal,'measured-Rmetal-woplasma')
                input.external.Rmetal=temp.Rmetal(idx);
            end
            input.external.Qcoil=temp.Qcoil(idx);
        otherwise
            input.external=get_input_external( flag, input.geometry, input.plasma.w_RF );
    end
    
    %% check flag for specified models
    flag_check=true;
    switch flag.electric_model
        case 'transformer-2018Jainb'
            disp('[INFO] The electric model transformer-2018Jainb needs:')
            flag_check=check_flag(flag.stoc_model,'stoc_model','2018Jainb-simplify');
            flag_check=check_flag(flag.skin_depth,'skin_depth','as-medium-simplified-finite-radius') && flag_check;
            flag_check=check_flag(flag.Rmetal,'Rmetal','measured-Rmetal-woplasma') && flag_check;
            flag_check=check_flag(flag.Lcoil,'Lcoil','calculated-Lcoil-woplasma') && flag_check;
            if flag_check
                disp('ok.')
            else
                warning('There are some differences.')
            end
        case 'analytical_base'
            disp('[INFO] The electric model analytical_base needs:')
            flag_check=check_flag(flag.stoc_model,'stoc_model','Vahedi-simplify');
            flag_check=check_flag(flag.skin_depth,'skin_depth','as-medium-simplified') && flag_check;
            flag_check=check_flag(flag.Rmetal,'Rmetal','calculated-Rcoil-woplasma') && flag_check;
            flag_check=check_flag(flag.Lcoil,'Lcoil','calculated-Lcoil-woplasma') && flag_check;
            if flag_check
                disp('ok.')
            else
                warning('There are some differences.')
            end
    end
    
    % get_output_json( input, 'input');
    
end
end


function flag_pass=check_flag(flag,name,excepted)
if ~strcmp(flag,excepted)
    disp([name ' : excepted ' excepted ', actual ' flag])
    flag_pass=false;
else
    flag_pass=true;
end
end