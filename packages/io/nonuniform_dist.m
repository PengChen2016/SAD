% TODO：目前使用一种rz解耦线平均密度的平均方法，这与诊断是一致的
% 下一步使用一种rz耦合体积平均密度的计算方法，偏差更小

classdef nonuniform_dist
    properties
        in_type
        ne_z
        Te_z
        ne_r
        Te_r
        n_r
        n_z
    end
    
    methods
        % 不建议使用非 无输入构造函数
%         function obj=nonuniform_dist(in_type)
%             % 在CHARLIE分布基础上，乘以绝对值
%             dist=obj.get_nonuniform_dist_CHARLIE(in_type);
%             if isfield(dist,'ne_z')
%                 obj.ne_z=dist.ne_z;
%             end
%                 
%             obj.in_type=in_type;
%         end
    end
    
    methods(Static)
        %% nonuniform distribution data of CHARLIE
        % interface of get CHARLIE data
        function obj=get_nonuniform_dist_CHARLIE(in_type)
            % get ne/Te of spatially nonuniform distribution, accroding to the
            % distribution of CHARLIE(1MHz, 1Pa, 520W). Details are in
            % fit_CHARLIE_nonuniform_plasma.sfit
            
            % io
            % in_type: str
            % meaning of characters
            % '': output all data
            % first character: 'z' for axially, or 'r' for radially
            % then:
            % 'a': average the whole region
            % number b: position, at b mm
            % 'b~c': average region from b mm to c mm.
            % 'pn': equally divided into n part and then average
            
            % obj: struct
            % obj.ne/Te: norm ne/Te distribution related to the input type

            % operation according to input type
            % all
            if isempty(in_type)
                obj=nonuniform_dist.get_ref_CHARLIE();
                return
            end
            
            % point value=value_z*value_r
            % average value=volume integral/volume = mean_z*mean_rφ
            % mean_z=@(f_z,z1,z2) integral(f_z,z1,z2)/(z2-z1);
            % mean_rphi=@(f_r,r1,r2) 2*integral(f_r(r)*r,r1,r2)/(r2*r2-r1*r1);
            type = strsplit(in_type, '_');
            for in_str=type
                % average the whole region
                if strcmp(in_str{1}(2),'a')
                    switch in_str{1}(1)
                        case 'r'
                            r=[0,45.5];
                            obj.ne_r=nonuniform_dist.get_ne_r_CHARLIE(r);
                            obj.Te_r=nonuniform_dist.get_Te_r_CHARLIE(r);
                        case 'z'
                            z=[-200,200];
                            obj.ne_z=nonuniform_dist.get_ne_z_CHARLIE(z);
                            obj.Te_z=nonuniform_dist.get_Te_z_CHARLIE(z);
                        otherwise
                            error('No such type!')
                    end
                    continue
                end
                % equally divided into n part and then average
                if strcmp(in_str{1}(2),'p')
                    temp=in_str{1}(3:end);
                    n=str2double(temp(isstrprop(temp,'digit')));
                    assert(1<=n && n<=5)
                    switch in_str{1}(1)
                        case 'r'
                            obj.n_r=n;
                            r=0:45.5/n:45.5;
                            r=[r(1:end-1)', r(2:end)'];
                            obj.ne_r=nonuniform_dist.get_ne_r_CHARLIE(r);
                            obj.Te_r=nonuniform_dist.get_Te_r_CHARLIE(r);
                        case 'z'
                            obj.n_z=n;
                            z=-200:400/n:200;
                            z=[z(1:end-1)', z(2:end)'];
                            obj.ne_z=nonuniform_dist.get_ne_z_CHARLIE(z);
                            obj.Te_z=nonuniform_dist.get_Te_z_CHARLIE(z);
                        otherwise
                            error('No such type!')
                    end
                    continue
                end
                
                idx=strfind(in_str{1},'~');
                if isempty(idx)
                    % get value at point b
                    b=str2double(in_str{1}(2:end));
                    assert(~isnan(b))
                    switch in_str{1}(1)
                        case 'r'
                            obj.ne_r=nonuniform_dist.get_ne_r_CHARLIE(b);
                            obj.Te_r=nonuniform_dist.get_Te_r_CHARLIE(b);
                        case 'z'
                            obj.ne_z=nonuniform_dist.get_ne_z_CHARLIE(b);
                            obj.Te_z=nonuniform_dist.get_Te_z_CHARLIE(b);
                        otherwise
                            error('No such type!')
                    end
                    continue
                elseif length(idx)==1
                    % average region from b to c
                    b=str2double(in_str{1}(2:idx-1));
                    assert(~isnan(b))
                    c=str2double(in_str{1}(idx+1:end));
                    assert(c>b)
                    switch in_str{1}(1)
                        case 'r'
                            r=[b,c];
                            obj.ne_r=nonuniform_dist.get_ne_r_CHARLIE(r);
                            obj.Te_r=nonuniform_dist.get_Te_r_CHARLIE(r);
                        case 'z'
                            z=[b,c];
                            obj.ne_z=nonuniform_dist.get_ne_z_CHARLIE(z);
                            obj.Te_z=nonuniform_dist.get_Te_z_CHARLIE(z);
                        otherwise
                            error('No such type!')
                    end
                    continue
                else
                    error('Too much ~!')
                end
            end
        end
        
        function obj=get_ref_CHARLIE()
            % experiment data(point values) as reference
            % origin data
            z_ne=[200, 130, 120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0];
            ne_z=[0, 0.196224, 0.23341, 0.286613, 0.340961, 0.405034, 0.476545, 0.553776, 0.639588, 0.727689, 0.816362, 0.891304, 0.947941, 0.98913, 1.00229];
            z_Te=[130, 120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0];
            Te_z=[0.750324414, 0.773111028, 0.780272704, 0.801105955, 0.833982816, 0.866860659, 0.89778514, 0.930337587, 0.954425788, 0.968097363, 0.986001062, 0.988930615, 0.992184582, 1];
            r_ne=[45.5, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0];
            ne_r=[0, 0.035130017, 0.180627079, 0.393114122, 0.5617439, 0.699048337, 0.798758892, 0.865696707, 0.899862222, 0.939810461, 1];
            r_Te=[45, 40, 35, 30, 25, 20, 15, 10, 5, 0];
            Te_r=[0.973151812, 1, 0.971381704, 0.94632216, 0.91439163, 0.879948658, 0.854767708, 0.837248269, 0.83310566, 0.829764696];
            
            % transform
            x_transform=@(half) [-half half(end-1:-1:1)];
            y_transform=@(half) [half half(end-1:-1:1)];
            obj.ne_z=[x_transform(z_ne)', y_transform(ne_z)'];
            obj.Te_z=[x_transform(z_Te)', y_transform(Te_z)'];
            obj.ne_r=[flip(r_ne)', flip(ne_r)'];
            obj.Te_r=[flip(r_Te)', flip(Te_r)'];
        end
        
        % get ne/Te (r/z): get point or average values according to the
        % fitting expressions
        function ne_z=get_ne_z_CHARLIE(z)
            % z: N X 1/2 matrix
            assert(isempty(find(z>200 | z<-200,1)))
            [n1,n2]=size(z);
            assert(n2<3)
            f_ne_z=@(z) 1.003+0.0002894*z-0.0001769*z.^2+1.397e-06*z.^3-3.226e-09*z.^4;
            f_ne_z_n=@(z) f_ne_z(-z);
            ne_z=zeros(n1,1);
            if 1==n2
                % at point z
                idx1=z>=0;
                ne_z(idx1)=f_ne_z(z(idx1));
                idx2=z<0;
                ne_z(idx2)=f_ne_z_n(z(idx2));
            elseif 2==n2
                % average z1<z<z2
                for n=1:n1
                    z1=z(n,1);
                    z2=z(n,2);
                    assert(z1<z2)
                    if z1<0 && z2>0
                        ne_z_1=integral(f_ne_z_n,z1,0);
                        ne_z_2=integral(f_ne_z,0,z2);
                        ne_z(n)=(ne_z_1+ne_z_2)./(z2-z1);
                    elseif z1<=0 && z2<=0
                        ne_z(n)=integral(f_ne_z_n,z1,z2)./(z2-z1);
                    elseif z1>=0 && z2>=0
                        ne_z(n)=integral(f_ne_z,z1,z2)./(z2-z1);
                    else
                        error('Bug!')
                    end
                end
            end
        end
        function Te_z=get_Te_z_CHARLIE(z)
            assert(isempty(find(z>200 | z<-200,1)))
            [n1,n2]=size(z);
            assert(n2<3)
            f_Te_z=@(z) 672.1*Lorentzian([0,214.4],z);
            f_Te_z_n=@(z) f_Te_z(-z);
            Te_z=zeros(n1,1);
            if 1==n2
                % at point z
                idx1=z>=0;
                Te_z(idx1)=f_Te_z(z(idx1));
                idx2=z<0;
                Te_z(idx2)=f_Te_z_n(z(idx2));
            elseif 2==n2
                % average z1<z<z2
                for n=1:n1
                    z1=z(n,1);
                    z2=z(n,2);
                    assert(z1<z2)
                    if z1<0 && z2>0
                        Te_z_1=integral(f_Te_z_n,z1,0);
                        Te_z_2=integral(f_Te_z,0,z2);
                        Te_z(n)=(Te_z_1+Te_z_2)./(z2-z1);
                    elseif z1<=0 && z2<=0
                        Te_z(n)=integral(f_Te_z_n,z1,z2)./(z2-z1);
                    elseif z1>=0 && z2>=0
                        Te_z(n)=integral(f_Te_z,z1,z2)./(z2-z1);
                    else
                        error('Bug!')
                    end
                end
            end
        end
        function ne_r=get_ne_r_CHARLIE(r)
            assert(isempty(find(r>45.5 | r<0,1)))
            [n1,n2]=size(r);
            assert(n2<3)
            f_ne_r=@(r) 1.002-0.01921*r+0.001517*r.^2-6.698e-05*r.^3+7.112e-07*r.^4;
            f_integrand=@(r) r.*f_ne_r(r);
            ne_r=zeros(n1,1);
            if 1==n2
                % at point r
                ne_r=f_ne_r(r);
            elseif 2==n2
                % average r1<r<r2
                for n=1:n1
                    r1=r(n,1);
                    r2=r(n,2);
                    assert(r1<r2)
                    ne_r(n)=2*integral(f_integrand,r1,r2)/(r2*r2-r1*r1);
                end
            end
        end
        function Te_r=get_Te_r_CHARLIE(r)
            assert(isempty(find(r>45.5 | r<0,1)))
            [n1,n2]=size(r);
            f_Te_r=@(r) 0.8296+0.0007247*r-8.095e-05*r.^2+1.27e-05*r.^3-2.147e-07*r.^4;
            f_integrand=@(r) r.*f_Te_r(r);
            Te_r=zeros(n1,1);
            if 1==n2
                % at point r
                Te_r=f_Te_r(r);
            elseif 2==n2
                % average r1<r<r2
                for n=1:n1
                    r1=r(n,1);
                    r2=r(n,2);
                    assert(r1<r2)
                    Te_r(n)=2*integral(f_integrand,r1,r2)/(r2*r2-r1*r1);
                end
            end
        end
        
        %% nonuniform distribution data of Tandem
%         % interface of get Tandem data
%         function obj=get_nonuniform_dist_Tandem(in_type)
%             % get ne/Te of spatially nonuniform distribution, accroding to the
%             % distribution of Tandem RF ion source. Details are in
%             % fit_Tandem_nonuniform_plasma.sfit
%             
%             % io
%             % in_type: str
%             % meaning of characters
%             % '': output all data
%             % first character: 'z' for axially, or 'r' for radially
%             % then:
%             % 'a': average the whole region
%             % number b: position, at b mm
%             % 'b~c': average region from b mm to c mm.
%             % 'pn': equally divided into n part and then average
%             
%             % obj: struct
%             % obj.ne/Te: norm ne/Te distribution related to the input type
% 
%             % operation according to input type
%             % all
%             if isempty(in_type)
%                 obj=nonuniform_dist.get_ref_Tandem();
%                 return
%             end
%             
%             % point value=value_z*value_r
%             % average value=volume integral/volume = mean_z*mean_rφ
%             % mean_z=@(f_z,z1,z2) integral(f_z,z1,z2)/(z2-z1);
%             % mean_rphi=@(f_r,r1,r2) 2*integral(f_r(r)*r,r1,r2)/(r2*r2-r1*r1);
%             type = strsplit(in_type, '_');
%             for in_str=type
%                 % average the whole region
%                 if strcmp(in_str{1}(2),'a')
%                     switch in_str{1}(1)
%                         case 'r'
%                             r=[0,1];
%                             obj.ne_r=nonuniform_dist.get_ne_r_Tandem(r);
%                             obj.Te_r=nonuniform_dist.get_Te_r_Tandem(r);
%                         case 'z'
%                             z=[0,1];
%                             obj.ne_z=nonuniform_dist.get_ne_z_Tandem(z);
%                             obj.Te_z=nonuniform_dist.get_Te_z_Tandem(z);
%                         otherwise
%                             error('No such type!')
%                     end
%                     continue
%                 end
%                 % equally divided into n part and then average
%                 if strcmp(in_str{1}(2),'p')
%                     temp=in_str{1}(3:end);
%                     n=str2double(temp(isstrprop(temp,'digit')));
%                     assert(1<=n && n<=5)
%                     switch in_str{1}(1)
%                         case 'r'
%                             obj.n_r=n;
%                             r=0:45.5/n:45.5;
%                             r=[r(1:end-1)', r(2:end)'];
%                             obj.ne_r=nonuniform_dist.get_ne_r_Tandem(r);
%                             obj.Te_r=nonuniform_dist.get_Te_r_Tandem(r);
%                         case 'z'
%                             obj.n_z=n;
%                             z=-200:400/n:200;
%                             z=[z(1:end-1)', z(2:end)'];
%                             obj.ne_z=nonuniform_dist.get_ne_z_Tandem(z);
%                             obj.Te_z=nonuniform_dist.get_Te_z_Tandem(z);
%                         otherwise
%                             error('No such type!')
%                     end
%                     continue
%                 end
%                 
%                 idx=strfind(in_str{1},'~');
%                 if isempty(idx)
%                     % get value at point b
%                     b=str2double(in_str{1}(2:end));
%                     assert(~isnan(b))
%                     switch in_str{1}(1)
%                         case 'r'
%                             obj.ne_r=nonuniform_dist.get_ne_r_Tandem(b);
%                             obj.Te_r=nonuniform_dist.get_Te_r_Tandem(b);
%                         case 'z'
%                             obj.ne_z=nonuniform_dist.get_ne_z_Tandem(b);
%                             obj.Te_z=nonuniform_dist.get_Te_z_Tandem(b);
%                         otherwise
%                             error('No such type!')
%                     end
%                     continue
%                 elseif length(idx)==1
%                     % average region from b to c
%                     b=str2double(in_str{1}(2:idx-1));
%                     assert(~isnan(b))
%                     c=str2double(in_str{1}(idx+1:end));
%                     assert(c>b)
%                     switch in_str{1}(1)
%                         case 'r'
%                             r=[b,c];
%                             obj.ne_r=nonuniform_dist.get_ne_r_Tandem(r);
%                             obj.Te_r=nonuniform_dist.get_Te_r_Tandem(r);
%                         case 'z'
%                             z=[b,c];
%                             obj.ne_z=nonuniform_dist.get_ne_z_Tandem(z);
%                             obj.Te_z=nonuniform_dist.get_Te_z_Tandem(z);
%                         otherwise
%                             error('No such type!')
%                     end
%                     continue
%                 else
%                     error('Too much ~!')
%                 end
%             end
%         end
        
        % 20210712 以下使用数据的相关说明见 eP-210607-01双驱串并联融合\analyze 过程.pptx
        function obj=get_ref_Tandem()
            % experiment data(point values) as reference
            % origin data 2009Mcneely-fig11 use 55kW, 0.27Pa 
            z_ne=[4.94972, 6.93855, 8.94972, 10.9832, 12.9721, 14.9497, 17.0056, 18.9721, 20.9609];
            ne_z=[0.673149, 3.87638, 5.43701, 6.49519, 7.43421, 6.55381, 6.41646, 6.04456, 4.77779];
            % normalization
            z_ne=z_ne-4; % 校正过程见 eP-210607-01 双驱串并联融合\analyze 过程.pptx
            z_ne=z_ne/(27-4); 
            ne_z=ne_z/max(ne_z);
           % old code
%             z_Te=[130, 120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0];
%             Te_z=[0.750324414, 0.773111028, 0.780272704, 0.801105955, 0.833982816, 0.866860659, 0.89778514, 0.930337587, 0.954425788, 0.968097363, 0.986001062, 0.988930615, 0.992184582, 1];
%             r_ne=[45.5, 45, 40, 35, 30, 25, 20, 15, 10, 5, 0];
%             ne_r=[0, 0.035130017, 0.180627079, 0.393114122, 0.5617439, 0.699048337, 0.798758892, 0.865696707, 0.899862222, 0.939810461, 1];
%             r_Te=[45, 40, 35, 30, 25, 20, 15, 10, 5, 0];
%             Te_r=[0.973151812, 1, 0.971381704, 0.94632216, 0.91439163, 0.879948658, 0.854767708, 0.837248269, 0.83310566, 0.829764696];
            
            obj.ne_z=[z_ne', ne_z'];
%             obj.Te_z=[x_transform(z_Te)', y_transform(Te_z)'];
%             obj.ne_r=[flip(r_ne)', flip(ne_r)'];
%             obj.Te_r=[flip(r_Te)', flip(Te_r)'];
        end
        
        % get ne/Te (r/z): get point or average values according to the
        % fitting expressions
        function ne_z=get_ne_z_Tandem(z)
            % z: N X 1/2 matrix
            assert(isempty(find(z>1 | z<0,1)))
            [n1,n2]=size(z);
            assert(n2<3)
            f_ne_z=@(z) 0.01681+3.841*z-1.113*z.^2-6.674*z.^3 ...
            +1.007*z.^4-0.3662*z.^5+8.776*z.^6-5.347*z.^7;
            ne_z=zeros(n1,1);
            if 1==n2
                % at point z
                ne_z=f_ne_z(z);
            elseif 2==n2
                % average z1<z<z2
                for n=1:n1
                    z1=z(n,1);
                    z2=z(n,2);
                    assert(z1<z2)
                    ne_z(n)=integral(f_ne_z,z1,z2)./(z2-z1);
                end
            end
        end
        
        
    end
end