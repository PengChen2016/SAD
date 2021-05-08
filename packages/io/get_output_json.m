function [ str ] = get_output_json( in_obj, filename )
% 输出json格式字符串。若filename非空，则写json文件

in_obj=complexStruct2str(in_obj);
str=jsonencode(in_obj); % 结构体转为json格式文本
str = prettyjson(str); % json格式化显示

if ~isempty(filename)
    % 写JSON文件
    % 注意当前文件夹
    dir_path=pwd;
    if strcmp(dir_path(end-1:end),'io')
        fileID = fopen(['../../others/' filename '.json'],'w');
    else
        fileID = fopen(['./others/' filename '.json'],'w');
    end
    fprintf(fileID,'%s\n',str);
    fclose(fileID);
end
end

function in_obj=complexStruct2str(in_obj)
  % jsoncode不支持复数，因此转成字符数组
  % 递归调用子函数，以适用于多级结构体
 field_names = fieldnames(in_obj);
for i=1:length(field_names)
    field_name= field_names{i};
    field_value=in_obj.(field_name);
    if isa(field_value,'double') && ~isreal(field_value)
        in_obj.(field_name)=num2str(field_value);
    elseif isa(field_value,'struct')
        in_obj.(field_name)=complexStruct2str(field_value);
    end
end
end