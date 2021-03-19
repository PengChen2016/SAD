function [ str ] = get_output_json( in_obj, filename )
  % 输出json格式字符串，及写json文件
  % jsoncode不支持复数，因此转成字符数组
 field_names = fieldnames(in_obj);
for i=1:length(field_names)
    field_name= field_names{i};
    field_value=in_obj.(field_name);
    if isa(field_value,'double') && ~isreal(field_value)
        in_obj.(field_name)=num2str(field_value);
    end
end
str=jsonencode(in_obj); % 结构体转为json格式文本
str = prettyjson(str); % json格式化显示
% 写JSON文件
if ~isempty(filename)
    % 注意当前文件夹
    fileID = fopen(['./others/' filename '.json'],'w');
    fprintf(fileID,'%s\n',str);
    fclose(fileID);
end
end