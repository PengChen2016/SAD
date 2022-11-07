% 手动单个：table(runtests(''))

clc
clear
% addpath(genpath('../packages')) % 添加./packages及其子文件夹到路径

% 搜索所有待测试test_前缀的文件，添加到list
cd packages
list=ls('**/test_*.m');
cd ..
testfile_name_list=cellstr(list);
idx=zeros(1,1);
idx(1)=find(strcmp(testfile_name_list, 'test_all.m' ));
testfile_name_list(idx)=[];
num_testfile=length(testfile_name_list);
%测试
i=1;
results = runtests(testfile_name_list{i});
result_tables=table(results);
for i=2:num_testfile
    results = runtests(testfile_name_list{i}); %合并结果table变量
    result_tables=[result_tables;table(results)];
end

log_name='.\others\test_all.log';
clc
diary(log_name) % 重定向方式输出日志
disp(result_tables);
idx=find(~result_tables.Passed);
if isempty(idx)
    disp('All unit test passed.')
else
    num_failed=length(idx);
    for i=1:num_failed
        warning([result_tables.Name{idx(i)} ' failed!'])
        record=result_tables.Details{idx(i)}.DiagnosticRecord;
        disp(record.Report)
    end
end
close all
diary off
disp('See the details in ./others/test_all.log.')