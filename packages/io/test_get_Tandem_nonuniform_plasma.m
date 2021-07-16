%% Main function to generate tests
function tests = test_get_Tandem_nonuniform_plasma
% test get_Tandem_nonuniform_plasma
tests = functiontests(localfunctions);
end

%% Test Functions
function test_get_data_Tandem(testCase)
% choose nonuniform case for FEM model
% fit from minisource_YJH.mph
z_line=(0:1/100:1)';
ne_z_fit=nonuniform_dist.get_ne_z_Tandem(z_line);
% experiment data from 2009Mcneely-fig11 55kW, 0.27Pa 
data_ref=nonuniform_dist.get_ref_Tandem();
% experiment data from 2009Mcneely-fig11 34kW, 0.75Pa
z_ne=[4.93855, 6.96089, 8.94972, 10.9832, 12.9609, 18.9944, 21.0056];
ne_z=[0.689181, 4.93681, 15.9185, 20.8496, 20.854, 14.9128, 12.623];
z_ne=(z_ne-4)/(27-4); % normalization
ne_z=ne_z/max(ne_z);

% TODO: 与2021Zielke - RF power transfer efficiency and plasma parameters of low pressure high power ICPs
% 中的center-edge数据做对比

figure
plot(data_ref.ne_z(:,1),data_ref.ne_z(:,2),'-k')
hold on 
plot(z_line,ne_z_fit,'-.r')
plot(z_ne,ne_z,'--b')
plot(z_line,log10(1e3*ne_z_fit)/max(log10(1e3*ne_z_fit)),'-.m')
axis([0,1,0,1]) 
xlabel('Normalized {\itz}')
ylabel('Normalized {\itn}_e')
grid on
plot([20 20]/131,[0 1],'-k','LineWidth',1)
plot([20+80 20+80]/131,[0 1],'-k','LineWidth',1)
L1=legend('BUG实验-55kW, 0.27Pa','小源流体仿真-YJH','BUG实验-34kW, 0.75Pa','lg(ne) from 仿真');
set(L1,'location','south');
set(L1,'box','off')

z_part=([0 20 20+80/2 20+80 131]/131)';
z_part=[z_part(1:end-1),z_part(2:end)];
ne_z_part=nonuniform_dist.get_ne_z_Tandem(z_part);
figure
plot(z_line,ne_z_fit,'-k')
hold on 
for i=1:4
plot(z_part(i,:),[ne_z_part(i) ne_z_part(i)],'-.r')
end
ne_z_part(end+1)=nonuniform_dist.get_ne_z_Tandem(1);
for i=1:3
plot([z_part(i,2) z_part(i,2)],[ne_z_part(i) ne_z_part(i+1)],'-.r')
end
plot([20 20]/131,[0 1],'-k','LineWidth',1)
plot([20+80 20+80]/131,[0 1],'-k','LineWidth',1)
axis([0,1,0,1.01]) 
xlabel('Normalized {\itz}')
ylabel('Normalized {\itn}_e')
grid on
L1=legend('连续-拟合表达式','离散-区域取平均');
set(L1,'location','south');
set(L1,'box','off')
% TODO 待参考之前论文，绘制论文用图——在专门.m文件中
end

%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
% set a new path, for example
end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example
end

%% Optional fresh fixtures  
function setup(testCase)  % do not change function name
% open a figure, for example
end

function teardown(testCase)  % do not change function name
% close figure, for example
end
