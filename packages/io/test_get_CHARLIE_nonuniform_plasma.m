%% Main function to generate tests
function tests = test_get_CHARLIE_nonuniform_plasma
% test get_CHARLIE_nonuniform_plasma
tests = functiontests(localfunctions);
end

%% Test Functions
function test_get_data_CHARLIE(testCase)
% test basic
% test get value at point
ref_dist=nonuniform_dist.get_ref_CHARLIE();
z=[130, 120, 110, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0]';
verifyEqual(testCase,nonuniform_dist.get_ne_z_CHARLIE([200;z]),ref_dist.ne_z(1:15,2),'Abs',0.01)
verifyEqual(testCase,nonuniform_dist.get_Te_z_CHARLIE(z),ref_dist.Te_z(1:14,2),'Abs',0.03)
r=flip([45, 40, 35, 30, 25, 20, 15, 10, 5, 0]');
verifyEqual(testCase,nonuniform_dist.get_ne_r_CHARLIE([r;45.5]),ref_dist.ne_r(1:11,2),'Abs',0.02)
verifyEqual(testCase,nonuniform_dist.get_Te_r_CHARLIE(r),ref_dist.Te_r(1:10,2),'Abs',0.01)
% test get average value
z=[0,200;-200,0;-200,200;-130,130];
actual_mean=nonuniform_dist.get_ne_z_CHARLIE(z);
excepted_mean=0.421*ones(3,1);
verifyEqual(testCase,actual_mean(1:3),excepted_mean,'Rel',0.04)
actual_mean=nonuniform_dist.get_Te_z_CHARLIE(z);
excepted_mean=10*trapz(ref_dist.Te_z(:,2))/260;
verifyEqual(testCase,actual_mean(4),excepted_mean,'Rel',0.01)
r=[0,45.5;0,45];
actual_mean=nonuniform_dist.get_ne_r_CHARLIE(r);
excepted_mean=2*5*trapz(ref_dist.ne_r(1:end-1,1).*ref_dist.ne_r(1:end-1,2))/(45*45);
verifyEqual(testCase,actual_mean(2),excepted_mean,'Rel',0.01)
actual_mean=nonuniform_dist.get_Te_r_CHARLIE(r);
excepted_mean=2*5*trapz(ref_dist.Te_r(1:end,1).*ref_dist.Te_r(1:end,2))/(45*45);
verifyEqual(testCase,actual_mean(2),excepted_mean,'Rel',0.01)
end
% test result: ok

function test_get_nonuniform_dist_CHARLIE(testCase)
% test get_nonuniform_dist_CHARLIE
dist0=nonuniform_dist.get_nonuniform_dist_CHARLIE('');
dist1=nonuniform_dist.get_nonuniform_dist_CHARLIE('r0_z0');
verifyEqual(testCase,dist1.ne_z,dist0.ne_z(dist0.ne_z(:,1)==0,2),'Rel',0.01)
verifyEqual(testCase,dist1.Te_r,dist0.Te_r(dist0.Te_r(:,1)==0,2),'Rel',0.01)

dist2=nonuniform_dist.get_nonuniform_dist_CHARLIE('ra_za');
dist3=nonuniform_dist.get_nonuniform_dist_CHARLIE('rp1_zp1');
verifyEqual(testCase,dist3.Te_z,dist2.Te_z,'Rel',0.01)
verifyEqual(testCase,dist3.ne_r,dist2.ne_r,'Rel',0.01)

dist4=nonuniform_dist.get_nonuniform_dist_CHARLIE('r0~9.1_z-200~-100');
dist5=nonuniform_dist.get_nonuniform_dist_CHARLIE('rp5_zp4');
verifyEqual(testCase,dist5.ne_r(1),dist4.ne_r,'Rel',0.01)
verifyEqual(testCase,dist5.ne_z(1),dist4.ne_z,'Rel',0.01)
end
% test result: ok

function test_electric_model_default_case(testCase)
% choose default case for electric models
dist_r0z0=nonuniform_dist.get_nonuniform_dist_CHARLIE('r0_z0')
dist_raza=nonuniform_dist.get_nonuniform_dist_CHARLIE('ra_za')
dist_origin=nonuniform_dist.get_nonuniform_dist_CHARLIE('r10_za') % origin data from experiments
ratio_origin2goal.ne_r=dist_raza.ne_r/dist_origin.ne_r;
ratio_origin2goal.Te_r=dist_raza.Te_r/dist_origin.Te_r;
disp('ratio_origin2goal for default case')
ratio_origin2goal
verifyEqual(testCase,nonuniform_dist.get_ne_r_CHARLIE([0,45.5])/nonuniform_dist.get_ne_r_CHARLIE(10),ratio_origin2goal.ne_r,'Rel',0.01)
end
% tset result: ne_raza=ne_origin*0.5344 as default case for electric models
% details are in ./other/Records of input and output data processing.xlsx

function test_FEM_model_nonuniform_case(testCase)
% choose nonuniform case for FEM model
norm_ne_r10=nonuniform_dist.get_ne_r_CHARLIE(10); % origin data from experiments
for i=1:5
    dist_rp(i)=nonuniform_dist.get_nonuniform_dist_CHARLIE(['rp' num2str(i)]);
    ratio_origin2goal(i).ne_r=dist_rp(i).ne_r/norm_ne_r10;
end
disp('ratio_origin2goal.ne_r for nonuniform case')
ratio_origin2goal.ne_r
data_ref=nonuniform_dist.get_ref_CHARLIE();
r_line=(0:45.5/100:45.5)';
ne_r_fit=nonuniform_dist.get_ne_r_CHARLIE(r_line);
Te_r_fit=nonuniform_dist.get_Te_r_CHARLIE(r_line);
% z的需要分150-100-150非均等三段
% z_line=(-200:400/100:200)';
% ne_z_fit=nonuniform_dist.get_ne_z_CHARLIE(z_line);
% Te_z_fit=nonuniform_dist.get_Te_z_CHARLIE(z_line);

figure
scatter(data_ref.ne_r(:,1),data_ref.ne_r(:,2),'o','MarkerEdgeColor','k')
hold on 
plot(r_line,ne_r_fit,'-y')
line_color={'r','g','b','c','m','y','k'};
for i=1:5
    r=0:45.5/i:45.5;
    ne_r=[dist_rp(i).ne_r; 0];
    for j=1:i
        line([r(j),r(j+1)],[ne_r(j),ne_r(j)],'Color',line_color{i},'linestyle','-.')
        line([r(j+1),r(j+1)],[ne_r(j),ne_r(j+1)],'Color',line_color{i},'linestyle','-.')
    end
end
axis([0,45.5,0,1]) 
xlabel('{\itr} [mm]')
ylabel('Normalized \itn_{\rme}')
grid on
L1=legend('PIC/MCC','Fit','Uniform case');
set(L1,'location','southwest');
set(L1,'box','off')

figure
scatter(data_ref.Te_r(:,1),data_ref.Te_r(:,2),'o','MarkerEdgeColor','k')
hold on 
plot(r_line,Te_r_fit,'-y')
line_color={'r','g','b','c','m','y','k'};
for i=1:5
    r=0:45.5/i:45.5;
    Te_r=[dist_rp(i).Te_r; 0];
    for j=1:i
        line([r(j),r(j+1)],[Te_r(j),Te_r(j)],'Color',line_color{i},'linestyle','-.')
        line([r(j+1),r(j+1)],[Te_r(j),Te_r(j+1)],'Color',line_color{i},'linestyle','-.')
    end
end
axis([0,45.5,0,1]) 
xlabel('{\itr} [mm]')
ylabel('Normalized \itT_{\rme}')
grid on
L1=legend('PIC/MCC','Fit','Uniform case');
set(L1,'location','southwest');
set(L1,'box','off')

% figure
% scatter(data_ref.ne_z(:,1),data_ref.ne_z(:,2),'o','MarkerEdgeColor','k')
% hold on 
% plot(z_line,ne_z_fit,'-y')
% line_color={'r','g','b','c','m','y','k'};
% for i=1:5
%     z=0:45.5/i:45.5;
%     ne_z=[dist_zp(i).ne_z; 0];
%     for j=1:i
%         line([z(j),z(j+1)],[ne_z(j),ne_z(j)],'Color',line_color{i},'linestyle','-.')
%         line([z(j+1),z(j+1)],[ne_z(j),ne_z(j+1)],'Color',line_color{i},'linestyle','-.')
%     end
% end
% axis([0,45.5,0,1]) 
% xlabel('{\itz} [mm]')
% ylabel('Normalized \itn_{\rme}')
% grid on
% L1=legend('PIC/MCC','Fit','Uniform case');
% set(L1,'location','southwest');
% set(L1,'box','off')

end
% test result: nonuniform case for FEM model
% details are in ./other/Records of input and output data processing.xlsx

%% aid function
% function line_xy=get_line_xy_p(dist_p_array)
% % get x/y to plot line, for dist of part-n array case
% len=length(dist_p_array);
% % ne_r
% for i=1:len
%     n=length(dist_p_array(i).ne_r);
%     
%     
%     r=0:45.5/n:45.5;
%     if n==1
%         line_xy(i).x=r;
%         line_xy(i).y=dist_p_array(i).ne_r*ones(1,2);
%     elseif n==2
%         
%     else
%         for j=1:n-1
%             line_xy(i)
%         end
% 
%         
%     end
% end
% end

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
