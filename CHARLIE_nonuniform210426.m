% CHARLIE_nonuniform case for paper
% 20210426 created by pengchen2016, Matlab R2017a

%% initialization
close all
clear
tic
addpath(genpath('./packages'))
now_str=datestr(now,'yyyymmdd_HHMMSS');
%% flag
solution_name='for_paper210415';
program_name='CHARLIE_sweepne210428';

addpath(genpath(['./others/' solution_name '/']))
flag.using_stored_data=false;
% flag.using_stored_data=true;
save_mat_name=['./others/' solution_name '/' program_name '.mat'];
%%%%%%%% plasma model
% flag.input_plasma='CHARLIE_10Pa1MHz520W_nonuniform';
flag.input_plasma='CHARLIE_10Pa1MHz520W_sweepne';

% stoc expression
% flag.stoc_model='';
% flag.stoc_model='Vahedi-simplify';
% flag.stoc_model='Cazzador-simplify';
% flag.stoc_model='2018Jainb-simplify';
flag.stoc_model='Cazzador-fit';

flag.medium_approximation='';
% flag.medium_approximation='consider_real(sigma)_only';
% flag.medium_approximation='sigma_dc';

flag.skin_depth='as-medium';
% flag.skin_depth='as-medium-simplified';
% flag.skin_depth='collisionless';
% flag.skin_depth='collision';
% flag.skin_depth='as-medium-simplified-finite-radius';

flag.output_plasma_model=true;
% flag.output_plasma_model=false;

%%%%%%%% electric model
flag.electric_model='';

%% solving
if ~flag.using_stored_data
    % calculate
    log_name=['./others/' solution_name '/' program_name '.log'];
    diary(log_name) % append to the end of the log file
    fprintf('\n-----%s %s-----\n\n',program_name,now_str)

    input=get_input_data( flag );
    input.plasma=plasma_model(flag, input.plasma);
    
    save(save_mat_name)
else
    % load data
    load(save_mat_name)
    disp(now_str)
    warning(['Using data stored in ' save_mat_name])
    if flag.output_plasma_model
        output_plasma_model(flag,input.plasma)
    end
end

%% post-processing
% check input
% fout=@(mat) [mat(1,1);mat(2,1:2)';mat(3,1:3)';mat(4,1:4)';mat(5,1:5)';mat(1,5)];
% 
% 
% % ref: test_get_CHARLIE_nonuniform_plasma / test_FEM_model_nonuniform_case
% norm_ne_r10=nonuniform_dist.get_ne_r(10); % origin data from experiments
% for i=1:5
%     dist_rp(i)=nonuniform_dist.get_nonuniform_dist_CHARLIE(['rp' num2str(i)]);
%     ratio_origin2goal(i).ne_r=dist_rp(i).ne_r/norm_ne_r10;
% end
% data_ref=nonuniform_dist.get_ref_CHARLIE();
% r_line=(0:45.5/100:45.5)';
% ne_r_fit=nonuniform_dist.get_ne_r(r_line);
% 
% figure
% scatter(data_ref.ne_r(:,1),data_ref.ne_r(:,2),'o','MarkerEdgeColor','k')
% hold on 
% plot(r_line,ne_r_fit,'-y')
% line_color={'r','g','b','c','m','y','k'};
% for i=1:5
%     r=0:45.5/i:45.5;
%     ne_r=[dist_rp(i).ne_r; 0];
%     for j=1:i
%         line([r(j),r(j+1)],[ne_r(j),ne_r(j)],'Color',line_color{i},'linestyle','-.')
%         line([r(j+1),r(j+1)],[ne_r(j),ne_r(j+1)],'Color',line_color{i},'linestyle','-.')
%     end
% end
% axis([0,45.5,0,1]) 
% xlabel('{\itr} [mm]')
% ylabel('Normalized \itn_{\rme}')
% grid on
% L1=legend('PIC/MCC','Fit','Uniform case');
% set(L1,'location','southwest');
% set(L1,'box','off')

% plot_1Y(input.plasma.ne, 'n_e [m^{-3}]');
% axis([0.3,10,5e16,5e17])
% plot_1Y(input.plasma.Te, 'T_e [eV]');
% axis([0.3,10,0,6])
% plot_1Y(input.external.Rmetal, 'R_{loss} [\Omega]');
% axis([0.3,10,0,0.3])
% plot_1Y(input.external.Icoil_rms, 'I_{rms} [A]');
% axis([0.3,10,0,60])
% 
% % check model
% % nu
% plot_2Yaxis(input.plasma.nu_m, '\nu_m [Hz]', input.plasma.nu_st, '\nu_{st} [Hz]');
% yyaxis left
% axis([0.3,10,2e6,2e8])
% line([0.3,10],[1e6*2*pi,1e6*2*pi],'linestyle','--')
% line([0.3,10],[4e6*2*pi,4e6*2*pi],'linestyle','--')
% yyaxis right
% axis([0.3,10,2e6,2e8])
% plot_1Y(input.plasma.nu_eff, '\nu_{eff} [Hz]');
% % skin_depth
% plot_2Yaxis(input.plasma.skin_depth, '\delta_{skin} [m]',input.plasma.wavelength, '\lambda [m]');
% yyaxis left
% line([0.3,8],[input.plasma.r,input.plasma.r],'linestyle','--')
% plot_2Yaxis(input.plasma.sigma_dc, '\sigma_{dc} [S/m]', input.plasma.skin_depth, '\delta_{skin} [m]');
% % power and impedance
% plot_2Yaxis(source.transformer.Rp, 'R_p [\Omega]',source.PER, 'R_{plasma} [\Omega]');
% plot_2Yaxis(source.Pplasma, 'P_{plasma} [W]',source.Psys, 'P_{sys} [W]');
% plot_2Yaxis(source.Rsys, 'R_{sys} [\Omega]',source.Xsys, 'X_{sys} [\Omega]');
% 
% plot_1Ylog(source.PER, 'R_{plasma} [\Omega]');
% axis([0.3,10,0.1,10])
% plot_1Y(source.PTE, '\eta');
% axis([0.3,10,0.5,1])
% % meidum
% plot_2Yaxis(input.plasma.sigma, '\sigma [S/m]', -input.plasma.w_RF.*input.plasma.eps_prime, '-\omega*\epsilon');
% yyaxis left
% axis([0.3,10,0,100])
% yyaxis right
% axis([0.3,10,0,100])
% plot_1Y(source.transformer.Lp, 'L');
% line([0.3,10],[source.transformer.Lmp,source.transformer.Lmp],'linestyle','--')
% legend('L_p, 1MHz','L_p, 4MHz','L_{mp}');
% % other
% plot_2Yaxis(input.plasma.wpe2wRF, '\omega_{pe}/\omega_{RF}',input.plasma.wpi2wRF, '\omega_{pi}/\omega_{RF}');

if ~flag.using_stored_data
    fprintf('\n-----END %s-----\n\n',now_str)
    diary off
end