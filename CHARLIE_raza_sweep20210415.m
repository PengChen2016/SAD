% CHARLIE_raza_sweep case for paper
% 20210415 created by pengchen2016, Matlab R2017a

%% initialization
close all
clear
tic
addpath(genpath('./packages'))
project_name='CHARLIE_raza_sweep210420';
now_str=datestr(now,'yyyymmdd_HHMMSS');
%% flag
flag.using_stored_data=false;
% flag.using_stored_data=true;
save_mat_name=['./others/' project_name '.mat'];
%%%%%%%% plasma model
flag.input_plasma='CHARLIE_raza_sweep';

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
% flag.electric_model='';
flag.electric_model='transformer-base';
% flag.electric_model='transformer-2011Chabert';
% flag.electric_model='transformer-2015Bandyopadhyay';
% flag.electric_model='transformer-2018Jainb';
% flag.electric_model='analytical_base';

flag.input_geometry='CHARLIE_base';

flag.Rmetal='measured-Rmetal-woplasma';
% flag.Rmetal=calculated-Rcoil-woplasma';
flag.Lcoil='measured-Lcoil-woplasma';
% flag.Lcoil='calculated-Lcoil-woplasma';

flag.output_electric_model=true;
% flag.output_electric_model=false;

%% solving
if ~flag.using_stored_data
    % calculate
    log_name=['.\others\' project_name '.log'];
    diary(log_name) % append to the end of the log file
    fprintf('\n-----%s-----\n\n',now_str)

    input=get_input_data( flag );
    
% modify geometry
% 应有flag.skin_depth='as-medium-simplified-finite-radius';
% 等离子体长度取为线圈长度，密度在线圈长度内平均
% geometry中的plasma尺寸需改变
% 区分source中的plasma尺寸（电模型等效尺寸），和plasma中的plasma尺寸（真实尺寸）
    
    input.plasma=plasma_model(flag, input.plasma);
    
    if isfield(flag,'electric_model') && ~isempty(flag.electric_model)
        source=electric_model( flag, input );
    end
    
    save(save_mat_name)
else
    % load data
    load(save_mat_name)
    disp(now_str)
    warning(['Using data stored in ' save_mat_name])
    if flag.output_plasma_model
        output_plasma_model(flag,input.plasma)
    end
    if flag.output_electric_model
        output_electric_model( flag, source )
    end
end

%% post-processing
% check input
plot_1Y(input.plasma.ne, 'n_e [m^{-3}]');
axis([0.3,10,5e16,5e17])
plot_1Y(input.plasma.Te, 'T_e [eV]');
axis([0.3,10,0,6])
plot_1Y(input.external.Rmetal, 'R_{loss} [\Omega]');
axis([0.3,10,0,0.3])
plot_1Y(input.external.Icoil_rms, 'I_{rms} [A]');
axis([0.3,10,0,60])

% check model
% nu
plot_2Yaxis(input.plasma.nu_m, '\nu_m [Hz]', input.plasma.nu_st, '\nu_{st} [Hz]');
yyaxis left
axis([0.3,10,2e6,2e8])
line([0.3,10],[1e6*2*pi,1e6*2*pi],'linestyle','--')
line([0.3,10],[4e6*2*pi,4e6*2*pi],'linestyle','--')
yyaxis right
axis([0.3,10,2e6,2e8])
plot_1Y(input.plasma.nu_eff, '\nu_{eff} [Hz]');
% skin_depth
plot_2Yaxis(input.plasma.skin_depth, '\delta_{skin} [m]',input.plasma.wavelength, '\lambda [m]');
yyaxis left
line([0.3,8],[input.plasma.r,input.plasma.r],'linestyle','--')
plot_2Yaxis(input.plasma.sigma_dc, '\sigma_{dc} [S/m]', input.plasma.skin_depth, '\delta_{skin} [m]');
% power and impedance
plot_2Yaxis(source.transform.Rp, 'R_p [\Omega]',source.PER, 'R_{plasma} [\Omega]');
plot_2Yaxis(source.Pplasma, 'P_{plasma} [W]',source.Psys, 'P_{sys} [W]');
plot_2Yaxis(source.Rsys, 'R_{sys} [\Omega]',source.Xsys, 'X_{sys} [\Omega]');

plot_1Ylog(source.PER, 'R_{plasma} [\Omega]');
axis([0.3,10,0.1,10])
plot_1Y(source.PTE, '\eta');
axis([0.3,10,0.5,1])
% meidum
plot_2Yaxis(input.plasma.sigma, '\sigma [S/m]', -input.plasma.w_RF.*input.plasma.eps_prime, '-\omega*\epsilon');
yyaxis left
axis([0.3,10,0,100])
yyaxis right
axis([0.3,10,0,100])
plot_1Y(source.transform.Lp, 'L');
line([0.3,10],[source.transform.Lmp,source.transform.Lmp],'linestyle','--')
legend('L_p, 1MHz','L_p, 4MHz','L_{mp}');
% other
plot_2Yaxis(input.plasma.wpe2wRF, '\omega_{pe}/\omega_{RF}',input.plasma.wpi2wRF, '\omega_{pi}/\omega_{RF}');

if ~flag.using_stored_data
    fprintf('\n-----END %s-----\n\n',now_str)
    diary off
end
%% aid function
% plot Y for specified X(size specified in this file: p×f)
function handle_fig=plot_1Y(Y, name_Y)
p=[0.3, 0.5, 1, 3, 5, 10]'; 
handle_fig=figure;
if max(Y(:))-min(Y(:))<100 && max(Y(:))/min(Y(:))<100
    semilogx(p,Y(:,1),'-.c');
    hold on
    semilogx(p,Y(:,2),'-.m');
else
    loglog(p,Y(:,1),'-.c');
    hold on
    loglog(p,Y(:,2),'-.m');
end
axis([0.3,10,-inf,inf])
xticks(p)
xlabel('{\itp} [Pa]');
ylabel(name_Y);
L1=legend('f=1MHz','f=4MHz');
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格
end

function handle_fig=plot_1Ylog(Y, name_Y)
p=[0.3, 0.5, 1, 3, 5, 10]'; 
handle_fig=figure;
loglog(p,Y(:,1),'-.c');
hold on
loglog(p,Y(:,2),'-.m');
axis([0.3,10,-inf,inf])
xticks(p)
xlabel('{\itp} [Pa]');
ylabel(name_Y);
L1=legend('f=1MHz','f=4MHz');
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格
end

function handle_fig=plot_2Yaxis(Y1, name_Y1, Y2, name_Y2)
p=[0.3, 0.5, 1, 3, 5, 10]'; 
handle_fig=figure;
if max(Y1(:))-min(Y1(:))<100 && max(Y1(:))/min(Y1(:))<100 ...
        && max(Y2(:))-min(Y2(:))<100 && max(Y2(:))/min(Y2(:))<100
    yyaxis left
    semilogx(p,Y1(:,1),'-.c');
    ylabel(name_Y1);
    axis([0.3,10,-inf,inf])
    yyaxis right
    semilogx(p,Y2(:,1),'--y');
    ylabel(name_Y2);
    axis([0.3,10,-inf,inf])
    hold on
    yyaxis left
    semilogx(p,Y1(:,2),'-.m');
    yyaxis right
    semilogx(p,Y2(:,2),'--g');
else
    yyaxis left
    loglog(p,Y1(:,1),'-.c');
    ylabel(name_Y1);
    axis([0.3,10,-inf,inf])
    yyaxis right
    loglog(p,Y2(:,1),'--y');
    ylabel(name_Y2);
    axis([0.3,10,-inf,inf])
    hold on
    yyaxis left
    loglog(p,Y1(:,2),'-.m');
    yyaxis right
    loglog(p,Y2(:,2),'--g');
end
xticks(p)
xlabel('{\itp} [Pa]');
L1=legend([name_Y1 ',f=1MHz'],[name_Y1 ',f=4MHz'],...
    [name_Y2 ',f=1MHz'],[name_Y2 ',f=4MHz']);
set(L1,'Location','best');
set(L1,'AutoUpdate','off');
grid on%显示网格
end