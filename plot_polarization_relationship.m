% This function is used to plot polaration relationship
% contain six panels 
% ***********************************************
%  sigma_m           theta_k_B0 
%  theta_db_B0       sigma_E' 
%  planarity         ellipse 
% ***********************************************
% for KAW
%  sigma_m > 1 or <1       theta_kB0 ~ 90
%  theta_db_B0 ~ 90        sigma_E' > 1 or <1
%  planarity ~ 1           ellipse ~ 1 
% ***********************************************
% for ICW
%  sigma_m > 1 or <1       theta_kB0 ~ 0
%  theta_db_B0 ~ 90        sigma_E' > 1 or <1
%  planarity ~ 1           ellipse ~ 1 
% ***********************************************
% Last change: 2021-03-08
% ***********************************************
clear; close all;
case_index = 1;
time_period = 0; % 0 show whole; 1 save inside; 2 save outside
save_or_not = 0;
%% list of cases
num_case = 12;
case_lst = 1 : 1 : num_case;
case_enc_lst = [5,5,5,5,5,6,6, ...
    6,8,8,8,8];
case_day_lst = ['20200529';'20200601';'20200604';'20200601';'20200608';'20200919';'20200920'; ...
    '20200925';'20210429';'20210429';'20210429';'20210429'];
plot_beg_lst = ['22:00:01';'15:35:01';'03:00:01';'19:48:01';'11:00:01';'23:00:01';'20:45:01'; ...
    '13:30:01';'00:50:01';'07:50:01';'09:20:01';'13:37:01'];
plot_end_lst = ['22:59:59';'16:24:59';'08:29:59';'19:56:59';'11:13:59';'23:05:59';'20:59:59'; ...
    '13:35:59';'00:59:59';'08:39:59';'09:29:59';'13:44:59'];
Vdata_lst = ['spc';'spi';'spc';'spi';'spi';'spc';'spc'; ...
    'spi';'spi';'spi';'spi';'spi'];
exht_beg_lst = ['22:26:00';'15:53:00';'04:00:00';'19:51:00';'11:05:30';'23:02:30';'20:53:40'; ...
    '13:33:10';'00:54:45';'08:14:30';'09:24:30';'13:40:40'];
exht_end_lst = ['22:50:00';'16:08:00';'05:58:00';'19:52:30';'11:06:30';'23:03:30';'20:54:10'; ...
    '13:34:00';'00:55:20';'08:28:00';'09:25:30';'13:41:40'];
back_beg_lst = ['22:00:00';'15:35:00';'06:30:00';'19:53:00';'11:04:00';'23:01:00';'20:52:00'; ...
    '13:31:10';'00:56:00';'08:00:00';'09:23:00';'13:39:00'];
back_end_lst = ['22:24:00';'15:50:00';'08:30:00';'19:54:30';'11:05:00';'23:02:00';'20:52:30'; ...
    '13:32:00';'00:56:35';'08:13:30';'09:24:00';'13:40:00'];
%% import data
psp_dir = ['E:\Research\Data\PSP\Encounter ',num2str(case_enc_lst(case_index)),'\'];
dTime = 0.5;
case_day = case_day_lst(case_index,:);
plot_beg_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',plot_beg_lst(case_index,:)]);
plot_end_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',plot_end_lst(case_index,:)]);
exht_beg_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',exht_beg_lst(case_index,:)]);
exht_end_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',exht_end_lst(case_index,:)]);
back_beg_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',back_beg_lst(case_index,:)]);
back_end_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',back_end_lst(case_index,:)]);

fld_list = ['00';'06';'12';'18'];
fld_file = ['psp_fld_l2_mag_rtn_',case_day,fld_list(1,:),'_v02.cdf'];
fld_dir = [psp_dir,fld_file];
fld_info = spdfcdfinfo(fld_dir);
% extract data
fld_epoch = spdfcdfread(fld_dir,'Variables','epoch_mag_RTN');
Brtn = spdfcdfread(fld_dir,'Variables','psp_fld_l2_mag_RTN');
for i_fld = 2 : length(fld_list)
    fld_file_sub = ['psp_fld_l2_mag_rtn_',case_day,fld_list(i_fld,:),'_v02.cdf'];
    fld_dir_sub = [psp_dir,fld_file_sub];
    fld_epoch_sub = spdfcdfread(fld_dir_sub,'Variables','epoch_mag_RTN');
    Brtn_sub = spdfcdfread(fld_dir_sub,'Variables','psp_fld_l2_mag_RTN');
    % joint data
    fld_epoch = cat(1,fld_epoch,fld_epoch_sub);
    Brtn = cat(1,Brtn,Brtn_sub);
end

% spc data
spc_file = ['psp_swp_spc_l3i_',case_day,'_v02.cdf'];
spc_dir = [psp_dir,spc_file];
spc_info = spdfcdfinfo(spc_dir);
spc_epoch = spdfcdfread(spc_dir,'Variables','Epoch');
Vrtn_spc = spdfcdfread(spc_dir,'Variables','vp_moment_RTN');
Np_spc = spdfcdfread(spc_dir,'Variables','np_moment');
sc_pos_HCI = spdfcdfread(spc_dir,'Variables','sc_pos_HCI');
sc_vel_HCI = spdfcdfread(spc_dir,'Variables','sc_vel_HCI');
sc_pos_HCIx = sc_pos_HCI(:,1); % [km]
sc_pos_HCIy = sc_pos_HCI(:,2); % [km]
sc_pos_HCIz = sc_pos_HCI(:,3); % [km]
sc_pos_dist = sqrt(sc_pos_HCIx.^2 + sc_pos_HCIy.^2 + sc_pos_HCIz.^2);
sc_vel_HCIx = sc_vel_HCI(:,1);
sc_vel_HCIy = sc_vel_HCI(:,2);
sc_vel_HCIz = sc_vel_HCI(:,3);
% switch sc_vel from HCI to scRTN (used for Vrtn_spi later)
[sc_vel_RTNr,sc_vel_RTNt,sc_vel_RTNn] = calc_HCI2SCRTN(sc_vel_HCIx,sc_vel_HCIy,sc_vel_HCIz,sc_pos_HCIx,sc_pos_HCIy,sc_pos_HCIz);
sc_vel_RTN = [sc_vel_RTNr,sc_vel_RTNt,sc_vel_RTNn];

% spi data
spi_file = ['psp_swp_spi_sf00_l3_mom_inst_',case_day,'_v03.cdf'];
spi_dir = [psp_dir,spi_file];
spi_info = spdfcdfinfo(spi_dir);
spi_epoch = spdfcdfread(spi_dir,'Variables','Epoch');
Vsc_spi = spdfcdfread(spi_dir,'Variables','VEL');
SC2RTN = spdfcdfread(spi_dir,'Variables','ROTMAT_SC_INST');
Np_spi = spdfcdfread(spi_dir,'Variables','DENS');
% switch Vsc_spi to Vrtn2sc_spi (spacecraft RTN frame)
Vrtn2sc_spi = SC2RTN * Vsc_spi.';
% switch Vrtn2sc_spi to Vrtn_spi (inertial RTN frame)
sc_vel_RTNr_spi = interp1(spc_epoch,sc_vel_RTNr,spi_epoch,'pchip');
sc_vel_RTNt_spi = interp1(spc_epoch,sc_vel_RTNt,spi_epoch,'pchip');
sc_vel_RTNn_spi = interp1(spc_epoch,sc_vel_RTNn,spi_epoch,'pchip');
Vr2sc_spi = -Vrtn2sc_spi(3,:); Vr_spi = Vr2sc_spi.' + sc_vel_RTNr_spi;
Vt2sc_spi =  Vrtn2sc_spi(1,:); Vt_spi = Vt2sc_spi.' + sc_vel_RTNt_spi;
Vn2sc_Spi = -Vrtn2sc_spi(2,:); Vn_spi = Vn2sc_Spi.' + sc_vel_RTNn_spi;
Vrtn_spi = [Vr_spi,Vt_spi,Vn_spi];
%% select spc or spi
if strcmp(Vdata_lst(case_index,:),'spc')
    VEpoch = spc_epoch; Np = Np_spc; Vrtn = Vrtn_spc;
elseif strcmp(Vdata_lst(case_index,:),'spi')
    VEpoch = spi_epoch; Np = Np_spi; Vrtn = Vrtn_spi;
end
%% extract time period
if time_period == 0
    data_beg_time = plot_beg_epoch; data_end_time = plot_end_epoch;
end
if time_period == 1
    data_beg_time = exht_beg_epoch; data_end_time = exht_end_epoch;
end
if time_period == 2
    data_beg_time = back_beg_epoch; data_end_time = back_end_epoch;
end
%% emliminate bad points
Brtn(abs(Brtn)>1e3)=nan;
% Vrtn(general_flag~=0,:)=nan;
Np(Np<0)=nan;
%% interp into standard time
% E_time = linspace(0,data_end_time - data_beg_time, floor((data_end_time - data_beg_time)*86400/0.034));
nTime = floor((data_end_time - data_beg_time)*86400/dTime);
std_time = linspace(data_beg_time - data_beg_time, data_end_time - data_beg_time, nTime);
Brtn_interp_vec = zeros(nTime,3);
Vrtn_interp_vec = zeros(nTime,3);

for i_col = 1:3
    Brtn_interp_vec(:,i_col) = interp1(fld_epoch - data_beg_time,Brtn(:,i_col),std_time,'pchip');
    Vrtn_interp_vec(:,i_col) = interp1(VEpoch - data_beg_time,Vrtn(:,i_col)',std_time,'pchip');
end
N_interp_vec = interp1(VEpoch - data_beg_time,Np,std_time,'pchip').';

% Et_interp = interp1(E_time,Et,std_time,'pchip');
% En_interp = interp1(E_time,En,std_time,'pchip');
Br_interp = Brtn_interp_vec(:,1).'; % size: 1*epoch
Bt_interp = Brtn_interp_vec(:,2).';
Bn_interp = Brtn_interp_vec(:,3).';
Vr_interp = Vrtn_interp_vec(:,1).'; % size: 1*epoch
Vt_interp = Vrtn_interp_vec(:,2).';
Vn_interp = Vrtn_interp_vec(:,3).';
% Br_interp = interp1(E_time,Br,std_time,'pchip');
% Bt_interp = interp1(E_time,Bt,std_time,'pchip');
% Bn_interp = interp1(E_time,Bn,std_time,'pchip');
fs=1./dTime; % sampling frequence;
%% calculate parameters
%     [wt_Et_interp,f1] = get_cwt(get_delta(Et_interp./1e3),fs,0.01,100);
%     [wt_En_interp, ~] = get_cwt(get_delta(En_interp./1e3),fs,0.01,100);
    [wt_Br,f1] = get_cwt(get_delta(Br_interp./1e9),fs,0.01,1);
    [wt_Bt,~] = get_cwt(get_delta(Bt_interp./1e9),fs,0.01,1);
    [wt_Bn,~] = get_cwt(get_delta(Bn_interp./1e9),fs,0.01,1);
%     [localBr,localBt,localBn] = calc_local_B(get_delta(Br_interp./1e9),get_delta(Bt_interp./1e9),get_delta(Bn_interp./1e9),f1,fs);
%     [local_Vr,local_Vt,local_Vn] = calc_local_V(Vrtn_interp_vec(:,1)'.*1e3,Vrtn_interp_vec(:,2)'.*1e3,Vrtn_interp_vec(:,3)'.*1e3,f1,fs);
%     [wt_delta_Et_prime,wt_delta_En_prime] = get_wt_E_prime(wt_Et_interp,wt_En_interp,local_Vr,local_Vt,local_Vn,wt_Br,wt_Bt,wt_Bn);
%% plot figures
freq_arr = get_wavelet_scale(Br_interp,fs); % unit:Hz [temp]
period_arr = 1./freq_arr; % unit:s [temp]

time_seq = std_time + data_beg_time;
linewidth = 2; fontsize = 12;
ylim_min = 0.02;% freq unit: Hz
ylim_max = 1;% freq unit: Hz
unit = '';
wt_delta_Et_prime = 0; wt_delta_En_prime = 0;
plot_figure(time_period,Br_interp,Bt_interp,Bn_interp,N_interp_vec,Vr_interp,Vt_interp,Vn_interp,wt_delta_Et_prime,wt_delta_En_prime,fs,time_seq, ...
    linewidth,fontsize,ylim_min,ylim_max,unit,case_day,exht_beg_epoch,exht_end_epoch,back_beg_epoch,back_end_epoch,save_or_not)
%% functions
function plot_figure(time_period,Br,Bt,Bn,Np,Vr,Vt,Vn,wt_delta_Et_prime,wt_delta_En_prime,fs,time_seq, ...
    linewidth,fontsize,ylim_min,ylim_max,unit,case_day,exht_beg,exht_end,back_beg,back_end,save_or_not)
% Br=ones(100,1)';
% Bt=sin(2*pi*2*[1:1:100]).^2;
% Bn=Br;
% fs=1/1;
% time_seq=[1:1:100];
% 
% linewidth = 1.2; fontsize = 14;
% ylim_min = 0.1;
% ylim_max = 1;
% *****************************

freq = get_wavelet_scale(Br,fs);
% ylim_min = min(freq); ylim_max = max(freq);
[localBr,localBt,localBn] = calc_local_B(Br,Bt,Bn,freq,fs); % size: freq*epoch
[num_freq,num_epoch] = size(localBr);

sigma_m = get_sigma_m(get_delta(Br), get_delta(Bt), get_delta(Bn), fs);
% sigma_E_prime = get_sigma_E_prime(wt_delta_Et_prime, wt_delta_En_prime);
[kr,kt,kn,dBr,dBt,dBn,planarity,ellipse,W_max,W_mid,W_min] = get_SVD(Br, Bt, Bn, fs);
% theta_db_B0_2D = get_theta_db_B0(dBr,dBt,dBn,localBr,localBt,localBn);
theta_k_B0 = get_theta_k_B0(kr,kt,kn,localBr,localBt,localBn);
%% calculate ratio of W_max to W_min
% fprintf('the ratio of W_max to W_mid is: %f \n',W_max./W_mid);
% fprintf('the ratio of W_mid to W_min is: %f \n',W_mid./W_min);
%% calculate LB-FAC coordinate
e_flow = zeros(3,num_epoch);
for i_epoch = 1 : num_epoch
    Vr_sub = Vr(i_epoch);
    Vt_sub = Vt(i_epoch);
    Vn_sub = Vn(i_epoch);
    e_flow(:,i_epoch) = [Vr_sub;Vt_sub;Vn_sub]/sqrt(Vr_sub.^2+Vt_sub.^2+Vn_sub.^2);
end
LBFAC = zeros(3,3,num_epoch,num_freq);
for i_epoch = 1 : num_epoch
    for i_freq = 1 : num_freq
        localBr_sub = localBr(i_freq,i_epoch);
        localBt_sub = localBt(i_freq,i_epoch);
        localBn_sub = localBn(i_freq,i_epoch);
        localB_magni_sub = sqrt(localBr_sub.^2+localBt_sub.^2+localBn_sub.^2);
        e_para_sub = [localBr_sub;localBt_sub;localBn_sub]/localB_magni_sub;
        e_perp1_sub = cross(e_para_sub,e_flow(:,i_epoch));
        e_perp2_sub = cross(e_para_sub,e_perp1_sub);
        LBFAC(:,1,i_epoch,i_freq) = e_para_sub;
        LBFAC(:,2,i_epoch,i_freq) = e_perp1_sub;
        LBFAC(:,3,i_epoch,i_freq) = e_perp2_sub;
    end
end
%% calculate para/perp1/perp2 component of B and V
B_LBFAC = zeros(3,num_epoch);
V_LBFAC = zeros(3,num_epoch);
for i_epoch = 1 : num_epoch
    for i_freq = 1 : num_freq
        LBFAC_sub = squeeze(LBFAC(:,:,i_epoch,i_freq));
        Brtn_sub = [Br(i_epoch),Bt(i_epoch),Bn(i_epoch)];
        B_LBFAC(:,i_epoch) = (Brtn_sub * LBFAC_sub).';
        Vrtn_sub = [Vr(i_epoch),Vt(i_epoch),Vn(i_epoch)];
        V_LBFAC(:,i_epoch) = (Vrtn_sub * LBFAC_sub).';
    end
end
B_para = B_LBFAC(1,:);
B_perp1 = B_LBFAC(2,:);
B_perp2 = B_LBFAC(3,:);
V_para = V_LBFAC(1,:);
V_perp1 = V_LBFAC(2,:);
V_perp2 = V_LBFAC(3,:);
%% calculate PSD
tlen = numel(Br);
fb = cwtfilterbank('SignalLength',tlen,'Wavelet','amor','SamplingFrequency',fs,'FrequencyLimits',[0.01 100]);
waver = wt(fb,Br);
wavet = wt(fb,Bt);
[waven,f] = wt(fb,Bn);
psdr = abs(waver).^2./f; psdt = abs(wavet).^2./f; psdn = abs(waven).^2./f;
psdtrace = psdr + psdt + psdn; % size: 90*36000

dBpara = waver.*localBr + wavet.*localBt + waven.*localBn;
dBpara = dBpara./sqrt(localBr.^2 + localBt.^2 + localBn.^2);
Bcompress = abs(dBpara).^2./f./psdtrace;
%% calculate cc(|B|,N)
B_magni = sqrt(Br.^2 + Bt.^2 + Bn.^2);
waveB_magni = wt(fb,B_magni);
waveNp = wt(fb,Np); 
% psdB = abs(waveB).^2./f;
% psdN = abs(waveN).^2./f;
cc_B_Np = calc_CC(waveB_magni,waveNp);
%% calculate cc(B_LBFAC_compoent,V_LBFAC_component)
waveB_para = wt(fb,B_para);
waveV_para = wt(fb,V_para);
cc_B_para_V_para = calc_CC(waveB_para,waveV_para);
waveB_perp1 = wt(fb,B_perp1);
waveV_perp1 = wt(fb,V_perp1);
cc_B_perp1_V_perp1 = calc_CC(waveB_perp1,waveV_perp1);
waveB_perp2 = wt(fb,B_perp2);
waveV_perp2 = wt(fb,V_perp2);
cc_B_perp2_V_perp2 = calc_CC(waveB_perp2,waveV_perp2);
%% calculate cc(N,V_para)
cc_Np_V_para = calc_CC(waveNp,waveV_para);
%% calculate V_{group}
% % period_bin = [0.10,0.13,0.17,0.22,0.28,0.37,0.48,0.62,0.81,1.05,1.36,1.76,2.29,2.97,3.85,5.00] % unit:s 
% % period_arr : [null,null,null,i= 2,i= 6,i=10,i=14,i=17,i=21,i=25,i=29,i=32,i=36,i=40,i=44,i=47]
% PoyntingFlux_r_file = importData('PoyntingFlux_r_arr.csv');
% PoyntingFlux_r_vec = PoyntingFlux_r_file.Data; % size: 36001*16
% PoyntingFlux_r = PoyntingFlux_r_vec';
% PoyntingFlux_r(:,36001) = []; % size: 16*36000
% mu0 = 1.25663706212e-6;
% 
% V_group_vec = zeros(13,36000);
% V_group_vec( 1,:) = PoyntingFlux_r( 4,:)./psdtrace( 2,:)*mu0*2e15; % unit: km/s
% V_group_vec( 2,:) = PoyntingFlux_r( 5,:)./psdtrace( 6,:)*mu0*2e15;
% V_group_vec( 3,:) = PoyntingFlux_r( 6,:)./psdtrace(10,:)*mu0*2e15;
% V_group_vec( 4,:) = PoyntingFlux_r( 7,:)./psdtrace(14,:)*mu0*2e15;
% V_group_vec( 5,:) = PoyntingFlux_r( 8,:)./psdtrace(17,:)*mu0*2e15;
% V_group_vec( 6,:) = PoyntingFlux_r( 9,:)./psdtrace(21,:)*mu0*2e15;
% V_group_vec( 7,:) = PoyntingFlux_r(10,:)./psdtrace(25,:)*mu0*2e15;
% V_group_vec( 8,:) = PoyntingFlux_r(11,:)./psdtrace(29,:)*mu0*2e15;
% V_group_vec( 9,:) = PoyntingFlux_r(12,:)./psdtrace(32,:)*mu0*2e15;
% V_group_vec(10,:) = PoyntingFlux_r(13,:)./psdtrace(36,:)*mu0*2e15;
% V_group_vec(11,:) = PoyntingFlux_r(14,:)./psdtrace(40,:)*mu0*2e15;
% V_group_vec(12,:) = PoyntingFlux_r(15,:)./psdtrace(44,:)*mu0*2e15;
% V_group_vec(13,:) = PoyntingFlux_r(16,:)./psdtrace(47,:)*mu0*2e15;
% period = [0.22,0.28,0.37,0.48,0.62,0.81,1.05,1.36,1.76,2.29,2.97,3.85,5.00];
%% plot figures
w_base = 0.12; h_base  = 0.05;
height = 0.14; width1 = 0.7; width2 = 0.7; space = 0.01;
%% plot sigma_m
% subplot('position',[w_base,h_base+height*2+space*2,width1,height]) % sigma_m
%     h = pcolor(time_seq,freq,sigma_m);
%     set(h,'LineStyle','none')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
%     colormap jet; colorbar;
% %     xlabel(['Time [' unit ']']);
%     yyaxis left
%     ylabel('Freq. [Hz]'); ylim([ylim_min ylim_max])
%     datetick('x','HH:MM')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     yyaxis right
%     ylabel('\sigma_m')
%     set(gca,'xticklabel',[],'yticklabel',[]); % 不显示y坐标轴刻度
%% plot cc(|B|,N)
subplot('position',[w_base,h_base+height*5+space*5,width1,height])
    h = pcolor(time_seq,freq,cc_B_Np);
    set(h,'LineStyle','none')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
    colormap jet; colorbar
    vertical_marks(exht_beg,exht_end,'-m');
    vertical_marks(back_beg,back_end,'-k');
%     xlabel(['Time [' unit ']']);
    yyaxis left
    ylabel('freq. [Hz]'); ylim([ylim_min ylim_max])
    yyaxis right
    ylabel('cc(|B|,Np)')
    datetick('x','HH:MM')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'xticklabel',[],'XMinorTick','on','ycolor','k'); % 不显示x坐标轴刻度
%% plot cc(B_para,V_para)
subplot('position',[w_base,h_base+height*4+space*4,width1,height])
    h = pcolor(time_seq,freq,cc_B_para_V_para);
    set(h,'LineStyle','none')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
    colormap jet; colorbar
    vertical_marks(exht_beg,exht_end,'-m');
    vertical_marks(back_beg,back_end,'-k');
%     xlabel(['Time [' unit ']']);
    yyaxis left
    ylabel('freq. [Hz]'); ylim([ylim_min ylim_max])
    yyaxis right
    ylabel('cc(B_{//},V_{//})')
    datetick('x','HH:MM')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'xticklabel',[],'XMinorTick','on','ycolor','k'); % 不显示x坐标轴刻度
%% plot cc(B_perp1,V_perp1)
subplot('position',[w_base,h_base+height*3+space*3,width1,height])
    h = pcolor(time_seq,freq,cc_B_perp1_V_perp1);
    set(h,'LineStyle','none')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
    colormap jet; colorbar
    vertical_marks(exht_beg,exht_end,'-m');
    vertical_marks(back_beg,back_end,'-k');
%     xlabel(['Time [' unit ']']);
    yyaxis left
    ylabel('freq. [Hz]'); ylim([ylim_min ylim_max])
    yyaxis right
    ylabel('cc(B_{\perp 1},V_{\perp1})')
    datetick('x','HH:MM')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'xticklabel',[],'XMinorTick','on','ycolor','k'); % 不显示x坐标轴刻度
%% plot cc(B_perp2,V_perp2)
subplot('position',[w_base,h_base+height*2+space*2,width1,height])
    h = pcolor(time_seq,freq,cc_B_perp2_V_perp2);
    set(h,'LineStyle','none')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
    colormap jet; colorbar
    vertical_marks(exht_beg,exht_end,'-m');
    vertical_marks(back_beg,back_end,'-k');
%     xlabel(['Time [' unit ']']);
    yyaxis left
    ylabel('freq. [Hz]'); ylim([ylim_min ylim_max])
    yyaxis right
    ylabel('cc(B_{\perp 2},V_{\perp 2})')
    datetick('x','HH:MM')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'xticklabel',[],'XMinorTick','on','ycolor','k'); % 不显示x坐标轴刻度
%% plot cc(N,V_para)
subplot('position',[w_base,h_base+height*1+space*1,width1,height])
    h = pcolor(time_seq,freq,cc_Np_V_para);
    set(h,'LineStyle','none')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
    colormap jet; colorbar
    vertical_marks(exht_beg,exht_end,'-m');
    vertical_marks(back_beg,back_end,'-k');
%     xlabel(['Time [' unit ']']);
    yyaxis left
    ylabel('freq. [Hz]'); ylim([ylim_min ylim_max])
    yyaxis right
    ylabel('cc(Np,V_{//})')
    datetick('x','HH:MM')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'xticklabel',[],'XMinorTick','on','ycolor','k'); % 不显示x坐标轴刻度
%% plot V_{group}
% subplot('position',[w_base,h_base+height*0+space*0,width2,height]) % V_group
%     h = pcolor(time_seq,period,V_group_vec);
%     set(h,'LineStyle','none')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
%     colormap jet; colorbar
% %     xlabel(['Time [' unit ']']);
%     ylabel('period. [s]'); ylim([min(period) max(period)])
%     hold on
%     plot(time_seq,ones(size(time_seq)),'w','linewidth',linewidth)
%     datetick('x','HH:MM')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     title('V_{group} [km/s]')
%% plot theta_{kB_0}
% subplot('position',[w_base,h_base+height*2+space*2,width1,height]) % theta_kB0
%     h = pcolor(time_seq,freq,theta_k_B0);
%     set(h,'LineStyle','none')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
%     colormap jet; colorbar
%     vertical_marks(exht_beg,exht_end,'-m');
%     vertical_marks(back_beg,back_end,'-k');
% %     xlabel(['Time [' unit ']']);
%     yyaxis left
%     ylabel('Freq. [Hz]'); ylim([ylim_min ylim_max])
%     yyaxis right
%     ylabel('$\theta_{\rm kB_0}$','interpreter','latex')
%     set(gca,'ytick',[],'yticklabel',[]); % 不显示y坐标轴刻度
%     grid off
%     datetick('x','HH:MM')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'xticklabel',[],'XMinorTick','on','ycolor','k'); % 不显示x坐标轴刻度
%% plot PSD(N)
% subplot('position',[w_base,h_base+height*1+space*1,width1,height]) % PSD(B_para)/PSD(B_trace)
%     h = pcolor(time_seq,freq,log(psdN));
%     set(h,'LineStyle','none')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
%     colormap jet; colorbar;
% %     xlabel(['Time [' unit ']']);
%     ylabel('Freq. [Hz]'); ylim([ylim_min ylim_max])
%     title('log[PSD(N)]')
%     datetick('x','HH:MM')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'xticklabel',[]); % 不显示x坐标轴刻度
%% plot PSD(B_{trace})
% subplot('position',[w_base,h_base+height*1+space*1,width1,height]) % PSD(B_para)/PSD(B_trace)
%     h = pcolor(time_seq,freq,log(psdtrace));
%     set(h,'LineStyle','none')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
%     colormap jet; colorbar;
%     vertical_marks(exht_beg,exht_end,'-m');
%     vertical_marks(back_beg,back_end,'-k');
% %     xlabel(['Time [' unit ']']);
%     yyaxis left
%     ylabel('Freq. [Hz]'); ylim([ylim_min ylim_max])
%     yyaxis right
%     ylabel('$log[{\rm PSD}({\rm B}_{trace})]$','interpreter','latex')
%     set(gca,'ytick',[],'yticklabel',[]); % 不显示y坐标轴刻度
%     grid off
%     datetick('x','HH:MM')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'xticklabel',[],'XMinorTick','on','ycolor','k'); % 不显示x坐标轴刻度
%% plot PSD(B_{para})/PSD(B_{trace}
subplot('position',[w_base,h_base+height*0+space*0,width2,height]) % PSD(B_para)/PSD(B_trace)
    h = pcolor(time_seq,freq,Bcompress);
    set(h,'LineStyle','none')
    xlim([time_seq(1),time_seq(length(time_seq))])
    set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
    colormap jet; colorbar;
    vertical_marks(exht_beg,exht_end,'-m');
    vertical_marks(back_beg,back_end,'-k');
%     xlabel(['Time [' unit ']']);
    yyaxis left
    ylabel('Freq. [Hz]'); ylim([ylim_min ylim_max])
    yyaxis right
    ylabel('$$\frac{{\rm PSD}({\rm B}_{\parallel})}{{\rm PSD}({\rm B}_{trace})}$','FontSize',fontsize/1.2,'interpreter','latex')
    set(gca,'ytick',[],'yticklabel',[],'XMinorTick','on','ycolor','k'); % 不显示y坐标轴刻度
    grid off
    datetick('x','HH:MM')
    xlim([time_seq(1),time_seq(length(time_seq))])
%% plot theta_{dbB_0}
% subplot(3,2,3) % theta_db_B0_2D
%     h = pcolor(time_seq,freq,theta_db_B0_2D);
%     set(h,'LineStyle','none')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
%     colormap jet; colorbar
%     xlabel(['Time [' unit ']']); ylabel('Freq. [Hz]'); ylim([ylim_min ylim_max])
%     title('\theta_{dbB_0}')
%     datetick('x','HH:MM')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%% plot sigma_{E''}
% subplot(3,2,4) % sigma_E'
%     h = pcolor(time_seq,freq,sigma_E_prime);
%     set(h,'LineStyle','none')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     xlabel(['Time [' unit ']']); ylabel('Freq. [Hz]'); ylim([ylim_min ylim_max])
%     title('\sigma_{E''}')
%     colormap jet; colorbar
%     set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
%     datetick('x','HH:MM')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%% plot planarity
% subplot(3,2,5) % planarity
%     h = pcolor(time_seq,freq,planarity);
%     set(h,'LineStyle','none')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
%     colormap jet; colorbar
%     xlabel(['Time [' unit ']']); ylabel('Freq. [Hz]'); ylim([ylim_min ylim_max])
%     title('planarity')
%     datetick('x','HH:MM')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%% plot ellipse
% subplot(3,2,6) % ellipse
%     h = pcolor(time_seq,freq,ellipse);
%     set(h,'LineStyle','none')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%     set(gca,'yscale','log','linewidth',linewidth,'fontsize',fontsize,'tickdir','out');
%     colormap jet; colorbar
%     xlabel(['Time [' unit ']']); ylabel('Freq. [Hz]'); ylim([ylim_min ylim_max])
%     title('ellipse')
%     datetick('x','HH:MM')
%     xlim([time_seq(1),time_seq(length(time_seq))])
%% sgtitle
% sgtitle(case_day);
%% save theta_{kB0}, PSD(B_para)/PSD(B_trace), cc(|B|,Np) and figure
if save_or_not == 1
    save_dir = 'E:\Research\Work\waves_in_exhaust_region\cases\';
    csvwrite([save_dir,'freq.csv'],freq);
    if time_period == 0
        saveas(gcf,['E:\Research\Work\waves_in_exhaust_region\cases\',case_day,'.png']);
    end
    if time_period == 1
        csvwrite([save_dir,'theta_k_B0_inside.csv'],theta_k_B0);
        csvwrite([save_dir,'PSD_ratio_inside.csv'],Bcompress);
        csvwrite([save_dir,'cc_inside.csv'],cc_B_Np);
    end
    if time_period == 2
        csvwrite([save_dir,'theta_k_B0_outside.csv'],theta_k_B0);
        csvwrite([save_dir,'PSD_ratio_outside.csv'],Bcompress);
        csvwrite([save_dir,'cc_outside.csv'],cc_B_Np);
    end
end
end

%% sub functions
function sigma_m = get_sigma_m(delta_Br, delta_Bt, delta_Bn, fs)
   [wt_delta_Br, ~] = get_cwt(double(delta_Br),fs,0.01,100);
   [wt_delta_Bt, ~] = get_cwt(double(delta_Bt),fs,0.01,100);
   [wt_delta_Bn, ~] = get_cwt(double(delta_Bn),fs,0.01,100);
    up = -2 .* imag(wt_delta_Bt .* conj(wt_delta_Bn));
    down = abs(wt_delta_Bt).^2 + abs(wt_delta_Bn).^2;
    sigma_m = up ./ down ;
end
function sigma_E = get_sigma_E_prime(wt_delta_Et_prime,wt_delta_En_prime)
    up = -2 .* imag(wt_delta_Et_prime .* conj(wt_delta_En_prime));
    down = abs(wt_delta_Et_prime).^2 + abs(wt_delta_En_prime).^2;
    sigma_E = up ./ down ;
end
function [kr,kt,kn,dBr,dBt,dBn,planarity,ellipse,W_max,W_mid,W_min] = get_SVD(Br, Bt, Bn, fs)
  % calculte wave vector k by SVD method
    [wt_Br, ~] = get_cwt(double(Br),fs,0.01,100);
    [wt_Bt, ~] = get_cwt(double(Bt),fs,0.01,100);
    [wt_Bn, ~] = get_cwt(double(Bn),fs,0.01,100);
    Bwt_size = size(wt_Br);
    for i = 1 : Bwt_size(1)
        for j = 1 : Bwt_size(2)
               S = get_S([wt_Br(i,j),wt_Bt(i,j), wt_Bn(i,j)]);
               A = [real(S);0,-imag(S(1,2)),-imag(S(1,3));imag(S(1,2)),0,-imag(S(2,3));imag(S(1,3)),imag(S(2,3)),0];
               [~,W,V] = svd(A);
               kr(i,j) = V(1,3);
               kt(i,j) = V(2,3);
               kn(i,j) = V(3,3);
               dBr(i,j) = V(1,1);
               dBt(i,j) = V(2,1);
               dBn(i,j) = V(3,1);
               W_min = W(3,3);
               W_mid = W(2,2);
               W_max = W(1,1);
               planarity(i,j) = 1 - sqrt(W_min./W_max);
               ellipse(i,j) = W_mid./W_max; 
        end
    end
end
function S = get_S(B)
    for i = 1 :3
        for j = 1:3
            S(i,j) = B(i) * conj(B(j));
        end
    end
end
function scale = get_wavelet_scale(Br,fs)
   [~, scale] = get_cwt(double(Br),fs,0.01,100);
end
function [localBx,localBy,localBz] = calc_local_B(Bx,By,Bz,scale,Fs)
% calculation the local mean magnetic field
%  Input: Bx,By,Bz(s,t),scale is the scale of the Gaussian window at each
%  scales
%  input scale is from the wavelet estimation procedure
% ref: Podesta, ApJ, 2009

time_num = numel(Bx);
freq_num = numel(scale);
Gaussian_scale = @(s,t0,t)exp(-(t0-t).^2./2./s.^2); 

localBx = zeros(size(Bx));
localBy = localBx;
localBz = localBx;

parfor i = 1:freq_num
    time_width = 3*scale(i);
    num_width = floor(time_width * Fs);
    time_arr = (-num_width:num_width) ./Fs;
    gaussin_window = Gaussian_scale(scale(i),time_arr,0);
    localBx(i,:) = conv(Bx,gaussin_window,'same');
    localBy(i,:) = conv(By,gaussin_window,'same');
    localBz(i,:) = conv(Bz,gaussin_window,'same');
end

end
function [localVx,localVy,localVz] = calc_local_V(Vx,Vy,Vz,scale,Fs)
% calculation the local mean magnetic field
%  Input: Bx,By,Bz(s,t),scale is the scale of the Gaussian window at each
%  scales
%  input scale is from the wavelet estimation procedure
% ref: Podesta, ApJ, 2009

time_num = numel(Vx);
freq_num = numel(scale);
Gaussian_scale = @(s,t0,t)exp(-(t0-t).^2./2./s.^2); 

localVx = zeros(size(Vx));
localVy = localVx;
localVz = localVx;

parfor i = 1:freq_num
    time_width = 3*scale(i);
    num_width = floor(time_width * Fs);
    time_arr = (-num_width:num_width) ./Fs;
    gaussin_window = Gaussian_scale(scale(i),time_arr,0);
    localVx(i,:) = conv(Vx,gaussin_window,'same');
    localVy(i,:) = conv(Vy,gaussin_window,'same');
    localVz(i,:) = conv(Vz,gaussin_window,'same');
end

end
function theta_db_B0_2D = get_theta_db_B0(dBr,dBt,dBn,localBr,localBt,localBn)
    wt_size = size(dBr);
    for i = 1:wt_size(1)
        for j = 1:wt_size(2)
            temp_db = [dBr(i,j), dBt(i,j), dBn(i,j)];
            temp_local_B0 = [localBr(i,j), localBt(i,j), localBn(i,j)];
%             temp_local_B0 = [nanmean(Br),nanmean(Bt),nanmean(Bn)];
            cos_theta_db_B0_2D(i,j) = dot(temp_db, temp_local_B0) ./ norm(temp_local_B0) ./ norm(temp_local_B0);        
        end
    end  
     theta_db_B0_2D = 90 - abs(acosd(cos_theta_db_B0_2D)-90);
end
function theta_k_B0 = get_theta_k_B0(kr,kt,kn,localBr,localBt,localBn)
    k_size = size(kr);
    for i = 1:k_size(1)
        for j = 1:k_size(2)
            temp_k = [kr(i,j), kt(i,j), kn(i,j)];
            temp_local_B0 = [localBr(i,j), localBt(i,j), localBn(i,j)];
            cos_theta_k_B0(i,j) = dot(temp_k, temp_local_B0)./ norm(temp_k)./norm(temp_local_B0);
        end
    end    
    theta_k_B0 = 90 - abs(acosd(cos_theta_k_B0)-90);
end
function delta_Data = get_delta(Data)
    delta_Data = Data - mean(Data,'omitnan');
end
function [wt_Et_prime,wt_En_prime] = get_wt_E_prime(wt_Et,wt_En,local_Vx,local_Vy,local_Vz,wt_Bx,wt_By,wt_Bz)
    size_wt_E = size(wt_Et);
    wt_Et_prime = wt_Et .*0;
    wt_En_prime = wt_En .*0;
    for i =1:size_wt_E(1)
        for j = 1:size_wt_E(2)
            local_V = [local_Vx(i,j),local_Vy(i,j),local_Vz(i,j)];
            local_B = [wt_Bx(i,j),wt_By(i,j),wt_Bz(i,j)];
            temp2 = local_V(3).*(local_B(1)) - local_V(1).*(local_B(3)) ;
            temp3 = local_V(1).*(local_B(2)) - local_V(2).*(local_B(1)) ;
            wt_Et_prime(i,j) = wt_Et(i,j) + temp2;
            wt_En_prime(i,j) = wt_En(i,j) + temp3;
        end
    end
end
function [result,f1] = get_cwt(Data,Fs,freq_min,freq_max)
    tlen = numel(Data);
    fb = cwtfilterbank('SignalLength',tlen,'Wavelet','amor','SamplingFrequency',Fs,'FrequencyLimits',[freq_min freq_max]);
    [result,f1] = wt(fb,Data);
end
function [CC_arr] = calc_CC(wave1,wave2)
    %This function used to calculate the coherency of two signals wave1 and
    %wave2 via CC = Re(WaWb*)/|Wa||Wb|
%   Detailed explanation goes here
    CC_arr = real(wave1.*conj(wave2))./abs(wave1)./abs(wave2);
end
function vertical_marks(mark_beg,mark_end,linestyle)
% Plot vertical lines
    LineWidth = 4;%2
    hold on;
    xline(mark_beg,linestyle,'LineWidth',LineWidth); hold on;
    xline(mark_end,linestyle,'LineWidth',LineWidth); hold on;
end
function [x_RTN,y_RTN,z_RTN] = calc_HCI2SCRTN(x_HCI,y_HCI,z_HCI,SC_HCIx,SC_HCIy,SC_HCIz)
% Change the coordiantes xyz in HCI frame to xyz in spacecraft RTN frame
%   input: x_HCI, y_HCI, z_HCI, the velocity in HCI frame (km/s)
%          SC_HCIx, SC_HCIy, SC_HCIz, the sapcecraft position in HCI frame (km)
%   output: x_RTN, y_RTN, z_RTN,the velocity in SC RTN frame (km/s)
% This function does not consider the move of the origin (using for velocity conversion)
    num = length(x_HCI);
    xyz_RTN = zeros(num,3);
    for i = 1:num
        Q = zeros(3,3);
        x1 = [1 0 0];
        y1 = [0 1 0];
        z1 = [0 0 1];
        x2 = [SC_HCIx(i),SC_HCIy(i),SC_HCIz(i)];
        if norm(x2)~= 0
            x2 = x2/norm(x2);
        end
        y2 = cross(z1,x2);
        if norm(y2)~= 0
            y2 = y2/norm(y2);
        end
        z2 = cross(x2,y2);
        if norm(z2)~= 0
            z2 = z2/norm(z2);
        end
        Q(1,1) = dot(x2,x1); Q(1,2) = dot(x2,y1); Q(1,3) = dot(x2,z1);
        Q(2,1) = dot(y2,x1); Q(2,2) = dot(y2,y1); Q(2,3) = dot(y2,z1);
        Q(3,1) = dot(z2,x1); Q(3,2) = dot(z2,y1); Q(3,3) = dot(z2,z1);
        xyz_RTN(i,:) = Q*[x_HCI(i);y_HCI(i);z_HCI(i)];
    end
    x_RTN = xyz_RTN(:,1); y_RTN = xyz_RTN(:,2); z_RTN = xyz_RTN(:,3);
end