clear; close all;
stamp = 'plotted by plot\_case\_statistics.m';
save_or_not = 0;
slt_or_not = 1;
%% import data
work_dir = 'E:\Research\Work\waves_in_exhaust_region\cases\statistics\';
qual_dir = 'theta_k_B0_mean';
data_file = [work_dir,qual_dir,'_delta_sign.csv'];
data_mat = importdata(data_file);
%% extract data
data_size = size(data_mat);

case_arr = 1:(data_size(1)-1);
case_slt_ind = [1:8,10:13];
case_slt = case_arr(case_slt_ind);

freq_arr = data_mat(1,:).';
freq_slt_ind = 21:37;
freq_slt = freq_arr(freq_slt_ind).';

qual_mat = data_mat(2:end,freq_slt_ind).';
qual_slt = qual_mat(:,case_slt_ind);
%% select plot
if slt_or_not == 0
    case_plot = 1 : length(case_arr);
    qual_plot = qual_mat;
elseif slt_or_not == 1
    case_plot = 1 : length(case_slt);
    qual_plot = qual_slt;
end
%% plot figure
figure();
FontSize = 20;
LineWidth = 2;

pc = pcolor(case_plot,freq_slt,qual_plot);
cb = colorbar;
set(cb,'Ticks',[-1,0,1]);
colormap gray;
xlim([1 length(case_plot)+1]);
xlabel('case index');
ylabel('freq [Hz]');
axis square
set(gca,'TickDir','out','XTick',[1:length(case_plot)],'YScale','log','LineWidth',LineWidth,'FontSize',FontSize);