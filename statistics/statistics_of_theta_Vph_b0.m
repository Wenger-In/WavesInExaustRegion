clear; close all;
stamp = 'plotted by statistics\_of\_theta\_Vph\_b0.m';
save_or_not = 0;
slt_freq = 0;
%% select case
work_dir = 'E:\Research\Work\waves_in_exhaust_region\cases\';
case_num = 12;
case_dir_lst = ['20200529-1';'20200601-1';'20200601-2';'20200604-1';'20200608-1';'20200919-1';'20200920-1'; ...
    '20200925-1';'20210429-1';'20210429-2';'20210429-3';'20210429-4'];
%% select frequency
freq_num = 16;
freq_lst = zeros(case_num,freq_num);
skew_inside_lst = zeros(case_num,freq_num); % skewness of inside theta_{VphB0} distribution
skew_outside_lst = zeros(case_num,freq_num);
kurt_inside_lst = zeros(case_num,freq_num); % kurtosis of inside theta_{VphB0} distribution
kurt_outside_lst = zeros(case_num,freq_num);
cm_inside_lst = zeros(case_num,freq_num); % central moment of inside theta_{VphB0} distribution
cm_outside_lst = zeros(case_num,freq_num);
%% analyses
if slt_freq == 0
    %% analyze case by case
    for i_case = 10 : 10%case_num
        case_dir = case_dir_lst(i_case,:);
        data_dir = [work_dir,case_dir,'\'];
        %% import data
        period_arr = importdata([data_dir,'period.csv']);
        theta_Vph_b0_inside_arr = importdata([data_dir,'theta_Vph_b0_inside.csv']);
        theta_Vph_b0_outside_arr = importdata([data_dir,'theta_Vph_b0_outside.csv']);
        theta_Vph_b0_inside = theta_Vph_b0_inside_arr'; % size: 16*points
        theta_Vph_b0_outside = theta_Vph_b0_outside_arr'; % size: 16*points
        period = period_arr';
        freq_arr = 1./period; % size: 1*16
        freq_slt = [0.10,0.20,0.30,0.40,0.50,0.60]; % 0.20 0.23
        mark_lst = ['a','b','c','d','e','f'];
        %% plot theta_Vph_B0
        figure();
        FontSize = 15;%12
        LinWidth = 2;
        displacement = 0;

        for i_freq = 1 : 6
            subplot(2,3,i_freq)
            index_sub = find(abs(freq_arr - freq_slt(i_freq))-min(abs(freq_arr - freq_slt(i_freq)))==0);
            freq_sub = freq_arr(abs(freq_arr - freq_slt(i_freq))-min(abs(freq_arr - freq_slt(i_freq)))==0);
            h11 = histogram(theta_Vph_b0_outside(index_sub,:),'BinLimits',[0,180],'NumBins',18,'Normalization','pdf','LineWidth',2,'FaceColor','none','EdgeColor','k');
            hold on
            h12 = histogram(theta_Vph_b0_inside(index_sub,:),'BinLimits',[0,180],'NumBins',18,'Normalization','pdf','LineWidth',2,'FaceColor','none','EdgeColor','m');
            h12.BinEdges = h12.BinEdges + displacement; % take displacement to distinguish edges
            hold on
            xline(90,':b','LineWidth',2)
            xlabel('\theta_{Vph,B0} [deg.]'); ylabel('PDF');
            title(['(',mark_lst(i_freq),') Freq = ',num2str(roundn(freq_sub,-3)),' Hz'])
            set(gca,'XminorTick','on','TickDir','out','LineWidth',LinWidth,'FontSize',FontSize);
        end

        legend('Outside','Inside');
        sgtitle(case_dir,'FontSize',FontSize*2);
        text(0,-0.005,stamp);
        if save_or_not == 1
            saveas(gcf,[data_dir,'theta_Vph_B0.png']);
        end
    end
    skew_delta_lst = abs(skew_inside_lst) - abs(skew_outside_lst); % consider first than kurt
    kurt_delta_lst = kurt_inside_lst - kurt_outside_lst;
    cm_delta_lst = cm_inside_lst - cm_outside_lst;
    cm_delta_sign = sign(cm_delta_lst);
    if save_or_not == 1
        csvwrite([work_dir,'theta_Vph_b0_cm_delta_sign.csv'],[freq_arr;cm_delta_sign]);
    end
elseif slt_freq == 1
    freq_slt = 0.15;
    %% plot histogram of select frequency
    figure();
    FontSize = 7;
    displacement = 0;
    mark_lst = ['a','b','c','d','e','f','g','h','i','j','k','l'];
    for i_case = 1 : case_num
        case_dir = case_dir_lst(i_case,:);
        data_dir = [work_dir,case_dir,'\'];
        stat_dir = [work_dir,'statistics\'];
        %% import data
        period_arr = importdata([data_dir,'period.csv']);
        theta_Vph_b0_inside_arr = importdata([data_dir,'theta_Vph_b0_inside.csv']);
        theta_Vph_b0_outside_arr = importdata([data_dir,'theta_Vph_b0_outside.csv']);
        theta_Vph_b0_inside = theta_Vph_b0_inside_arr'; % size: 16*points
        theta_Vph_b0_outside = theta_Vph_b0_outside_arr'; % size: 16*points
        period = period_arr';
        freq_arr = 1./period; % size: 1*16
        %% plot cc(|B|,Np); [require freq < 0.5 Hz]
        subplot(3,4,i_case)
        index_sub = find(abs(freq_arr - freq_slt)-min(abs(freq_arr - freq_slt))==0);
        freq_sub = freq_arr(abs(freq_arr - freq_slt)-min(abs(freq_arr - freq_slt))==0);
        h1 = histogram(theta_Vph_b0_outside(index_sub,:),'BinLimits',[0,180],'NumBins',20,'Normalization','pdf','LineWidth',1,'FaceColor','none','EdgeColor','k');
        hold on
        h2 = histogram(theta_Vph_b0_inside(index_sub,:),'BinLimits',[0,180],'NumBins',20,'Normalization','pdf','LineWidth',1,'FaceColor','none','EdgeColor','m');
        h2.BinEdges = h2.BinEdges + displacement;
        hold on
        xline(90,':b','LineWidth',1)
        xlabel('\theta_{Vph,B0} [deg.]'); ylabel('PDF');
        title(['(',mark_lst(i_case),') case ',num2str(i_case)])
        set(gca,'FontSize',FontSize);
        %     legend('Outside','Inside');
    end
    sgtitle(['\theta_{Vph,B0} distribution at ',num2str(roundn(freq_sub,-3)),' Hz'],'FontSize',FontSize*2);
    text(0,-0.15,stamp);
    if save_or_not == 1
        saveas(gcf,[stat_dir,'theta_Vph_b0_',num2str(roundn(freq_sub,-3)),'_Hz.png']);
    end
end

%% function
function cm = central_moment(array,center,k)
    n = length(array);
    cm = sum((array - center).^k)/(n - 1);
end