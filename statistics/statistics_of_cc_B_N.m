clear; close all;
stamp = 'plotted by statistics\_of\_cc\_B\_N.m';
save_or_not = 0;
slt_freq = 1;
%% select case
work_dir = 'E:\Research\Work\waves_in_exhaust_region\cases\';
case_num = 12;
case_dir_lst = ['20200529-1';'20200601-1';'20200601-2';'20200604-1';'20200608-1';'20200919-1';'20200920-1'; ...
    '20200925-1';'20210429-1';'20210429-2';'20210429-3';'20210429-4'];
%% select frequency
freq_num = 67;
freq_lst = zeros(case_num,freq_num);
mean_inside_lst = zeros(case_num,freq_num); % mean value of inside cc(|B|,Np) distribution
mean_outside_lst = zeros(case_num,freq_num); % mean value of outside cc(|B|,Np) distribution
cm_inside_lst = zeros(case_num,freq_num);
cm_outside_lst = zeros(case_num,freq_num);
%% analyses
if slt_freq == 0
    %% analyze case by case
    for i_case = 1 : case_num
        case_dir = case_dir_lst(i_case,:);
        data_dir = [work_dir,case_dir,'\'];
        %% import data
        freq_arr = importdata([data_dir,'freq.csv']);
        cc_inside = importdata([data_dir,'cc_inside.csv']);
        cc_outside = importdata([data_dir,'cc_outside.csv']);
        freq_slt = [0.01,0.05,0.10,0.15,0.20,0.25]; % 0.20 0.23
        mark_lst = ['a','b','c','d','f','g'];
        %% plot cc(|B|,Np); [require freq < 0.5 Hz]
        figure();
        FontSize = 12;
        displacement = 0;

        for i_freq = 1 : 6
            subplot(2,3,i_freq)
            index_sub = find(abs(freq_arr - freq_slt(i_freq))-min(abs(freq_arr - freq_slt(i_freq)))==0);
            freq_sub = freq_arr(abs(freq_arr - freq_slt(i_freq))-min(abs(freq_arr - freq_slt(i_freq)))==0);
            h61 = histogram(cc_outside(index_sub,:),'BinLimits',[-1,1],'NumBins',20,'Normalization','pdf','LineWidth',2,'FaceColor','none','EdgeColor','m');
            hold on
            h62 = histogram(cc_inside(index_sub,:),'BinLimits',[-1,1],'NumBins',20,'Normalization','pdf','LineWidth',2,'FaceColor','none','EdgeColor','k');
            h62.BinEdges = h62.BinEdges + displacement;
            hold on
            xline(0,':b','LineWidth',2)
            xlabel('CC(|B|,N_p)'); ylabel('PDF');
            title(['(',mark_lst(i_freq),') Freq = ',num2str(roundn(freq_sub,-3)),' Hz'])
            set(gca,'FontSize',FontSize);
        end
        %     legend('Outside','Inside');
        sgtitle(case_dir,'FontSize',FontSize*2);
        text(0,-2,stamp);
        if save_or_not == 1
            saveas(gcf,[data_dir,'cc_B_Np.png']);
        end
        %% calculate mean value
        for i_freq = 1 : freq_num
            if i_freq <= length(freq_arr)
                freq_lst(i_case,i_freq) = freq_arr(i_freq);
                mean_inside_lst(i_case,i_freq) = mean(cc_inside(i_freq,:));
                mean_outside_lst(i_case,i_freq) = mean(cc_outside(i_freq,:));
                cm_inside_lst(i_case,i_freq) = central_moment(cc_inside(i_freq,:),1,2);
                cm_outside_lst(i_case,i_freq) = central_moment(cc_outside(i_freq,:),1,2);
            end
        end
    end
    mean_delta_lst = mean_inside_lst - mean_outside_lst;
    cm_delta_lst = cm_inside_lst - cm_outside_lst;
    mean_delta_sign = sign(mean_delta_lst);
    cm_delta_sign = sign(cm_delta_lst);
    if save_or_not == 1
        freq_arr_comp = importdata([work_dir,case_dir_lst(1,:),'\','freq.csv']);
        csvwrite([work_dir,'cc_B_N_mean_delta_sign.csv'],[freq_arr_comp.';mean_delta_sign]);
        csvwrite([work_dir,'cc_B_N_cm_delta_sign.csv'],[freq_arr_comp.';cm_delta_sign]);
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
        freq_arr = importdata([data_dir,'freq.csv']);
        cc_inside = importdata([data_dir,'cc_inside.csv']);
        cc_outside = importdata([data_dir,'cc_outside.csv']);
        %% plot cc(|B|,Np); [require freq < 0.5 Hz]
        subplot(3,4,i_case)
        index_sub = find(abs(freq_arr - freq_slt)-min(abs(freq_arr - freq_slt))==0);
        freq_sub = freq_arr(abs(freq_arr - freq_slt)-min(abs(freq_arr - freq_slt))==0);
        h1 = histogram(cc_outside(index_sub,:),'BinLimits',[-1,1],'NumBins',20,'Normalization','pdf','LineWidth',1,'FaceColor','none','EdgeColor','m');
        hold on
        h2 = histogram(cc_inside(index_sub,:),'BinLimits',[-1,1],'NumBins',20,'Normalization','pdf','LineWidth',1,'FaceColor','none','EdgeColor','k');
        h2.BinEdges = h2.BinEdges + displacement;
        hold on
        xline(0,':b','LineWidth',1)
        xlabel('CC(|B|,N_p)'); ylabel('PDF');
        title(['(',mark_lst(i_case),') case ',num2str(i_case)])
        set(gca,'FontSize',FontSize);
        %     legend('Outside','Inside');
    end
    sgtitle(['CC(|B|,N_p) distribution at ',num2str(roundn(freq_sub,-3)),' Hz'],'FontSize',FontSize*2);
    text(0,-0.15,stamp);
    if save_or_not == 1
        saveas(gcf,[stat_dir,'cc_B_Np_',num2str(roundn(freq_sub,-3)),'_Hz.png']);
    end
end

%% function
function cm = central_moment(array,center,k)
    n = length(array);
    cm = sum((array - center).^k)/(n - 1);
end