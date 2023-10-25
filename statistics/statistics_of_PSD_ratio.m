clear; close all;
stamp = 'plotted by statistics\_of\_PSD\_ratio.m';
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
mean_inside_lst = zeros(case_num,freq_num); % mean value of inside PSD(B_para)/PSD(B_trace) distribution
mean_outside_lst = zeros(case_num,freq_num); % mean value of outside PSD(B_para)/PSD(B_trace) distribution
%% analyses
if slt_freq == 0
    %% analyze case by case
    for i_case = 1 : case_num
        case_dir = case_dir_lst(i_case,:);
        data_dir = [work_dir,case_dir,'\'];
        %% import data
        freq_arr = importdata([data_dir,'freq.csv']);
        PSD_ratio_inside = importdata([data_dir,'PSD_ratio_inside.csv']);
        PSD_ratio_outside = importdata([data_dir,'PSD_ratio_outside.csv']);
        freq_slt = [0.01,0.05,0.10,0.15,0.20,0.25];
        mark_lst = ['a','b','c','d','f','g'];
        %% plot cc(|B|,Np); [require freq < 0.5 Hz]
        figure();
        FontSize = 10;
        displacement = 0;

        for i_freq = 1 : 6
            subplot(2,3,i_freq)
            index_sub = find(abs(freq_arr - freq_slt(i_freq))-min(abs(freq_arr - freq_slt(i_freq)))==0);
            freq_sub = freq_arr(abs(freq_arr - freq_slt(i_freq))-min(abs(freq_arr - freq_slt(i_freq)))==0);
            h61 = histogram(PSD_ratio_outside(index_sub,:),'BinLimits',[0,1],'NumBins',20,'Normalization','pdf','LineWidth',2,'FaceColor','none','EdgeColor','k');
            hold on
            h62 = histogram(PSD_ratio_inside(index_sub,:),'BinLimits',[0,1],'NumBins',20,'Normalization','pdf','LineWidth',2,'FaceColor','none','EdgeColor','m');
            h62.BinEdges = h62.BinEdges + displacement;
            hold on
            xline(0.5,':b','LineWidth',2)
            xlabel('$$\frac{{\rm PSD}({\rm B}_{\parallel})}{{\rm PSD}({\rm B}_{trace})}$','interpreter','latex'); ylabel('PDF');
            title(['(',mark_lst(i_freq),') Freq = ',num2str(roundn(freq_sub,-3)),' Hz'])
            set(gca,'FontSize',FontSize);
        end
        %     legend('Outside','Inside');
        sgtitle(case_dir,'FontSize',FontSize*2);
        text(0,-5,stamp);
        if save_or_not == 1
            saveas(gcf,[data_dir,'PSD_ratio.png']);
        end
        %% calculate mean value
        for i_freq = 1 : freq_num
            if i_freq <= length(freq_arr)
                freq_lst(i_case,i_freq) = freq_arr(i_freq);
                mean_inside_lst(i_case,i_freq) = mean(PSD_ratio_inside(i_freq,:));
                mean_outside_lst(i_case,i_freq) = mean(PSD_ratio_outside(i_freq,:));
            end
        end
    end
    mean_delta_lst = mean_inside_lst - mean_outside_lst;
    mean_delta_sign = sign(mean_delta_lst);
    if save_or_not == 1
        freq_arr_comp = importdata([work_dir,case_dir_lst(1,:),'\','freq.csv']);
        csvwrite([work_dir,'PSD_ratio_mean_delta_sign.csv'],[freq_arr_comp.';mean_delta_sign]);
    end
elseif slt_freq == 1
    freq_slt = 0.1;
    %% plot histogram of select frequency
    fig = figure();
    FontSize = 14;%7
    LineWidth = 2;
    m_color = [1,0,1];
    k_color = [0,0,0];
    set(fig,'defaultAxesColorOrder',[k_color;m_color]);
    displacement = 0;
    mark_lst = ['a','b','c','d','e','f','g','h','i','j','k','l'];
    for i_case = 1 : case_num
        case_dir = case_dir_lst(i_case,:);
        data_dir = [work_dir,case_dir,'\'];
        stat_dir = [work_dir,'statistics\'];
        %% import data
        freq_arr = importdata([data_dir,'freq.csv']);
        PSD_ratio_inside = importdata([data_dir,'PSD_ratio_inside.csv']);
        PSD_ratio_outside = importdata([data_dir,'PSD_ratio_outside.csv']);
        %% plot cc(|B|,Np); [require freq < 0.5 Hz]
        subplot(3,4,i_case);
        index_sub = find(abs(freq_arr - freq_slt) - min(abs(freq_arr - freq_slt)) == 0);
        freq_sub = freq_arr(abs(freq_arr - freq_slt) - min(abs(freq_arr - freq_slt)) == 0);
        yyaxis left
        h1 = histogram(PSD_ratio_outside(index_sub,:),'BinLimits',[0,1],'NumBins',20,'Normalization','pdf','LineWidth',1,'FaceColor','none','EdgeColor','k');
        ylabel('PDF(outside)');
        hold on
        yyaxis right
        h2 = histogram(PSD_ratio_inside(index_sub,:),'BinLimits',[0,1],'NumBins',20,'Normalization','pdf','LineWidth',1,'FaceColor','none','EdgeColor','m');
        h2.BinEdges = h2.BinEdges + displacement;
        ylabel('PDF(inside)');
        hold on
        xline(0.5,':b','LineWidth',1)
        xlabel('$${\rm PSD}({\rm B}_{\parallel})/{\rm PSD}({\rm B}_{trace})$','interpreter','latex');
        title(['(',mark_lst(i_case),') case ',num2str(i_case)])
        set(gca,'Xcolor','b','FontSize',FontSize,'LineWidth',LineWidth,'TickDir','out','XminorTick','on');
        %     legend('Outside','Inside');
    end
%     sgtitle(['PSD(B_{para})/PSD(B_{trace}) distribution at ',num2str(roundn(freq_sub,-2)),' Hz'],'FontSize',FontSize*2);
    sgtitle(['$${\rm PSD}({\rm B}_{\parallel})/{\rm PSD}({\rm B}_{trace})$ distribution at ',num2str(roundn(freq_sub,-2)),' Hz'],'FontSize',FontSize*1.5,'interpreter','latex')
    text(0,-0.15,stamp);
    if save_or_not == 1
        saveas(gcf,[stat_dir,'PSD_ratio_',num2str(roundn(freq_sub,-3)),'_Hz(norm).png']);
    end
end