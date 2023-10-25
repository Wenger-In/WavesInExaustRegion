clear; close all;
stamp = 'plotted by statistics\_of\_Vph\_div\_Va.m';
save_or_not = 0;
slt_freq = 0;
%% list of cases
num_case = 12;
case_lst = 1 : 1 : num_case;
case_enc_lst = [5,5,5,5,5,6,6, ...
    6,8,8,8,8];
case_day_lst = ['20200529';'20200601';'20200601';'20200604';'20200608';'20200919';'20200920'; ...
    '20200925';'20210429';'20210429';'20210429';'20210429'];
Vdata_lst = ['spc';'spi';'spi';'spc';'spi';'spc';'spc'; ...
    'spi';'spi';'spi';'spi';'spi'];
exht_beg_lst = ['22:26:00';'15:53:00';'19:51:00';'04:00:00';'11:05:30';'23:02:30';'20:53:40'; ...
    '13:33:10';'00:54:45';'08:14:30';'09:24:30';'13:40:40'];
exht_end_lst = ['22:50:00';'16:08:00';'19:52:30';'05:58:00';'11:06:30';'23:03:30';'20:54:10'; ...
    '13:34:00';'00:55:20';'08:28:00';'09:25:30';'13:41:40'];
back_beg_lst = ['22:00:00';'15:35:00';'19:53:00';'06:30:00';'11:04:00';'23:01:00';'20:52:00'; ...
    '13:31:10';'00:56:00';'08:00:00';'09:23:00';'13:39:00'];
back_end_lst = ['22:24:00';'15:50:00';'19:54:30';'08:30:00';'11:05:00';'23:02:00';'20:52:30'; ...
    '13:32:00';'00:56:35';'08:13:30';'09:24:00';'13:40:00'];
%% select case
work_dir = 'E:\Research\Work\waves_in_exhaust_region\cases\';
case_num = 12;
case_dir_lst = ['20200529-1';'20200601-1';'20200601-2';'20200604-1';'20200608-1';'20200919-1';'20200920-1'; ...
    '20200925-1';'20210429-1';'20210429-2';'20210429-3';'20210429-4'];
%% select frequency
freq_num = 16;
freq_lst = zeros(case_num,freq_num);
%% analyses
if slt_freq == 0
    %% analyze case by case
    for i_case = 1 : 1%case_num
        case_dir = case_dir_lst(i_case,:);
        data_dir = [work_dir,case_dir,'\'];
        %% import data
        period_arr = importdata([data_dir,'period.csv']);
        Vph_inside_arr = importdata([data_dir,'Vphase_inside.csv']);
        Vph_outside_arr = importdata([data_dir,'Vphase_outside.csv']);
        Vph_inside = Vph_inside_arr'; % size: 16*points
        Vph_outside = Vph_outside_arr'; % size: 16*points
        [nTime_inside,~] = size(Vph_inside_arr);
        [nTime_outside,~] = size(Vph_outside_arr);
        period = period_arr';
        freq_arr = 1./period; % size: 1*16
        freq_slt = [0.10,0.20,0.30,0.40,0.50,0.60]; % 0.20 0.23
        mark_lst = ['a','b','c','d','f','g'];
        %% import PSP data
        psp_dir = ['E:\Research\Data\PSP\Encounter ',num2str(case_enc_lst(i_case)),'\'];
        case_day = case_day_lst(i_case,:);
        exht_beg_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',exht_beg_lst(i_case,:)]);
        exht_end_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',exht_end_lst(i_case,:)]);
        back_beg_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',back_beg_lst(i_case,:)]);
        back_end_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',back_end_lst(i_case,:)]);

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
        Np_spc = spdfcdfread(spc_dir,'Variables','np_moment');

        % spi data
        spi_file = ['psp_swp_spi_sf00_l3_mom_inst_',case_day,'_v03.cdf'];
        spi_dir = [psp_dir,spi_file];
        spi_info = spdfcdfinfo(spi_dir);
        spi_epoch = spdfcdfread(spi_dir,'Variables','Epoch');
        Np_spi = spdfcdfread(spi_dir,'Variables','DENS');
        %% select spc or spi
        if strcmp(Vdata_lst(i_case,:),'spc')
            VEpoch = spc_epoch; Np = Np_spc;
        elseif strcmp(Vdata_lst(i_case,:),'spi')
            VEpoch = spi_epoch; Np = Np_spi;
        end
        %% emliminate bad points
        Brtn(abs(Brtn)>1e3)=nan;
        Np(Np<0)=nan;
        %% calculate for inside 
        std_epoch_inside = linspace(exht_beg_epoch - exht_beg_epoch, exht_end_epoch - exht_beg_epoch, nTime_inside);
        Brtn_interp_inside = zeros(nTime_inside,3);
        Np_interp_inside = zeros(nTime_inside,1);
        for i_col = 1:3
            Brtn_interp_inside(:,i_col) = interp1(fld_epoch - exht_beg_epoch,Brtn(:,i_col),std_epoch_inside,'pchip');
        end
        Np_interp_inside(:,1) = interp1(VEpoch - exht_beg_epoch,Np,std_epoch_inside,'pchip');
        Br_interp_inside = Brtn_interp_inside(:,1);
        Bt_interp_inside = Brtn_interp_inside(:,2);
        Bn_interp_inside = Brtn_interp_inside(:,3);
        B_mod_interp_inside = sqrt(Br_interp_inside.^2 + Bt_interp_inside.^2 + Bn_interp_inside.^2);
        % calculate Va
        mu0 = 4*pi*1e-7; mp = 1.67262192e-27; 
        Va_interp_inside = B_mod_interp_inside ./ sqrt(mu0 .* mp .* Np_interp_inside) * 1e-15; % unit: km/s
        Vph_div_Va_inside = Vph_inside_arr./Va_interp_inside;
        %% calculate for outside 
        std_epoch_outside = linspace(back_beg_epoch - back_beg_epoch, back_end_epoch - back_beg_epoch, nTime_outside);
        Brtn_interp_outside = zeros(nTime_outside,3);
        Np_interp_outside = zeros(nTime_outside,1);
        for i_col = 1:3
            Brtn_interp_outside(:,i_col) = interp1(fld_epoch - back_beg_epoch,Brtn(:,i_col),std_epoch_outside,'pchip');
        end
        Np_interp_outside(:,1) = interp1(VEpoch - back_beg_epoch,Np,std_epoch_outside,'pchip');
        Br_interp_outside = Brtn_interp_outside(:,1);
        Bt_interp_outside = Brtn_interp_outside(:,2);
        Bn_interp_outside = Brtn_interp_outside(:,3);
        B_mod_interp_outside = sqrt(Br_interp_outside.^2 + Bt_interp_outside.^2 + Bn_interp_outside.^2);
        % calculate Va
        mu0 = 4*pi*1e-7; mp = 1.67262192e-27; 
        Va_interp_outside = B_mod_interp_outside ./ sqrt(mu0 .* mp .* Np_interp_outside) * 1e-15; % unit: km/s
        Vph_div_Va_outside = Vph_outside_arr./Va_interp_outside;
        %% plot Vph_div_Va histogram
        figure();
        FontSize = 12;
        displacement = 0;

        for i_freq = 1 : 6
            subplot(2,3,i_freq)
            index_sub = find(abs(freq_arr - freq_slt(i_freq))-min(abs(freq_arr - freq_slt(i_freq)))==0);
            freq_sub = freq_arr(abs(freq_arr - freq_slt(i_freq))-min(abs(freq_arr - freq_slt(i_freq)))==0);
            h11 = histogram(Vph_div_Va_outside(index_sub,:),'BinLimits',[0,2],'NumBins',18,'Normalization','pdf','LineWidth',2,'FaceColor','none','EdgeColor','m');
            hold on
            h12 = histogram(Vph_div_Va_inside(index_sub,:),'BinLimits',[0,2],'NumBins',18,'Normalization','pdf','LineWidth',2,'FaceColor','none','EdgeColor','k');
            h12.BinEdges = h12.BinEdges + displacement; % take displacement to distinguish edges
            hold on
            xline(1,':b','LineWidth',2)
            xlabel('Vph/Va'); ylabel('PDF');
            title(['(',mark_lst(i_freq),') Freq = ',num2str(roundn(freq_sub,-3)),' Hz'])
            set(gca,'FontSize',FontSize);
        end

        %     legend('Outside','Inside');
        sgtitle(case_dir,'FontSize',FontSize*2);
        text(0,-0.005,stamp);
        if save_or_not == 1
            saveas(gcf,[data_dir,'Vph_div_Va.png']);
        end
        %% plot Vph_div_Va pcolor
        figure();
        subplot(2,1,1);
        p1 = pcolor(std_epoch_inside,freq_arr,Vph_div_Va_inside.');
        set(p1,'LineStyle','none');
        colorbar; colormap jet;
        set(gca,'Ycolor','k');
        subplot(2,1,2);
        p2 = pcolor(std_epoch_outside,freq_arr,Vph_div_Va_outside.');
        set(p2,'LineStyle','none');
        colorbar; colormap jet;
        set(gca,'Ycolor','m');
    end
elseif slt_freq == 1
    freq_slt = 0.1;
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
        Vph_inside_arr = importdata([data_dir,'Vphase_inside.csv']);
        Vph_outside_arr = importdata([data_dir,'Vphase_outside.csv']);
        Vph_inside = Vph_inside_arr'; % size: 16*points
        Vph_outside = Vph_outside_arr'; % size: 16*points
        [nTime_inside,~] = size(Vph_inside_arr);
        [nTime_outside,~] = size(Vph_outside_arr);
        period = period_arr';
        freq_arr = 1./period; % size: 1*16
        %% import PSP data
        psp_dir = ['E:\Research\Data\PSP\Encounter ',num2str(case_enc_lst(i_case)),'\'];
        case_day = case_day_lst(i_case,:);
        exht_beg_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',exht_beg_lst(i_case,:)]);
        exht_end_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',exht_end_lst(i_case,:)]);
        back_beg_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',back_beg_lst(i_case,:)]);
        back_end_epoch = datenum([case_day(1:4),'-',case_day(5:6),'-',case_day(7:8),' ',back_end_lst(i_case,:)]);

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
        Np_spc = spdfcdfread(spc_dir,'Variables','np_moment');

        % spi data
        spi_file = ['psp_swp_spi_sf00_l3_mom_inst_',case_day,'_v03.cdf'];
        spi_dir = [psp_dir,spi_file];
        spi_info = spdfcdfinfo(spi_dir);
        spi_epoch = spdfcdfread(spi_dir,'Variables','Epoch');
        Np_spi = spdfcdfread(spi_dir,'Variables','DENS');
        %% select spc or spi
        if strcmp(Vdata_lst(i_case,:),'spc')
            VEpoch = spc_epoch; Np = Np_spc;
        elseif strcmp(Vdata_lst(i_case,:),'spi')
            VEpoch = spi_epoch; Np = Np_spi;
        end
        %% emliminate bad points
        Brtn(abs(Brtn)>1e3)=nan;
        Np(Np<0)=nan;
        %% calculate for inside 
        std_epoch_inside = linspace(exht_beg_epoch - exht_beg_epoch, exht_end_epoch - exht_beg_epoch, nTime_inside);
        Brtn_interp_inside = zeros(nTime_inside,3);
        Np_interp_inside = zeros(nTime_inside,1);
        for i_col = 1:3
            Brtn_interp_inside(:,i_col) = interp1(fld_epoch - exht_beg_epoch,Brtn(:,i_col),std_epoch_inside,'pchip');
        end
        Np_interp_inside(:,1) = interp1(VEpoch - exht_beg_epoch,Np,std_epoch_inside,'pchip');
        Br_interp_inside = Brtn_interp_inside(:,1);
        Bt_interp_inside = Brtn_interp_inside(:,2);
        Bn_interp_inside = Brtn_interp_inside(:,3);
        B_mod_interp_inside = sqrt(Br_interp_inside.^2 + Bt_interp_inside.^2 + Bn_interp_inside.^2);
        % calculate Va
        mu0 = 4*pi*1e-7; mp = 1.67262192e-27; 
        Va_interp_inside = B_mod_interp_inside ./ sqrt(mu0 .* mp .* Np_interp_inside) * 1e-15; % unit: km/s
        Vph_div_Va_inside = Vph_inside_arr./Va_interp_inside;
        %% calculate for outside 
        std_epoch_outside = linspace(back_beg_epoch - back_beg_epoch, back_end_epoch - back_beg_epoch, nTime_outside);
        Brtn_interp_outside = zeros(nTime_outside,3);
        Np_interp_outside = zeros(nTime_outside,1);
        for i_col = 1:3
            Brtn_interp_outside(:,i_col) = interp1(fld_epoch - back_beg_epoch,Brtn(:,i_col),std_epoch_outside,'pchip');
        end
        Np_interp_outside(:,1) = interp1(VEpoch - back_beg_epoch,Np,std_epoch_outside,'pchip');
        Br_interp_outside = Brtn_interp_outside(:,1);
        Bt_interp_outside = Brtn_interp_outside(:,2);
        Bn_interp_outside = Brtn_interp_outside(:,3);
        B_mod_interp_outside = sqrt(Br_interp_outside.^2 + Bt_interp_outside.^2 + Bn_interp_outside.^2);
        % calculate Va
        mu0 = 4*pi*1e-7; mp = 1.67262192e-27; 
        Va_interp_outside = B_mod_interp_outside ./ sqrt(mu0 .* mp .* Np_interp_outside) * 1e-15; % unit: km/s
        Vph_div_Va_outside = Vph_outside_arr./Va_interp_outside;
        %% plot cc(|B|,Np); [require freq < 0.5 Hz]
        subplot(3,4,i_case)
        index_sub = find(abs(freq_arr - freq_slt)-min(abs(freq_arr - freq_slt))==0);
        freq_sub = freq_arr(abs(freq_arr - freq_slt)-min(abs(freq_arr - freq_slt))==0);
        h1 = histogram(Vph_div_Va_outside(index_sub,:),'BinLimits',[0,2],'NumBins',20,'Normalization','pdf','LineWidth',1,'FaceColor','none','EdgeColor','m');
        hold on
        h2 = histogram(Vph_div_Va_inside(index_sub,:),'BinLimits',[0,2],'NumBins',20,'Normalization','pdf','LineWidth',1,'FaceColor','none','EdgeColor','k');
        h2.BinEdges = h2.BinEdges + displacement;
        hold on
        xline(1,':b','LineWidth',1)
        xlabel('Vph/Va'); ylabel('PDF');
        title(['(',mark_lst(i_case),') case ',num2str(i_case)])
        set(gca,'FontSize',FontSize);
        %     legend('Outside','Inside');
    end
    sgtitle(['Vph/Va distribution at ',num2str(roundn(freq_sub,-3)),' Hz'],'FontSize',FontSize*2);
    text(0,-0.15,stamp);
    if save_or_not == 1
        saveas(gcf,[stat_dir,'Vph_div_Va_',num2str(roundn(freq_sub,-3)),'_Hz.png']);
    end
end