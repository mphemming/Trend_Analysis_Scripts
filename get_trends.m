function [trends] = get_trends(data_struct, data_struct_server, trend_data_struct,trend_data_struct_server,multiplier)

% get number of depths depending on input data
if iscell(trend_data_struct.EEMD_trend) ~= 1
    n_depths = size(trend_data_struct.EEMD_trend,1);
else
    n_depths = size(trend_data_struct.EEMD_trend,2);
end

% Loop to get trends
for n_depth = 1:n_depths
    
    % get significance boundaries
    sig = trend_data_struct_server.EEMD_conf_std_limit(1,:);
    sig_EAC = trend_data_struct_server.EEMD_conf_std_limit_EAC(1,:);
    tt_sig = data_struct.t_conv(n_depth).t;
    
    % if tt_sig time is not sam length get monthly binned version
    if length(tt_sig) ~= length(sig)
        % Need to interpolate to get dates
        min_t = nanmin(tt_sig); min_t =  datestr(min_t); min_t = datenum([num2str(14),min_t(3:end)]);
        max_t = nanmax(tt_sig); max_t =  datestr(max_t); max_t = datenum([num2str(14),max_t(3:end)]);
        [y_min,~,~,~,~,~] = datevec(min_t); [y_max,~,~,~,~,~] = datevec(max_t);
        n_index = 0;
        for n_yr = y_min:1:y_max
            for n_mn = 1:12
                n_index = n_index+1;
                time(n_index) = datenum(n_yr,n_mn,14,12,00,00);
            end
        end
        tt_sig = time;
    end
    
    % Get the trends 
    if iscell(trend_data_struct.EEMD_trend) ~= 1
        tr = trend_data_struct.EEMD_trend(n_depth,:) / nanstd(trend_data_struct.EEMD_T(n_depth,:)); % normalise trend 
        tr_EAC = trend_data_struct.EEMD_trend_EAC(n_depth,:)  / nanstd(trend_data_struct.EEMD_T(n_depth,:)); % normalise trend
        tr_0 = tr-tr(1);  
        tr_EAC_0 = tr_EAC-tr_EAC(1);  
        a = squeeze(trend_data_struct.EEMD_t(n_depth,:,:));
    else
        tr = trend_data_struct.EEMD_trend{n_depth} / nanstd(trend_data_struct.EEMD_T{n_depth}); % normalise trend 
        tr_EAC = trend_data_struct.EEMD_trend_EAC{n_depth}  / nanstd(trend_data_struct.EEMD_T{n_depth}); % normalise trend
        tr_0 = tr-tr(1);  
        tr_EAC_0 = tr_EAC-tr_EAC(1);  
        a = squeeze(trend_data_struct.EEMD_t{n_depth});
    end
    clear aa
    for nt = 1:size(a,1)
        aa(nt) = convertCharsToStrings(a(nt,:));
    end
    tt = datenum(aa);

    % if significance was done monthly, interpolate onto daily grid
    if length(tt) ~= length(sig)
        sig = interp1(tt_sig,sig,tt,'Linear');
        sig_EAC = interp1(tt_sig,sig_EAC,tt,'Linear');
    end

    % get significance from start of time period of interest
    sig = sig(1:numel(tt));
    sig_EAC = sig_EAC(1:numel(tt));
    % flag as NaN if insignificant trend value
    tr_0_sig = ones(size(tr_0));
    tr_EAC_0_sig = ones(size(tr_EAC_0));
    clear tr_0_sig tr_EAC_0_sig
    for n = 1:numel(tt)
        if tr_0(n) <= sig(n)
            tr_0_sig(n) = 0;
        else
            tr_0_sig(n) = 1;
        end
        if tr_EAC_0(n) <= sig_EAC(n)
            tr_EAC_0_sig(n) = 0;
        else
           tr_EAC_0_sig(n) = 1; 
        end
    end
    % 1950s
    check = tt >= datenum(1950,01,01) & tt < datenum(1960,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t1950s_sig = nanmean(T_rate_sig)*multiplier;   
        trends(n_depth).t1950s_sig_std = nanstd(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t1950s_sig = NaN;   
        trends(n_depth).t1950s_sig_std = NaN;   
    end
    trends(n_depth).t1950s_insig = nanmean(T_rate_insig)*multiplier;    
    trends(n_depth).t1950s_insig_std = nanstd(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t1950s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
        trends(n_depth).t1950s_EAC_sig_std = nanstd(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t1950s_EAC_sig = NaN;   
        trends(n_depth).t1950s_EAC_sig_std = NaN;   
    end
    trends(n_depth).t1950s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier; 
    trends(n_depth).t1950s_EAC_insig_std = nanstd(T_rate_EAC_insig)*multiplier;      
    % 1960s
    check = tt >= datenum(1960,01,01) & tt < datenum(1970,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t1960s_sig = nanmean(T_rate_sig)*multiplier;   
        trends(n_depth).t1960s_sig_std = nanstd(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t1960s_sig = NaN;   
        trends(n_depth).t1960s_sig_std = NaN;  
    end
    trends(n_depth).t1960s_insig = nanmean(T_rate_insig)*multiplier;    
    trends(n_depth).t1960s_insig_std = nanstd(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t1960s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
        trends(n_depth).t1960s_EAC_sig_std = nanstd(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t1960s_EAC_sig = NaN;  
        trends(n_depth).t1960s_EAC_sig_std = NaN;   
    end
    trends(n_depth).t1960s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;  
    trends(n_depth).t1960s_EAC_insig_std = nanstd(T_rate_EAC_insig)*multiplier;    
    % 1970s
    check = tt >= datenum(1970,01,01) & tt < datenum(1980,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t1970s_sig = nanmean(T_rate_sig)*multiplier; 
        trends(n_depth).t1970s_sig_std = nanstd(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t1970s_sig = NaN;   
        trends(n_depth).t1970s_sig_std = NaN;   
    end
    trends(n_depth).t1970s_insig = nanmean(T_rate_insig)*multiplier;    
    trends(n_depth).t1970s_insig_std = nanstd(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t1970s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;
        trends(n_depth).t1970s_EAC_sig_std = nanstd(T_rate_EAC_sig)*multiplier;
    else
        trends(n_depth).t1970s_EAC_sig = NaN;   
        trends(n_depth).t1970s_EAC_sig_std = NaN;   
    end
    trends(n_depth).t1970s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;    
    trends(n_depth).t1970s_EAC_insig_std = nanstd(T_rate_EAC_insig)*multiplier;    
    % 1980s
    check = tt >= datenum(1980,01,01) & tt < datenum(1990,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t1980s_sig = nanmean(T_rate_sig)*multiplier;   
        trends(n_depth).t1980s_sig_std = nanstd(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t1980s_sig = NaN;   
        trends(n_depth).t1980s_sig_std = NaN;   
    end
    trends(n_depth).t1980s_insig = nanmean(T_rate_insig)*multiplier; 
    trends(n_depth).t1980s_insig_std = nanstd(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t1980s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
        trends(n_depth).t1980s_EAC_sig_std = nanstd(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t1980s_EAC_sig = NaN;   
        trends(n_depth).t1980s_EAC_sig_std = NaN;   
    end
    trends(n_depth).t1980s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;      
    trends(n_depth).t1980s_EAC_insig_std = nanstd(T_rate_EAC_insig)*multiplier;    
    % 1990s
    check = tt >= datenum(1990,01,01) & tt < datenum(2000,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t1990s_sig = nanmean(T_rate_sig)*multiplier;   
        trends(n_depth).t1990s_sig_std = nanstd(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t1990s_sig = NaN;   
        trends(n_depth).t1990s_sig_std = NaN;   
    end
    trends(n_depth).t1990s_insig = nanmean(T_rate_insig)*multiplier;    
    trends(n_depth).t1990s_insig_std = nanstd(T_rate_insig)*multiplier;  
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t1990s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
        trends(n_depth).t1990s_EAC_sig_std = nanstd(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t1990s_EAC_sig = NaN;   
        trends(n_depth).t1990s_EAC_sig_std = NaN;   
    end
    trends(n_depth).t1990s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;   
    trends(n_depth).t1990s_EAC_insig_std = nanmean(T_rate_EAC_insig)*multiplier;  
    % 2000s
    check = tt >= datenum(2000,01,01) & tt < datenum(2010,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t2000s_sig = nanmean(T_rate_sig)*multiplier;   
        trends(n_depth).t2000s_sig_std = nanstd(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t2000s_sig = NaN;   
        trends(n_depth).t2000s_sig_std = NaN;   
    end
    trends(n_depth).t2000s_insig = nanmean(T_rate_insig)*multiplier;    
    trends(n_depth).t2000s_insig_std = nanstd(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t2000s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
        trends(n_depth).t2000s_EAC_sig_std = nanstd(T_rate_EAC_sig)*multiplier;  
    else
        trends(n_depth).t2000s_EAC_sig = NaN;   
        trends(n_depth).t2000s_EAC_sig_std = NaN;   
    end
    trends(n_depth).t2000s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;      
    trends(n_depth).t2000s_EAC_insig_std = nanstd(T_rate_EAC_insig)*multiplier;      
    % 2010s
    check = tt >= datenum(2010,01,01) & tt < datenum(2020,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t2010s_sig = nanmean(T_rate_sig)*multiplier; 
        trends(n_depth).t2010s_sig_std = nanstd(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t2010s_sig = NaN;   
        trends(n_depth).t2010s_sig_std = NaN;   
    end
    trends(n_depth).t2010s_insig = nanmean(T_rate_insig)*multiplier; 
    trends(n_depth).t2010s_insig_std = nanstd(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t2010s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
        trends(n_depth).t2010s_EAC_sig_std = nanstd(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t2010s_EAC_sig = NaN;   
        trends(n_depth).t2010s_EAC_sig_std = NaN;   
    end
    trends(n_depth).t2010s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;      
    trends(n_depth).t2010s_EAC_insig_std = nanstd(T_rate_EAC_insig)*multiplier;      
    % total change   
    trends(n_depth).total_insig = abs(tr(end)-tr(1));
    tr_sig = tr(tr_0_sig' == 1);
    if ~isempty(tr_sig)
        trends(n_depth).total_sig = abs(tr_sig(end)-tr_sig(1));
    else
        trends(n_depth).total_sig = NaN;
    end
    trends(n_depth).total_EAC_insig = abs(tr_EAC(end)-tr_EAC(1));
    tr_EAC_sig = tr_EAC(tr_EAC_0_sig' == 1);
    if ~isempty(tr_EAC_sig)
        trends(n_depth).total_EAC_sig = abs(tr_EAC_sig(end)-tr_EAC_sig(1));    
    else
        trends(n_depth).total_EAC_sig = NaN;
    end
    % significant trends from 1990s
    sig_trends = [trends(n_depth).t1990s_sig , trends(n_depth).t2000s_sig, trends(n_depth).t2010s_sig];
    sig_trends_EAC = [trends(n_depth).t1990s_EAC_sig , trends(n_depth).t2000s_EAC_sig, trends(n_depth).t2010s_EAC_sig];
    if sum(isnan(sig_trends)) == 0
        check = tt' > datenum(1990,01,01) & tr_0_sig == 1;
        tR = tr(check);
        trends(n_depth).total_sig_1990s = tR(end)-tR(1);
    else
        trends(n_depth).total_sig_1990s = NaN;
    end
    if sum(isnan(sig_trends_EAC)) == 0
        check = tt' > datenum(1990,01,01) & tr_EAC_0_sig == 1;
        tR = tr_EAC(check);
        trends(n_depth).total_EAC_sig_1990s = tR(end)-tR(1);
    else
        trends(n_depth).total_EAC_sig_1990s = NaN;
    end   
    % determine if > 10% of data is gap-filled
    trends(n_depth).data_coverage = 0;
    
    if 1 - sum(isfinite(data_struct.Tbin(n_depths,:))) / sum(isfinite(data_struct_server.Tbin(n_depths,:))) > 0.1
        trends(n_depth).data_coverage = 0.1;
    end    
    if 1 - sum(isfinite(data_struct.Tbin(n_depths,:))) / sum(isfinite(data_struct_server.Tbin(n_depths,:))) > 0.25
        trends(n_depth).data_coverage = 0.25;
    end
    if 1 - sum(isfinite(data_struct.Tbin(n_depths,:))) / sum(isfinite(data_struct_server.Tbin(n_depths,:))) > 0.5
        trends(n_depth).data_coverage = 0.50;
    end      
end
end