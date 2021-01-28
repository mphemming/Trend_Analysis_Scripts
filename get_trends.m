function [trends] = get_trends(data_struct, trend_data_struct,trend_data_struct_server,multiplier)

n_depths = numel(trend_data_struct.EEMD_trend);


for n_depth = 1:n_depths
    
    sig = trend_data_struct_server.EEMD_conf_std_limit(1,:);
    tr = trend_data_struct.EEMD_trend{n_depth}; 
    tr_EAC = trend_data_struct.EEMD_trend_EAC{n_depth}; 
    tr_0 = tr-tr(1);    
    tr_0 = tr_0 / nanstd(trend_data_struct.EEMD_T{n_depth}); % normalise trend
    tr_EAC_0 = tr_EAC-tr_EAC(1);  
    tt = datenum(cell2mat(trend_data_struct.EEMD_t(n_depth)));
    % flag as NaN if insignificant trend value
    tr_0_sig = ones(size(tr_0));
    tr_EAC_0_sig = ones(size(tr_EAC_0));
    for n = 1:numel(tt)
        check = data_struct.t_conv(1).t == tt(n);
        if tr_0(n) < sig(check)
            tr_0_sig(n) = 0;
        end
        if tr_EAC_0(n) < sig(check)
            tr_EAC_0_sig(n) = 0;
        end        
    end
    % 1960s
    check = tt >= datenum(1960,01,01) & tt < datenum(1970,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t1960s_sig = nanmean(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t1960s_sig = NaN;   
    end
    trends(n_depth).t1960s_insig = nanmean(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t1960s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t1960s_EAC_sig = NaN;   
    end
    trends(n_depth).t1960s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;    
    % 1970s
    check = tt >= datenum(1970,01,01) & tt < datenum(1980,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t1970s_sig = nanmean(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t1970s_sig = NaN;   
    end
    trends(n_depth).t1970s_insig = nanmean(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t1970s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t1970s_EAC_sig = NaN;   
    end
    trends(n_depth).t1970s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;     
    % 1980s
    check = tt >= datenum(1980,01,01) & tt < datenum(1990,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t1980s_sig = nanmean(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t1980s_sig = NaN;   
    end
    trends(n_depth).t1980s_insig = nanmean(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t1980s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t1980s_EAC_sig = NaN;   
    end
    trends(n_depth).t1980s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;      
    % 1990s
    check = tt >= datenum(1990,01,01) & tt < datenum(2000,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t1990s_sig = nanmean(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t1990s_sig = NaN;   
    end
    trends(n_depth).t1990s_insig = nanmean(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t1990s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t1990s_EAC_sig = NaN;   
    end
    trends(n_depth).t1990s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;          
    % 2000s
    check = tt >= datenum(2000,01,01) & tt < datenum(2010,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t2000s_sig = nanmean(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t2000s_sig = NaN;   
    end
    trends(n_depth).t2000s_insig = nanmean(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t2000s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t2000s_EAC_sig = NaN;   
    end
    trends(n_depth).t2000s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;      
    % 2010s
    check = tt >= datenum(2010,01,01) & tt < datenum(2020,01,01);
    T_rate_sig = diff(tr(check & tr_0_sig' == 1));
    T_rate_insig = diff(tr(check));
    T_rate_EAC_sig = diff(tr_EAC(check & tr_EAC_0_sig' == 1));
    T_rate_EAC_insig = diff(tr_EAC(check));
    if sum(isfinite(T_rate_sig)) > 90
        trends(n_depth).t2010s_sig = nanmean(T_rate_sig)*multiplier;   
    else
        trends(n_depth).t2010s_sig = NaN;   
    end
    trends(n_depth).t2010s_insig = nanmean(T_rate_insig)*multiplier;    
    if sum(isfinite(T_rate_EAC_sig)) > 90
        trends(n_depth).t2010s_EAC_sig = nanmean(T_rate_EAC_sig)*multiplier;   
    else
        trends(n_depth).t2010s_EAC_sig = NaN;   
    end
    trends(n_depth).t2010s_EAC_insig = nanmean(T_rate_EAC_insig)*multiplier;      
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
end


end