function [Tchange, Tchange_EAC] = get_Tchange(data_struct, data_struct_server, trend_data_struct,trend_data_struct_server,start_year,end_year)

n_depths = numel(trend_data_struct.EEMD_trend);


for n_depth = 1:n_depths
    
    sig = trend_data_struct_server.EEMD_conf_std_limit(1,:);
    sig_EAC = trend_data_struct_server.EEMD_conf_std_limit_EAC(1,:);
    tt_sig = data_struct.t_conv(n_depth).t;
    tr = trend_data_struct.EEMD_trend{n_depth} / nanstd(trend_data_struct.EEMD_T{n_depth}); % normalise trend 
    tr_EAC = trend_data_struct.EEMD_trend_EAC{n_depth}  / nanstd(trend_data_struct.EEMD_T{n_depth}); % normalise trend
    tr_0 = tr-tr(1);  
    tr_EAC_0 = tr_EAC-tr_EAC(1);  
    tt = datenum(cell2mat(trend_data_struct.EEMD_t(n_depth)));
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
    % T change background
    check = tt' > datenum(start_year,01,01) & tr_0_sig == 1;
    tR = tr(check & isfinite(tr));
    if ~isempty(tR)
        Tchange(n_depth) = tR(end)-tR(1);
    else
        Tchange(n_depth) = NaN;
    end
    % T change EAC
    check = tt' > datenum(start_year,01,01) & tr_EAC_0_sig == 1;
    tR = tr_EAC(check & isfinite(tr_EAC));
    if ~isempty(tR)
        Tchange_EAC(n_depth) =  tR(end)-tR(1);
    else
        Tchange(n_depth) = NaN;
    end    
end
end