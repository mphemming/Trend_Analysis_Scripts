function [data, comb_file] = combine_files(site_str,data)

%% load in all available file splits

path = 'C:\Users\mphem\Documents\Work\UNSW\Trends\Data\Split_data_server_deseason\';
files = dir([path,site_str,'*trends*.mat']);

for n_file = 1:numel(files)
    fname = files(n_file).name;
    mat_files(n_file).data = load([path,fname]);
end

%% concatenate the files together

comb_file = mat_files(1).data;
% add conf limits for seasons
comb_file.EEMD_conf_std_limit = mat_files(2).data.EEMD_conf_std_limit;
comb_file.EEMD_conf_std_limit_EAC = mat_files(2).data.EEMD_conf_std_limit_EAC;
comb_file.EEMD_std_array = mat_files(2).data.EEMD_std_array;
comb_file.EEMD_std_array_EAC = mat_files(2).data.EEMD_std_array_EAC;
comb_file.EEMD_trend_sims = mat_files(2).data.EEMD_trend_sims;
comb_file.EEMD_trend_sims_EAC = mat_files(2).data.EEMD_trend_sims_EAC;
comb_file.EEMD_sims = mat_files(2).data.EEMD_sims;
comb_file.EEMD_conf_std_limit_Su = mat_files(2).data.EEMD_conf_std_limit_Su;
comb_file.EEMD_conf_std_limit_EAC_Su = mat_files(2).data.EEMD_conf_std_limit_EAC_Su;
comb_file.EEMD_conf_std_limit_Au = mat_files(3).data.EEMD_conf_std_limit_Au;
comb_file.EEMD_conf_std_limit_EAC_Au = mat_files(3).data.EEMD_conf_std_limit_EAC_Au;
comb_file.EEMD_conf_std_limit_Wi = mat_files(4).data.EEMD_conf_std_limit_Wi;
comb_file.EEMD_conf_std_limit_EAC_Wi = mat_files(4).data.EEMD_conf_std_limit_EAC_Wi;
comb_file.EEMD_conf_std_limit_Sp = mat_files(5).data.EEMD_conf_std_limit_Sp;
comb_file.EEMD_conf_std_limit_EAC_Sp = mat_files(5).data.EEMD_conf_std_limit_EAC_Sp;

%% Add / convert times

% EEMD time
if iscell(comb_file.EEMD_t)
    
    for nn = 1:numel(comb_file.EEMD_trend)
        a = comb_file.EEMD_t{nn};
        for t = 1:size(a,1)
            b = a(t,:);
            comb_file.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
        end
    end
    
    nt = size(data.t,2);
    for nn = 1:numel(comb_file.EEMD_trend)
        for t = 1:nt
            a = squeeze(vertcat(data.t(nn,t,:)))';
            data.t_conv(nn).t(t) = datenum(convertCharsToStrings(a));
        end
    end    
      
else
    
    for nn = 1:size(comb_file.EEMD_trend,1)
        a = squeeze(comb_file.EEMD_t(nn,:,:));
        for t = 1:size(a,1)
            b = a(t,:);
            comb_file.EEMD_t_conv(nn).t(t) = datenum(convertCharsToStrings(b));
        end
    end    
    
    nt = size(data.t,2);
    for nn = 1:size(data.t,1)
        for t = 1:nt
            a = squeeze(vertcat(data.t(nn,t,:)))';
            data.t_conv(nn).t(t) = datenum(convertCharsToStrings(a));
        end
    end
end


% get times for seasonal trends
[~,mn,dy,~,~,~] = datevec(comb_file.EEMD_t_conv(1).t);
% Summer
check_15_Su = mn == 12 | mn <= 2;
comb_file.EEMD_t_Su = data.t_conv(1).t(check_15_Su);
% Autumn
check_15_Au = mn > 2 & mn <= 5;
comb_file.EEMD_t_Au = data.t_conv(1).t(check_15_Au);
% Winter
check_15_Wi = mn > 5 & mn <= 8;
comb_file.EEMD_t_Wi = data.t_conv(1).t(check_15_Wi);
% Spring
check_15_Sp = mn > 8 & mn <= 11;
comb_file.EEMD_t_Sp = data.t_conv(1).t(check_15_Sp);

% % get time for std_array / conf interval
% [~,mn,dy,~,~,~] = datevec(data.t_conv(1).t))

end