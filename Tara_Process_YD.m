% Variables to add manually

% change out dev file:
% acs043_0310 03/01/10-07/16/10
% acs007_0910 09/05/10-03/10/11
% acs091 03/11/11-08/12/11
% acs007_0611 08/14/11-10/26/11
% acs091_1028 11/24/11-12/22/11
% acs057 12/29/11-01/19/12

% if rejecting too many dissolved (green), raise <= value (begin at 0.001)
% plot to check this on is titled:
% checking to make sure we get good filtered measurements

% run unconstrained rows (all data)
% if need to remove some measurements
% split day into several intervals, run, then move completed files
% so they won't be overrun and repeat for rest of day's measurements, then
% put files together.

close all
clear all

yd = 196 % yearday
yd_str='196'; % yearday
year = 2010; % year
cruise='Mayotte-CT'; % cruise name
file_tsg = ['TSGdailyhex\DailyFull\TARATSG2010195.hex'];
mkdir('E:\Tara\Output',[num2str(year) '_' num2str(yd)]);
% Filter flush time (for 10 minute cycle, 6 gives 5 minutes of filtered data)
flushtime = 6 % dev file 007
% flushtime = 3
% flushtime = 6

% Convert yd,year to matlab datenum
dn = datenum(year, 0, yd);
datestr(dn)
[xx, month, day] = datevec(dn)

% Generate array of minutes for data binning
bin_dn = dn + datenum(0,0,0,0,0:(60*24-1),0)';

% Find files corresponding to this yearday
% **Making lots of assumptions about locations/naming of files, but should be ok.
% Note also that we look for both ACS .bin and .dat.

% files_flow = dir(['Flow/*' num2str(year) num2str(yd) '*.log']);
files_flow = dir(['Flow/*' num2str(year) yd_str '*.log']);

files_flow = {files_flow.name}';

files_acsd = dir(['ACS/*_' num2str(year) num2str(month, '%02u') num2str(day, '%02u') '*.dat']);
files_acsd = {files_acsd.name}';

files_acsb = dir(['ACS/*_' num2str(year) num2str(month, '%02u') num2str(day, '%02u') '*.bin']);
files_acsb = {files_acsb.name}';

files_acs = [files_acsb; files_acsd];

files_acs_isbin = logical([repmat(1,length(files_acsb),1); repmat(0,length(files_acsd),1)]);

% =========================================================================
% Flow file format
% 2009-09-28 13:02:00 UTC	-1	00000000	21.00	0.00	0.00	0.00	3.360	0.000	0.000	0.000
% 2009-09-28 13:02:01 UTC	-1	00000000	21.00	0.00	0.00	0.00	3.360	0.000	0.000	0.000
% 2009-09-28 13:02:02 UTC	-1	00000000	21.00	0.00	0.00	0.00	3.360	0.000	0.000	0.000
% 2009-09-28 13:02:03 UTC	-1	00000000	20.83	0.00	0.00	0.00	3.333	0.000	0.000	0.000
% 2009-09-28 13:02:04 UTC	-1	00000000	21.17	0.00	0.00	0.00	3.387	0.000	0.000	0.000
% =========================================================================

bin_flow = nan(length(bin_dn),2);
bin_flow_std = nan(length(bin_dn),2);
bin_valve = nan(length(bin_dn),1);

for kf = 1:length(files_flow)
    files_flow{kf}

    % Load in data file
    % ** Note use of HeaderLines to ignore first row of data in each file
    % ** Added because first date was showing up with three strange characters
    fid = fopen(['Flow/' files_flow{kf}]);
    dat = textscan(fid, ' %s %s UTC %f %s %f %f %f %f %f %f %f %f', 'CollectOutput',0, 'HeaderLines',1);
    %dat = textscan(fid, '%d-%d-%d %d:%d:%d UTC %f %s %f %f %f %f %f %f %f %f', 'CollectOutput',0);
    fclose(fid);
    numgoodrows= length(dat{12})
    % Need to override strcat whitespace behavior
    dat_dn = strcat(dat{1}(1:numgoodrows), 'x', dat{2}(1:numgoodrows));
    dat_dn = strrep(dat_dn, 'x', ' ');
    % Convert string date to numeric matlab datenum
    dat_dn = datenum(dat_dn);
    % Get valve state and flow
    dat_flow = [dat{5}(1:numgoodrows) dat{9}(1:numgoodrows)];
    dat_valve = dat{3}(1:numgoodrows);
    
    % Find first and last indices in to bin_dn where we'll place data from current file
    v = datevec(min(dat_dn));
    idx1 = find(bin_dn >= datenum([v(1:5) 0]), 1, 'first');
    v = datevec(max(dat_dn));
    idx2 = find(bin_dn <= datenum([v(1:5) 0]), 1, 'last');
    
    for kb = idx1:idx2
        idx_b = (dat_dn>=bin_dn(kb) & dat_dn<bin_dn(kb)+datenum(0,0,0,0,1,0));
        bin_flow(kb,:) = median(dat_flow(idx_b,:));
        bin_flow_std(kb,:) = std(dat_flow(idx_b,:));
        bin_valve(kb,1) = mean(dat_valve(idx_b));
    end
 end

% =========================================================================
% AC-S Data Files
% =========================================================================

% Read first file to get dimensions of cwl,awl
if files_acs_isbin(1)
    [xx,xx,xx,bin_rawcwl,bin_rawawl] = Read_ACSBin043_0310(['ACS/' files_acs{1}]);
else
    [xx,xx,xx,bin_rawcwl,bin_rawawl] = Read_ACSDat043_0310(['ACS/' files_acs{1}]);
end

bin_rawa = nan(length(bin_dn),length(bin_rawawl));
bin_rawc = nan(length(bin_dn),length(bin_rawcwl));
bin_rawa_std = nan(length(bin_dn),length(bin_rawawl));
bin_rawc_std = nan(length(bin_dn),length(bin_rawcwl));

for kf = 1:length(files_acs)
    files_acs{kf}

    % Load in data file
    if files_acs_isbin(kf)
        [dat_msec,dat_c,dat_a,xx,xx] = Read_ACSBin043_0310(['ACS/' files_acs{kf}]);
    else
        [dat_msec,dat_c,dat_a,xx,xx] = Read_ACSDat043_0310(['ACS/' files_acs{kf}]);
    end
    
    if isempty(dat_msec)|isempty(dat_c)|isempty(dat_a)
        continue
    end
    
    % Construct dn time vector using acs filename and msec timestamp in
    % file
    if files_acs_isbin(kf)
        idx = strfind(files_acs{kf}, '.bin');
    else
        idx = strfind(files_acs{kf}, '.dat');
    end
    dat_dn0 = datenum(str2num(files_acs{kf}((idx-14):(idx-11))), ...
        str2num(files_acs{kf}((idx-10):(idx-9))), str2num(files_acs{kf}((idx-8):(idx-7))), ...
        str2num(files_acs{kf}((idx-6):(idx-5))), str2num(files_acs{kf}((idx-4):(idx-3))), ...
        str2num(files_acs{kf}((idx-2):(idx-1))));
    dat_dn = dat_dn0 + datenum(0,0,0,0,0,(dat_msec-dat_msec(1))/1000);
    
    % Kludgie way to deal with the filter-wheel discontinuity--similar way I
    % dealt with it for Philex
    % ** a side discontinuity is between awl(46) and awl(47)
    % ** c side discontinuity is between cwl(41) and cwl(42)
    del_a = nan(size(dat_dn));
    del_c = nan(size(dat_dn));
    for k = 1:length(dat_dn)
        p = polyfit(bin_rawawl(37:41), dat_a(k,37:41), 1);
        del_a(k) = p(1)*bin_rawawl(42)+p(2) - dat_a(k,42);
        p = polyfit(bin_rawcwl(37:41), dat_c(k,37:41), 1);
        del_c(k) = p(1)*bin_rawcwl(42)+p(2) - dat_c(k,42);
    end
    del_a = del_a - nanmedian(del_a);
    del_c = del_c - nanmedian(del_c);
    dat_a(abs(del_a)>prctile(abs(del_a),95),:) = NaN;
    dat_c(abs(del_c)>prctile(abs(del_c),95),:) = NaN;
    
    % Find first and last indices in to bin_dn where we'll place data from current file 
    v = datevec(min(dat_dn));
    idx1 = find(bin_dn >= datenum([v(1:5) 0]), 1, 'first');
    v = datevec(max(dat_dn));
    idx2 = find(bin_dn <= datenum([v(1:5) 0]), 1, 'last');
    
    for kb = idx1:idx2
        idx_b = (dat_dn>=bin_dn(kb) & dat_dn<bin_dn(kb)+datenum(0,0,0,0,1,0));
        bin_rawa(kb,:) = nanmedian(dat_a(idx_b,:));
        bin_rawa_std(kb,:) = nanstd(dat_a(idx_b,:));
        bin_rawc(kb,:) =nanmedian(dat_c(idx_b,:));
        bin_rawc_std(kb,:) = nanstd(dat_c(idx_b,:));
    end
end

% acs data now loaded

% plots bin_rawc for assessing threshold
% figure, plot(bin_dn, bin_rawc(:,35), 'k.-'), 
% hold on, plot(bin_dn, bin_rawc_std(:,35), 'r.-'), 
% legend('raw c(536)','STD raw c(536)')
% datetick

% figure, plot(bin_dn, bin_rawa(:,12), 'b.-'), 
% hold on, plot(bin_dn, bin_rawa_std(:,12), 'r.-'), 
% legend('raw a(442)','STD raw a(442)')
% datetick

% save loaded data
% fname = ['Temp_Before_QC_' num2str(year) '_' num2str(yd) '.mat']
% save(fname, 'bin_*')

thresh_a = 0.015;
thresh_c = 0.03;

rejects = (bin_rawa_std(:,12) > thresh_a) | (bin_rawc_std(:,35) > thresh_c);

% shows rejected data (in squares)
% figure, hold on,
% plot(bin_dn, bin_rawa(:,12), 'b.-'), 
% plot(bin_dn, bin_rawc(:,35), 'k.-'), 
% plot(bin_dn(rejects), bin_rawa(rejects,12), 'rs'), 
% plot(bin_dn(rejects), bin_rawc(rejects,35), 'rs'), 
% legend('raw a(442)','raw c(536)', 'REJECT')
% datetick

bin_rawa(rejects,:) = nan;
bin_rawc(rejects,:) = nan;

bin_rawa(isnan(bin_valve),:) = nan;
bin_rawc(isnan(bin_valve),:) = nan;

% =========================================================================
% Truncate binned data to time range where we have ACS data
% =========================================================================

nodata = all(isnan(bin_rawa),2);
bin_dn(nodata,:) = [];
bin_flow(nodata,:) = [];
bin_flow_std(nodata,:) = [];
bin_valve(nodata,:) = [];
bin_rawa(nodata,:) = [];
bin_rawc(nodata,:) = [];
bin_rawa_std(nodata,:) = [];
bin_rawc_std(nodata,:) = [];

bin_part = true(size(bin_dn));
ixfiltstart = find(diff(bin_valve==-1)==1);
ixfiltend = find(diff(bin_valve==-1)==-1)+1;
 
 for ki = 1:(length(ixfiltstart)-45)

    % Don't consider first few minutes of filtered cycle (flushing)
     ixthisfiltstart = ixfiltstart(ki);
    
    % Find end of this filtered cycle, deal with possibility that the
    % data day ends during filtered cycle
     ixthisfiltend = ixfiltend(find(ixfiltend>(ixthisfiltstart+flushtime),1));
     if isempty(ixthisfiltend)
         ixthisfiltend = length(bin_dn);
     end
    
    if ixthisfiltend < 1
        ixthisfiltend = 1;
    end
    
    bin_filt((ixthisfiltstart+flushtime):ixthisfiltend) = 1;
    
    % Exclude filtered cycle and valve transition time from particle measurements
    ixexcludestart = ixthisfiltstart - 1;
    ixexcludeend = ixthisfiltend + 2;
    
    if ixexcludestart < 1
        ixexcludestart = 1;
    end
    
    if ixexcludeend > length(bin_dn)
        ixexcludeend = length(bin_dn);
    end
    
    bin_part(ixexcludestart:ixexcludeend) = 0;
    
 end
 
% =========================================================================
% BEGIN RERUNS HERE
% =========================================================================

% For ease of automated processing, also exclude particle measurements
% before first and after last filtered cycles
bin_part(1:ixfiltstart(1)) = 0;
bin_part(ixfiltend(end):end) = 0;

% =========================================================================
% Can we determine filtered measurements using std?
% =========================================================================

bin_rawag = nan(length(bin_dn),length(bin_rawawl));
bin_rawcg = nan(length(bin_dn),length(bin_rawcwl));

figure(1)
subplot(121),hold on
    plot(bin_dn, bin_rawa_std(:,20), 'b.-')
    plot(bin_dn(bin_valve==-1), bin_rawa_std(bin_valve==-1,20), 'ro')
    ylabel('a-side std')
title('Selecting stdev threshold for finding filtered measurements...')
subplot(122),hold on
    plot(bin_dn, bin_rawc_std(:,20), 'b.-')
    plot(bin_dn(bin_valve==-1), bin_rawc_std(bin_valve==-1,20), 'ro')
    ylabel('c-side std')


% ADJUST you can adjust <= value here if necessary

idx_ag = find(bin_valve==-1 & bin_rawa_std(:,20)<=0.004);
idx_cg = find(bin_valve==-1 & bin_rawc_std(:,20)<=0.004);

bin_rawag(idx_ag,:) = bin_rawa(idx_ag,:);
bin_rawcg(idx_cg,:) = bin_rawc(idx_cg,:);

figure(2) % good filtered data
subplot(121),hold on
    plot(bin_dn, bin_rawa(:,20), 'b.-')
    plot(bin_dn(idx_ag), bin_rawa(idx_ag,20), 'go')
    ylabel('a-side')
title('Checking to make sure we get good filtered measurements...')
subplot(122),hold on
    plot(bin_dn, bin_rawc(:,20), 'b.-')
    plot(bin_dn(idx_cg), bin_rawc(idx_cg,20), 'go')
    ylabel('c-side')

    print(figure(2),'-djpeg',['E:\Tara\Output\' num2str(year) '_' num2str(yd) '\filtered.jpeg'])
% if rejecting too many dissolved (green) in next plot, raise <= value
% and go back to ADJUST and rerun, otherwise stay with default of 0.004
% for both

save(['E:\Tara\Output\' num2str(year) '_' num2str(yd) '\' num2str(year) '_' num2str(yd) '.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUN TO HERE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IF data look good and you don't need to run in parts, 
% skip to RUN FIGURE 3

% PARTS ONLY -------------------------------------
% IF data look good but you need to run in parts, do this section,
% otherwise skip to FILTER ONLY:

% to find a specific point
[x,y]=ginput(1)
 
% to find all points below this value
a_below=find(bin_dn<x)
Time_division_point=a_below(end)
 
% to find all points above this value
a_after=find(bin_dn>x)
Time_division_point=a_after(1)
 
% to find points in between
a_after=find(bin_dn>x) %points above this
Time_division_point_end=a_after(1)
 
a_below=find(bin_dn<x) %points below this
Time_division_point_start=a_below(end)

% reject anything below this point 
rejects_ag = find(bin_dn<bin_dn(Time_division_point));
rejects_cg = find(bin_dn<bin_dn(Time_division_point));

% OR reject anything above this point
rejects_ag = find(bin_dn>bin_dn(Time_division_point));
rejects_cg = find(bin_dn>bin_dn(Time_division_point));

% OR reject anything below < and above > these points 
rejects_ag = find(bin_dn>bin_dn(Time_division_point_start) | bin_dn<bin_dn(Time_division_point_end)); 
rejects_cg = find(bin_dn>bin_dn(Time_division_point_start) | bin_dn<bin_dn(Time_division_point_end));

% GO TO RERUN FIGURE 2 -------------------------------------

% FILTER ONLY -------------------------------------
% OR if data need filtering (too many green in non-dissovled fraction), and
% you don't need to run in parts, run this section (adjusting for the
% > values), otherwise skip to FILTER & PARTS

rejects_ag = find(bin_rawag(:,20)>0.12); %edit the value based on filtered plot
rejects_cg = find(bin_rawcg(:,20)>0.04); %edit the value based on filtered plot

% GO TO RERUN FIGURE 2 -------------------------------------

% FILTER & PARTS -----------------------------------------
% OR if you need to filter the data (adjusting for the > values) AND run
% in parts, do this section then go to RERUN FIGURE 2

% to find a specific point
[x,y]=ginput(1)
 
% to find all points below this value
a_below=find(bin_dn<x)
Time_division_point=a_below(end)
 
% to find all points above this value
a_after=find(bin_dn>x)
Time_division_point=a_after(1)
 
% to find points in between
a_after=find(bin_dn>x) %points above this
Time_division_point_end=a_after(1)
 
a_below=find(bin_dn<x) %points below this
Time_division_point_start=a_below(end)

% reject anything below this point 
rejects_ag = find(bin_rawag(:,20)>0.12 | bin_dn<bin_dn(Time_division_point));
rejects_cg = find(bin_rawcg(:,20)>0.06 | bin_dn<bin_dn(Time_division_point));

% OR reject anything above this point
rejects_ag = find(bin_rawag(:,20)>0.12 | bin_dn>bin_dn(Time_division_point));
rejects_cg = find(bin_rawcg(:,20)>0.06 | bin_dn>bin_dn(Time_division_point));

% OR reject anything below < and above > these points 
rejects_ag = find(bin_rawag(:,20)>0.12 | bin_dn>bin_dn(Time_division_point_start) | bin_dn<bin_dn(Time_division_point_end)); 
rejects_cg = find(bin_rawcg(:,20)>0.06 | bin_dn>bin_dn(Time_division_point_start) | bin_dn<bin_dn(Time_division_point_end));

% RERUN FIGURE 2 -------------------------------------

close figure 2
figure(2) % good filtered data
subplot(121),hold on
    plot(bin_dn, bin_rawa(:,20), 'b.-')
    plot(bin_dn(idx_ag), bin_rawa(idx_ag,20), 'go')
    plot(bin_dn(rejects_ag), bin_rawag(rejects_ag,20), 'ro')
    ylabel('a-side')
    title('Checking to make sure we get good filtered measurements...')
subplot(122),hold on
    plot(bin_dn, bin_rawc(:,20), 'b.-')
    plot(bin_dn(idx_cg), bin_rawc(idx_cg,20), 'go')
    plot(bin_dn(rejects_cg), bin_rawcg(rejects_cg,20), 'ro')
    ylabel('c-side')

bin_rawag(rejects_ag,:) = nan;
bin_rawcg(rejects_cg,:) = nan;

% RUN FIGURE 3
filt_a = all(isfinite(bin_rawag),2);
filt_c = all(isfinite(bin_rawcg),2);
cgi = interp1(bin_dn(filt_c), bin_rawcg(filt_c,:), bin_dn, 'linear');
agi = interp1(bin_dn(filt_a), bin_rawag(filt_a,:), bin_dn, 'linear');

% NOTE: if rerunning this figure with different parameters, you must close
% it first
% close figure 3
figure(3)
subplot(121),hold on
    plot(bin_dn, bin_rawa(:,20), 'k.-')
    plot(bin_dn, agi(:,20), 'b.-')
    ylabel('a-side')
    datetick
    title('Make sure things are peachy...')
subplot(122),hold on
    plot(bin_dn, bin_rawc(:,20), 'k.-')
    plot(bin_dn, cgi(:,20), 'b.-')
    ylabel('c-side')
    datetick
    print(figure(3),'-djpeg',['E:\Tara\Output\' num2str(year) '_' num2str(yd) '\peachy.jpeg'])
% =========================================================================
% Load SBE 45 TSG data
% ** Note that filename YD seems to be one off from our yearday. Check this
% out with Herve Le Goff
% =========================================================================

dat = textread(file_tsg,'%s','delimiter','\n');
numheader = strmatch('*END*', dat); 

if isempty(numheader)
    numheader = 0;
end

fid = fopen(file_tsg);
dat = textscan(fid, 't1=%f, c1=%f, s=%f, t2=%f, lat=%u %f %c, lon=%u %f %c, hms=%6c, dmy=%6c', ...
    'HeaderLines',numheader);
fclose(fid);

tsg_dn = datenum(str2num(strcat('20', dat{12}(:,5:6))), str2num(dat{12}(:,3:4)), str2num(dat{12}(:,1:2)), ...
        str2num(dat{11}(:,1:2)), str2num(dat{11}(:,3:4)), str2num(dat{11}(:,5:6)));

tsg_tempc = [dat{1} dat{4}];
tsg_salin = dat{3};
tsg_lat = double(dat{5})+ dat{6}/60;
tsg_lon = double(dat{8})+ dat{9}/60;
% set sign appropriate for N/S and E/W
tsg_lat(dat{7}=='S') = -tsg_lat(dat{7}=='S');
tsg_lon(dat{10}=='W') = -tsg_lon(dat{10}=='W');

[b1, m1, n1] = unique(tsg_dn, 'first'); 
bin_tempc=interp1(tsg_dn(m1),tsg_tempc(m1),bin_dn,'linear');
bin_salin=interp1(tsg_dn(m1),tsg_salin(m1),bin_dn,'linear');
if(isempty(find(isfinite(bin_salin)))==1) 
    error('problem with TSG');
end
bin_lat=interp1(tsg_dn(m1),tsg_lat(m1),bin_dn,'linear');
bin_lon=interp1(tsg_dn(m1),tsg_lon(m1),bin_dn,'linear');

bin_cp = bin_rawc - cgi;
bin_ap_uncorr = bin_rawa - agi;

bin_cp(~bin_part,:) = NaN;
bin_ap_uncorr(~bin_part,:) = NaN;

% Deal with mismatch in spectral band positions between a and c measurements.
% Interpolate a onto c, limiting range to where the two overlap
bin_ap_uncorr = interp1(bin_rawawl, bin_ap_uncorr', bin_rawcwl, 'linear')';

bin_ap_std = interp1(bin_rawawl, bin_rawa_std', bin_rawcwl, 'linear')';
bin_cp_std = bin_rawc_std;

keeperwl = ~all(isnan(bin_ap_uncorr));
bin_wl = bin_rawcwl(keeperwl);
bin_cp = bin_cp(:,keeperwl);
bin_ap_uncorr = bin_ap_uncorr(:,keeperwl);

bin_cp_std=bin_cp_std(:,keeperwl);
bin_ap_std=bin_ap_std(:,keeperwl);

% Residual TS Correction (Slade et al 2009)
bin_ap = ACS_ResidTempScatCorr(bin_ap_uncorr,bin_cp, bin_wl);

% =========================================================================
% Save binned data to MAT file
% =========================================================================

% fname = ['Tara_ACS_' num2str(year) '_' num2str(yd) '.mat']
% save(fname, 'bin_*')

% =========================================================================
% Save subset to XLS and TXT files
% =========================================================================

sb_fname = ['Tara_ACS_apcp' num2str(year) '_' num2str(yd) '.xls'];
sb_fname_ascii = ['Tara_ACS_apcp' num2str(year) '_' num2str(yd)];

sb_hdr_xls = {'date','time','lat','lon','Wt','sal'};
sb_hdr_xls = [sb_hdr_xls strcat('ap', cellstr(num2str(bin_wl'))')];
sb_hdr_xls = [sb_hdr_xls strcat('cp', cellstr(num2str(bin_wl'))')];

% remove spaces in header fields
sb_hdr_xls=strrep(sb_hdr_xls, ' ', '');

%goodrows = logical(all(isfinite(bin_ap),2) .* isfinite(bin_salin) .* all(isfinite(bin_cp),2));
goodrows = logical(isfinite(bin_ap(:,20)) .* isfinite(bin_salin) .* isfinite(bin_cp(:,20)));
%goodrows = logical(isfinite(bin_ap(:,20)).* isfinite(bin_cp(:,20)));

% convert sb_datetime to yyyymmdd and HH:MM:SS for seabass submittal
sb_datetime = cellstr(datestr(bin_dn(goodrows)));

sb_datetime_date_xls=datestr(datenum(sb_datetime),'dd-mmm-yyyy');
sb_datetime_date=datestr(datenum(sb_datetime),'yyyymmdd');
sb_datetime_time=datestr(datenum(sb_datetime),'HH:MM:SS');

% =========================================================================
% export to excel
% =========================================================================

sb_dat_xls = [bin_lat bin_lon bin_tempc bin_salin bin_ap bin_cp];
sb_dat_xls = sb_dat_xls(goodrows,:);

[m,n]=size(sb_dat_xls);
tmp_xls = cell(m+1,n+2); %two more than n
tmp_xls(1,:) = sb_hdr_xls;

% add separated and reformatted date and time to tmp
tmp_xls(2:end,1) = cellstr(sb_datetime_date_xls);
tmp_xls(2:end,2) = cellstr(sb_datetime_time);
tmp_xls(2:end,3:end) = num2cell(sb_dat_xls);

% xlswrite(sb_fname,tmp_xls)

% =========================================================================
% generate seabass file
% =========================================================================

sb_hdr_ap = {'date','time','lat','lon','Wt','sal'};
sb_hdr_ap = [sb_hdr_ap strcat('ap', cellstr(num2str(bin_wl'))')];

sb_hdr_cp = {'date','time','lat','lon','Wt','sal'};
sb_hdr_cp = [sb_hdr_cp strcat('cp', cellstr(num2str(bin_wl'))')];

% remove spaces in header fields
sb_hdr_cp=strrep(sb_hdr_cp, ' ', '');
sb_hdr_ap=strrep(sb_hdr_ap, ' ', '');

sb_dat_ap = [bin_lat bin_lon bin_tempc bin_salin bin_ap];
sb_dat_ap = sb_dat_ap(goodrows,:);

sb_dat_cp = [bin_lat bin_lon bin_tempc bin_salin bin_cp];
sb_dat_cp = sb_dat_cp(goodrows,:);

[m,n]=size(sb_dat_ap);
tmp_ap = cell(m+1,n+2); % ap only
tmp_ap(1,:) = sb_hdr_ap;

[m,n]=size(sb_dat_cp);
tmp_cp = cell(m+1,n+2); % cp only
tmp_cp(1,:) = sb_hdr_cp;

tmp_date = cellstr(sb_datetime_date);
tmp_time = cellstr(sb_datetime_time);

% add separated and reformatted date and time to tmp
tmp_ap(2:end,1) = cellstr(sb_datetime_date);
tmp_ap(2:end,2) = cellstr(sb_datetime_time);
tmp_ap(2:end,3:end) = num2cell(sb_dat_ap);

tmp_cp(2:end,1) = cellstr(sb_datetime_date);
tmp_cp(2:end,2) = cellstr(sb_datetime_time);
tmp_cp(2:end,3:end) = num2cell(sb_dat_cp);

% get start and end dates and times
[x,y]=size(sb_datetime);
start_date=datestr(datenum(sb_datetime(1,1)),'yyyymmdd');
end_date=datestr(datenum(sb_datetime(x,1)),'yyyymmdd');
start_time=datestr(datenum(sb_datetime(1,1)),'HH:MM:SS');
end_time=datestr(datenum(sb_datetime(x,1)),'HH:MM:SS');

% =========================================================================
% make ap txt file
% =========================================================================

% name file
% fpath='d:\misclab\tara\Tara Processing\working\';
fpath=['E:\Tara\Output\' num2str(year) '_' num2str(yd) '\'];
fname=strcat(fpath,sb_fname_ascii);
extension = 'ap.txt';
fid = fopen(strcat(fname,extension), 'wt');

fprintf(fid,'/begin_header\n');
fprintf(fid,'/investigators=Emmanuel_Boss,Wayne_Slade,Lisa_Taylor\n');
fprintf(fid,'/affiliations=UMaine-MISC_Lab,UMaine-MISC_Lab,UMaine-MISC_Lab\n');
fprintf(fid,'/contact=emmanuel.boss@maine.edu\n');
fprintf(fid,'/experiment=TARA_expedition\n');
fprintf(fid,'/cruise=');
fprintf(fid,cruise);
fprintf(fid,'\n');
fprintf(fid,'/station=NA\n');
fprintf(fid,'/data_file_name=');
fprintf(fid,strcat(sb_fname_ascii,extension));
fprintf(fid,'\n');
fprintf(fid,'/documents=TARA_doc.txt\n');
fprintf(fid,'/calibration_files=TARA_doc.txt\n');
fprintf(fid,'/data_type=flow_thru\n');
fprintf(fid,'/data_status=final\n');
fprintf(fid, '/start_date=%s\n',start_date);
fprintf(fid, '/end_date=%s\n',end_date);
fprintf(fid,'/start_time=%s[GMT]\n',start_time);
fprintf(fid,'/end_time=%s[GMT]\n',end_time);
fprintf(fid,'/north_latitude=%5.3f[DEG]\n',max(sb_dat_ap(:,1)));
fprintf(fid,'/south_latitude=%5.3f[DEG]\n',min(sb_dat_ap(:,1)));
fprintf(fid,'/east_longitude=%5.3f[DEG]\n',max(sb_dat_ap(:,2)));
fprintf(fid,'/west_longitude=%5.3f[DEG]\n',min(sb_dat_ap(:,2)));
fprintf(fid,'/water_depth=NA\n');
fprintf(fid,'/measurement_depth=1.5\n'); % not allowed if depth is in the data
fprintf(fid,'/secchi_depth=NA\n');
fprintf(fid,'/cloud_percent=NA\n');
fprintf(fid,'/wind_speed=NA\n');
fprintf(fid,'/wave_height=NA\n');
fprintf(fid,'/missing=-9999\n');
fprintf(fid,'/delimiter=space\n');

fprintf(fid,'/fields=');
  for j=1:length(sb_hdr_ap)
      fprintf(fid, '%s,',cell2mat(sb_hdr_ap(j)));
  end
fprintf(fid,'\n');
fprintf(fid,'/units=yyyymmdd,hh:mm:ss,degrees,degrees,degreesC,PSU,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m\n');
fprintf(fid,'/end_header\n');    

[m, n]=size(sb_dat_ap);
for i=1:m
     fprintf(fid, '%8s ',cell2mat(tmp_date(i,1)));
     fprintf(fid, '%8s ',cell2mat(tmp_time(i,1)));
     % fprintf(fid, '%20s ',cell2mat(sb_datetime(i)));
     fprintf(fid,'%6.4f ',sb_dat_ap(i,:));
     fprintf(fid,'\n');
end     
fclose(fid);

% =========================================================================
% make cp txt file
% =========================================================================

% name file
fpath=['E:\Tara\Output\' num2str(year) '_' num2str(yd) '\'];
fname=strcat(fpath,sb_fname_ascii);
extension = 'cp.txt';
fid = fopen(strcat(fname,extension), 'wt');

fprintf(fid,'/begin_header\n');
fprintf(fid,'/investigators=Emmanuel_Boss,Wayne_Slade,Lisa_Taylor\n');
fprintf(fid,'/affiliations=UMaine-MISC_Lab,UMaine-MISC_Lab,UMaine-MISC_Lab\n');
fprintf(fid,'/contact=emmanuel.boss@maine.edu\n');
fprintf(fid,'/experiment=TARA_expedition\n');
fprintf(fid,'/cruise=');
fprintf(fid,cruise);
fprintf(fid,'\n');
fprintf(fid,'/station=NA\n');
fprintf(fid,'/data_file_name=');
fprintf(fid,strcat(sb_fname_ascii,extension));
fprintf(fid,'\n');
fprintf(fid,'/documents=TARA_doc.txt\n');
fprintf(fid,'/calibration_files=TARA_doc.txt\n');
fprintf(fid,'/data_type=flow_thru\n');
fprintf(fid,'/data_status=final\n');
fprintf(fid, '/start_date=%s\n',start_date);
fprintf(fid, '/end_date=%s\n',end_date);
fprintf(fid,'/start_time=%s[GMT]\n',start_time);
fprintf(fid,'/end_time=%s[GMT]\n',end_time);
fprintf(fid,'/north_latitude=%5.3f[DEG]\n',max(sb_dat_cp(:,1)));
fprintf(fid,'/south_latitude=%5.3f[DEG]\n',min(sb_dat_cp(:,1)));
fprintf(fid,'/east_longitude=%5.3f[DEG]\n',max(sb_dat_cp(:,2)));
fprintf(fid,'/west_longitude=%5.3f[DEG]\n',min(sb_dat_cp(:,2)));
fprintf(fid,'/water_depth=NA\n');
fprintf(fid,'/measurement_depth=1.5\n'); % not allowed if depth is in the data
fprintf(fid,'/secchi_depth=NA\n');
fprintf(fid,'/cloud_percent=NA\n');
fprintf(fid,'/wind_speed=NA\n');
fprintf(fid,'/wave_height=NA\n');
fprintf(fid,'/missing=-9999\n');
fprintf(fid,'/delimiter=space\n');

fprintf(fid,'/fields=');
  for j=1:length(sb_hdr_cp)
      fprintf(fid, '%s,',cell2mat(sb_hdr_cp(j)));
  end
fprintf(fid,'\n');
fprintf(fid,'/units=yyyymmdd,hh:mm:ss,degrees,degrees,degreesC,PSU,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m\n');
fprintf(fid,'/end_header\n');    

[m, n]=size(sb_dat_cp);
for i=1:m
     fprintf(fid, '%8s ',cell2mat(tmp_date(i,1)));
     fprintf(fid, '%8s ',cell2mat(tmp_time(i,1)));
     fprintf(fid,'%6.4f ',sb_dat_cp(i,:));
     fprintf(fid,'\n');
end     
fclose(fid);

% =========================================================================
% generate seabass error files
% =========================================================================

sb_hdr_ap = {'date','time','lat','lon','Wt','sal'};
sb_hdr_ap = [sb_hdr_ap strcat('ap', cellstr(num2str(bin_wl'))'), 'ap730'];

sb_hdr_cp = {'date','time','lat','lon','Wt','sal'};
sb_hdr_cp = [sb_hdr_cp strcat('cp', cellstr(num2str(bin_wl'))')];

% remove spaces in header fields
sb_hdr_ap=strrep(sb_hdr_ap, ' ', '');
sb_hdr_cp=strrep(sb_hdr_cp, ' ', '');

bin_ap_730=interp1(bin_wl,bin_ap',730,'linear');

sb_dat_ap_std = [bin_lat bin_lon bin_tempc bin_salin bin_ap_std bin_ap_730'];
sb_dat_ap_std = sb_dat_ap_std(goodrows,:);

sb_dat_cp_std = [bin_lat bin_lon bin_tempc bin_salin bin_cp_std];
sb_dat_cp_std = sb_dat_cp_std(goodrows,:);

[m,n]=size(sb_dat_ap_std);
tmp_ap = cell(m+1,n+2); % ap only
tmp_ap(1,:) = sb_hdr_ap;

[m,n]=size(sb_dat_cp_std);
tmp_cp = cell(m+1,n+2); % cp only
tmp_cp(1,:) = sb_hdr_cp;

tmp_date = cellstr(sb_datetime_date);
tmp_time = cellstr(sb_datetime_time);

% add separated and reformatted date and time to tmp
tmp_ap(2:end,1) = cellstr(sb_datetime_date);
tmp_ap(2:end,2) = cellstr(sb_datetime_time);
tmp_ap(2:end,3:end) = num2cell(sb_dat_ap_std);

tmp_cp(2:end,1) = cellstr(sb_datetime_date);
tmp_cp(2:end,2) = cellstr(sb_datetime_time);
tmp_cp(2:end,3:end) = num2cell(sb_dat_cp_std);

% get start and end dates and times
[x,y]=size(sb_datetime);
start_date=datestr(datenum(sb_datetime(1,1)),'yyyymmdd');
end_date=datestr(datenum(sb_datetime(x,1)),'yyyymmdd');
start_time=datestr(datenum(sb_datetime(1,1)),'HH:MM:SS');
end_time=datestr(datenum(sb_datetime(x,1)),'HH:MM:SS');

% =========================================================================
% make ap error txt file
% =========================================================================

% name file
fpath=['E:\Tara\Output\' num2str(year) '_' num2str(yd) '\'];
fname=strcat(fpath,sb_fname_ascii);
extension = 'ap_uncertainty.txt';
fid = fopen(strcat(fname,extension), 'wt');

fprintf(fid,'/begin_header\n');
fprintf(fid,'/investigators=Emmanuel_Boss,Wayne_Slade,Lisa_Taylor\n');
fprintf(fid,'/affiliations=UMaine-MISC_Lab,UMaine-MISC_Lab,UMaine-MISC_Lab\n');
fprintf(fid,'/contact=emmanuel.boss@maine.edu\n');
fprintf(fid,'/experiment=TARA_expedition\n');
fprintf(fid,'/cruise=');
fprintf(fid,cruise);
fprintf(fid,'\n');
fprintf(fid,'/station=NA\n');
fprintf(fid,'/data_file_name=');
fprintf(fid,strcat(sb_fname_ascii,extension));
fprintf(fid,'\n');
fprintf(fid,'/documents=TARA_doc.txt\n');
fprintf(fid,'/calibration_files=TARA_doc.txt\n');
fprintf(fid,'/data_type=flow_thru\n');
fprintf(fid,'/data_status=final\n');
fprintf(fid, '/start_date=%s\n',start_date);
fprintf(fid, '/end_date=%s\n',end_date);
fprintf(fid,'/start_time=%s[GMT]\n',start_time);
fprintf(fid,'/end_time=%s[GMT]\n',end_time);
fprintf(fid,'/north_latitude=%5.3f[DEG]\n',max(sb_dat_ap(:,1)));
fprintf(fid,'/south_latitude=%5.3f[DEG]\n',min(sb_dat_ap(:,1)));
fprintf(fid,'/east_longitude=%5.3f[DEG]\n',max(sb_dat_ap(:,2)));
fprintf(fid,'/west_longitude=%5.3f[DEG]\n',min(sb_dat_ap(:,2)));
fprintf(fid,'/water_depth=NA\n');
fprintf(fid,'/measurement_depth=1.5\n'); % not allowed if depth is in the data
fprintf(fid,'/secchi_depth=NA\n');
fprintf(fid,'/cloud_percent=NA\n');
fprintf(fid,'/wind_speed=NA\n');
fprintf(fid,'/wave_height=NA\n');
fprintf(fid,'/missing=-9999\n');
fprintf(fid,'/delimiter=space\n');

fprintf(fid,'/fields=');
  for j=1:length(sb_hdr_ap)
      fprintf(fid, '%s,',cell2mat(sb_hdr_ap(j)));
  end
fprintf(fid,'\n');
fprintf(fid,'/units=yyyymmdd,hh:mm:ss,degrees,degrees,degreesC,PSU,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m\n');
fprintf(fid,'/end_header\n');    

[m, n]=size(sb_dat_ap);
for i=1:m
     fprintf(fid, '%8s ',cell2mat(tmp_date(i,1)));
     fprintf(fid, '%8s ',cell2mat(tmp_time(i,1)));
     % fprintf(fid, '%20s ',cell2mat(sb_datetime(i)));
     fprintf(fid,'%6.4f ',sb_dat_ap_std(i,:));
     fprintf(fid,'\n');
end     
fclose(fid);

% =========================================================================
% make cp error txt file
% =========================================================================

% name file
fpath=['E:\Tara\Output\' num2str(year) '_' num2str(yd) '\'];
fname=strcat(fpath,sb_fname_ascii);
extension = 'cp_uncertainty.txt';
fid = fopen(strcat(fname,extension), 'wt');

fprintf(fid,'/begin_header\n');
fprintf(fid,'/investigators=Emmanuel_Boss,Wayne_Slade,Lisa_Taylor\n');
fprintf(fid,'/affiliations=UMaine-MISC_Lab,UMaine-MISC_Lab,UMaine-MISC_Lab\n');
fprintf(fid,'/contact=emmanuel.boss@maine.edu\n');
fprintf(fid,'/experiment=TARA_expedition\n');
fprintf(fid,'/cruise=');
fprintf(fid,cruise);
fprintf(fid,'\n');
fprintf(fid,'/station=NA\n');
fprintf(fid,'/data_file_name=');
fprintf(fid,strcat(sb_fname_ascii,extension));
fprintf(fid,'\n');
fprintf(fid,'/documents=TARA_doc.txt\n');
fprintf(fid,'/calibration_files=TARA_doc.txt\n');
fprintf(fid,'/data_type=flow_thru\n');
fprintf(fid,'/data_status=final\n');
fprintf(fid, '/start_date=%s\n',start_date);
fprintf(fid, '/end_date=%s\n',end_date);
fprintf(fid,'/start_time=%s[GMT]\n',start_time);
fprintf(fid,'/end_time=%s[GMT]\n',end_time);
fprintf(fid,'/north_latitude=%5.3f[DEG]\n',max(sb_dat_cp(:,1)));
fprintf(fid,'/south_latitude=%5.3f[DEG]\n',min(sb_dat_cp(:,1)));
fprintf(fid,'/east_longitude=%5.3f[DEG]\n',max(sb_dat_cp(:,2)));
fprintf(fid,'/west_longitude=%5.3f[DEG]\n',min(sb_dat_cp(:,2)));
fprintf(fid,'/water_depth=NA\n');
fprintf(fid,'/measurement_depth=1.5\n'); % not allowed if depth is in the data
fprintf(fid,'/secchi_depth=NA\n');
fprintf(fid,'/cloud_percent=NA\n');
fprintf(fid,'/wind_speed=NA\n');
fprintf(fid,'/wave_height=NA\n');
fprintf(fid,'/missing=-9999\n');
fprintf(fid,'/delimiter=space\n');

fprintf(fid,'/fields=');
  for j=1:length(sb_hdr_cp)
      fprintf(fid, '%s,',cell2mat(sb_hdr_cp(j)));
  end
fprintf(fid,'\n');
fprintf(fid,'/units=yyyymmdd,hh:mm:ss,degrees,degrees,degreesC,PSU,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m,1/m\n');
fprintf(fid,'/end_header\n');    

[m, n]=size(sb_dat_cp);
for i=1:m
     fprintf(fid, '%8s ',cell2mat(tmp_date(i,1)));
     fprintf(fid, '%8s ',cell2mat(tmp_time(i,1)));
     fprintf(fid,'%6.4f ',sb_dat_cp_std(i,:));
     fprintf(fid,'\n');
end     
fclose(fid);

figure(4)
plot(bin_wl, bin_ap)
title('wl & ap');

print(figure(4),'-djpeg',['E:\Tara\Output\' num2str(year) '_' num2str(yd) '\wlap.jpeg'])
 
figure(5)
plot(bin_wl, bin_cp)
title('wl & cp');

print(figure(5),'-djpeg',['E:\Tara\Output\' num2str(year) '_' num2str(yd) '\wlcp.jpeg'])

figure(6)
plot(bin_dn, bin_flow)
title('dn & flow');

