%% =============================================
%  HYCOM vs ARGO VALIDATION SCRIPT
%
%  Validates the HYCOM GLBy0.08 upper-ocean state against
%  quality-controlled Argo delayed-mode profiles for the four
%  North Indian Ocean tropical cyclones analysed in
%  Singh & Behera (2026).
%
%  Set STORM_ID below to one of: KYARR | AMPHAN | FANI | TAUKTAE
%
%  Method : Spatiotemporal nearest-neighbour colocation (50 km / 3 h)
%  Metrics: Bias, RMSE, MAE, Correlation, Murphy Skill Score
%  Fields : Temperature, Salinity, Mixed-layer depth (MLD)
%
%  MLD is computed identically for Argo and HYCOM from each
%  temperature profile using a 0.2 degC threshold referenced to
%  5 m (de Boyer Montegut et al., 2004), so the comparison is
%  strictly like-for-like.
%
%  Output (per storm): <out_dir>/validation_full_report.txt
%                      <out_dir>/stats_temperature.csv
%                      <out_dir>/stats_salinity.csv
% ============================================================
clear; clc; close all;

%% ── USER SETTINGS ───────────
STORM_ID = 'TAUKTAE';      % KYARR | AMPHAN | FANI | TAUKTAE

% Root folder holding the input data, organised per the README
% "Expected Data Layout":
%   <DATA_ROOT>/Argo/<Storm>/   daily   YYYYMMDD_prof.nc
%   <DATA_ROOT>/HYCOM/<Storm>/  hourly  <Storm>_DATA_YYYYMMDD_HH.nc
% Use the relative 'data' folder for the archived layout, or set an
% absolute path here for a local run.
DATA_ROOT = fullfile('data');

max_dist_km  = 50;     % spatial colocation radius (km)
max_dt_hours = 3;      % time colocation window (hours)
stat_depths  = [5 10 20 30 50 75 100 150 200 300 500]; % standard depths (m)
argo_ref     = datenum(1950,1,1);   % Argo JULD reference epoch
FILL_VAL     = 99999;               % Argo standard fill value

% Per-storm metadata: HYCOM filename prefix + active-period label
switch upper(STORM_ID)
    case 'KYARR',   cfg.prefix = 'Kyarr';   cfg.label = 'October-November 2019';
    case 'AMPHAN',  cfg.prefix = 'Amphan';  cfg.label = 'May 2020';
    case 'FANI',    cfg.prefix = 'Fani';    cfg.label = 'April-May 2019';
    case 'TAUKTAE', cfg.prefix = 'Tauktae'; cfg.label = 'May 2021';
    otherwise, error('Unknown STORM_ID: %s', STORM_ID);
end

argo_dir  = fullfile(DATA_ROOT, 'Argo',  cfg.prefix);
hycom_dir = fullfile(DATA_ROOT, 'HYCOM', cfg.prefix);
out_dir   = fullfile('Validation', lower(cfg.prefix));
%% ─────────────────────────────────────────────────────────

if ~exist(out_dir,'dir'), mkdir(out_dir); end

fprintf('\n================  %s  (%s)  ================\n\n', upper(STORM_ID), cfg.label);

%% ══════════════════════════════
%  STEP 1 – SCAN HYCOM FILES & READ STATIC GRID
%  (must come first so the grid is available during Argo QC)
%% ═══════════════════════════════════════════════
fprintf('[1/5] Scanning HYCOM files and reading grid...\n');

hpattern = [cfg.prefix '_DATA_????????_??.nc'];
flist_h = dir(fullfile(hycom_dir, hpattern));
if isempty(flist_h)
    error('No HYCOM files found in: %s\nExpected pattern: %s', hycom_dir, hpattern);
end

% Parse timestamps from filenames
nf = numel(flist_h);
hycom_dnum = NaN(nf, 1);
for k = 1:nf
    tok = regexp(flist_h(k).name, '(\d{8})_(\d{2})', 'tokens', 'once');
    if isempty(tok)
        warning('Cannot parse timestamp from: %s', flist_h(k).name);
        continue;
    end
    yr = str2double(tok{1}(1:4));
    mo = str2double(tok{1}(5:6));
    dy = str2double(tok{1}(7:8));
    hr = str2double(tok{2});
    hycom_dnum(k) = datenum(yr, mo, dy, hr, 0, 0);
end
[hycom_dnum, sidx] = sort(hycom_dnum);
flist_h = flist_h(sidx);

fprintf('   HYCOM files  : %d\n', nf);
fprintf('   HYCOM period : %s  ->  %s\n\n', ...
    datestr(hycom_dnum(1),   'yyyy-mmm-dd HH:MM'), ...
    datestr(hycom_dnum(end), 'yyyy-mmm-dd HH:MM'));

% Read static HYCOM grid from first file
f1    = fullfile(hycom_dir, flist_h(1).name);
lon_h = double(ncread(f1, 'lon'));
lat_h = double(ncread(f1, 'lat'));
dep_h = double(ncread(f1, 'depth'));
[LON_H, LAT_H] = meshgrid(lon_h, lat_h);   % [N_lat x N_lon]

lon_h_min = min(lon_h);  lon_h_max = max(lon_h);
lat_h_min = min(lat_h);  lat_h_max = max(lat_h);

fprintf('   HYCOM domain : Lon %.2f-%.2f E  |  Lat %.2f-%.2f N\n\n', ...
    lon_h_min, lon_h_max, lat_h_min, lat_h_max);

%% ══════════════════════════════════════════════════
%  STEP 2 – READ & MERGE DAILY ARGO FILES WITH FULL QC
%% ═══════════════════════════════════════
fprintf('[2/5] Reading daily Argo files...\n');

flist_a = dir(fullfile(argo_dir, '????????_prof.nc'));
if isempty(flist_a)
    error('No daily Argo files found in: %s\nExpected names like 20210510_prof.nc', argo_dir);
end
fprintf('   Found %d daily file(s)\n', numel(flist_a));

all_juld = [];
all_lat  = [];
all_lon  = [];
all_pres = [];
all_temp = [];
all_psal = [];

for kf = 1:numel(flist_a)
    fn = fullfile(argo_dir, flist_a(kf).name);
    fprintf('   Reading %s\n', flist_a(kf).name);

    % ── Position / time ──
    juld_k  = double(ncread(fn, 'JULD'));
    lat_k   = double(ncread(fn, 'LATITUDE'));
    lon_k   = double(ncread(fn, 'LONGITUDE'));
    pos_qc  = ncread(fn, 'POSITION_QC');   % char [N_PROF]
    juld_qc = ncread(fn, 'JULD_QC');       % char [N_PROF]

    % ── Prefer ADJUSTED fields; fall back to raw ──
    try
        pres_k  = double(ncread(fn, 'PRES_ADJUSTED'));
        temp_k  = double(ncread(fn, 'TEMP_ADJUSTED'));
        psal_k  = double(ncread(fn, 'PSAL_ADJUSTED'));
        pres_qc = ncread(fn, 'PRES_ADJUSTED_QC');
        temp_qc = ncread(fn, 'TEMP_ADJUSTED_QC');
        psal_qc = ncread(fn, 'PSAL_ADJUSTED_QC');
        fprintf('     -> using ADJUSTED fields\n');
    catch
        pres_k  = double(ncread(fn, 'PRES'));
        temp_k  = double(ncread(fn, 'TEMP'));
        psal_k  = double(ncread(fn, 'PSAL'));
        pres_qc = ncread(fn, 'PRES_QC');
        temp_qc = ncread(fn, 'TEMP_QC');
        psal_qc = ncread(fn, 'PSAL_QC');
        fprintf('     -> using RAW fields (adjusted unavailable)\n');
    end

    % ── Replace fill values with NaN ──
    juld_k(juld_k >= FILL_VAL) = NaN;
    lat_k( lat_k  >= FILL_VAL) = NaN;
    lon_k( lon_k  >= FILL_VAL) = NaN;
    pres_k(pres_k >= FILL_VAL) = NaN;
    temp_k(temp_k >= FILL_VAL) = NaN;
    psal_k(psal_k >= FILL_VAL) = NaN;

    % ── Level QC: keep only flags '1' (good) and '2' (probably good) ──
    pres_k(~(pres_qc == '1' | pres_qc == '2')) = NaN;
    temp_k(~(temp_qc == '1' | temp_qc == '2')) = NaN;
    psal_k(~(psal_qc == '1' | psal_qc == '2')) = NaN;

    % ── Profile-level QC ──
    good_pos  = (pos_qc  == '1' | pos_qc  == '2') & ~isnan(lat_k) & ~isnan(lon_k);
    good_time = (juld_qc == '1' | juld_qc == '2') & ~isnan(juld_k);
    n_valid_T = sum(~isnan(temp_k), 1)' >= 5;   % at least 5 valid T levels
    good_prof = good_pos & good_time & n_valid_T;

    % ── Keep only profiles inside HYCOM domain ──
    in_domain = lat_k >= lat_h_min & lat_k <= lat_h_max & ...
                lon_k >= lon_h_min & lon_k <= lon_h_max;
    good_prof = good_prof & in_domain;

    if sum(good_prof) == 0
        fprintf('     -> 0 valid profiles in HYCOM domain; skipping.\n');
        continue;
    end

    new_pres = pres_k(:, good_prof);
    new_temp = temp_k(:, good_prof);
    new_psal = psal_k(:, good_prof);

    % ── Concatenate — pad rows if level count differs ──
    n_cur = size(all_pres, 1);
    n_new = size(new_pres, 1);
    if n_cur == 0
        all_pres = new_pres;
        all_temp = new_temp;
        all_psal = new_psal;
    elseif n_new > n_cur
        pad = NaN(n_new - n_cur, size(all_pres, 2));
        all_pres = [[all_pres; pad], new_pres];
        all_temp = [[all_temp; pad], new_temp];
        all_psal = [[all_psal; pad], new_psal];
    elseif n_new < n_cur
        pad = NaN(n_cur - n_new, size(new_pres, 2));
        all_pres = [all_pres, [new_pres; pad]];
        all_temp = [all_temp, [new_temp; pad]];
        all_psal = [all_psal, [new_psal; pad]];
    else
        all_pres = [all_pres, new_pres];
        all_temp = [all_temp, new_temp];
        all_psal = [all_psal, new_psal];
    end

    all_juld = [all_juld; juld_k(good_prof)];   %#ok<*AGROW>
    all_lat  = [all_lat;  lat_k(good_prof)];
    all_lon  = [all_lon;  lon_k(good_prof)];
end

if isempty(all_juld)
    error('No valid Argo profiles found inside the HYCOM domain. Check file paths and domain bounds.');
end

% Working variables
lat_a     = all_lat;
lon_a     = all_lon;
pres      = all_pres;   % [N_LEVELS x N_PROF]
temp_a    = all_temp;
psal_a    = all_psal;
argo_dnum = argo_ref + all_juld;
nprof     = numel(all_juld);

fprintf('\n   Total profiles (QC + in domain) : %d\n', nprof);
fprintf('   Date range : %s  ->  %s\n', ...
    datestr(min(argo_dnum), 'yyyy-mmm-dd'), datestr(max(argo_dnum), 'yyyy-mmm-dd'));
fprintf('   Lat : %.3f  ->  %.3f N\n', min(lat_a), max(lat_a));
fprintf('   Lon : %.3f  ->  %.3f E\n\n', min(lon_a), max(lon_a));

%% ═════════════════════════════════════════
%  STEP 3 – COLOCATION
%% ═════════════════════════════════════════
fprintf('[3/5] Colocating Argo profiles with HYCOM...\n');

nz = length(stat_depths);
co_temp_argo  = NaN(nz, nprof);
co_temp_hycom = NaN(nz, nprof);
co_psal_argo  = NaN(nz, nprof);
co_psal_hycom = NaN(nz, nprof);
co_mld_argo   = NaN(nprof, 1);
co_mld_hycom  = NaN(nprof, 1);
co_dist_km    = NaN(nprof, 1);
co_dt_hours   = NaN(nprof, 1);
co_matched    = false(nprof, 1);

% Haversine distance function (returns km)
hav = @(la1, lo1, la2, lo2) 2 * 6371 * asin(sqrt( ...
    sind((la2 - la1) / 2).^2 + ...
    cosd(la1) .* cosd(la2) .* sind((lo2 - lo1) / 2).^2));

for ip = 1:nprof
    if mod(ip, 50) == 0 || ip == nprof
        fprintf('   Processing profile %d / %d\r', ip, nprof);
    end

    % ── Find nearest HYCOM time step ──
    [dt_days, fi] = min(abs(hycom_dnum - argo_dnum(ip)));
    dt_hrs = dt_days * 24;
    if dt_hrs > max_dt_hours, continue; end

    % ── Find nearest HYCOM grid point ──
    dist = hav(lat_a(ip), lon_a(ip), LAT_H, LON_H);   % [N_lat x N_lon]
    [dmin, idx] = min(dist(:));
    if dmin > max_dist_km, continue; end
    [irow, icol] = ind2sub(size(dist), idx);   % irow = lat index, icol = lon index

    % ── Read HYCOM profile at matched grid point ──
    fn_h  = fullfile(hycom_dir, flist_h(fi).name);
    wt_h  = double(squeeze(ncread(fn_h, 'water_temp', [icol, irow, 1], [1, 1, Inf])));
    ws_h  = double(squeeze(ncread(fn_h, 'salinity',   [icol, irow, 1], [1, 1, Inf])));

    co_dist_km(ip)   = dmin;
    co_dt_hours(ip)  = dt_hrs;
    co_matched(ip)   = true;

    % ── Extract valid Argo levels ──
    pa = pres(:, ip);
    ta = temp_a(:, ip);
    sa = psal_a(:, ip);
    ok_a = ~isnan(pa) & ~isnan(ta) & ~isnan(sa) & pa > 0;
    if sum(ok_a) < 3, continue; end
    pa = pa(ok_a);  ta = ta(ok_a);  sa = sa(ok_a);

    % Deduplicate pressure levels
    [pa, ui] = unique(pa, 'stable');
    ta = ta(ui);  sa = sa(ui);
    if length(pa) < 3, continue; end

    % ── Argo MLD ──
    co_mld_argo(ip) = calc_mld(pa, ta);

    % ── Clean HYCOM depth vector ──
    ok_h = ~isnan(dep_h) & ~isnan(wt_h) & ~isnan(ws_h);
    dh   = dep_h(ok_h);
    wth  = wt_h(ok_h);
    wsh  = ws_h(ok_h);
    [dh, ui_h] = unique(dh, 'stable');
    wth = wth(ui_h);  wsh = wsh(ui_h);
    if length(dh) < 3, continue; end

    % ── HYCOM MLD via the same 0.2 C / 5 m routine used for Argo (like-for-like) ──
    co_mld_hycom(ip) = calc_mld(dh, wth);

    % ── Interpolate both profiles to standard depths ──
    for iz = 1:nz
        z = stat_depths(iz);

        % Argo: interpolate only within data range
        if z >= min(pa) && z <= max(pa)
            co_temp_argo(iz, ip) = interp1(pa, ta, z, 'linear');
            co_psal_argo(iz, ip) = interp1(pa, sa, z, 'linear');
        end

        % HYCOM: interpolate only within data range
        if z >= min(dh) && z <= max(dh)
            co_temp_hycom(iz, ip) = interp1(dh, wth, z, 'linear');
            co_psal_hycom(iz, ip) = interp1(dh, wsh, z, 'linear');
        end
    end
end

fprintf('\n   Matched profiles : %d / %d  (%.1f%%)\n\n', ...
    sum(co_matched), nprof, 100 * sum(co_matched) / nprof);

%% ═══════════════════════════════════════
%  STEP 4 – COMPUTE STATISTICS
%% ═══════════════════════════════════════
fprintf('[4/5] Computing validation statistics...\n');

stat_names = {'Bias','RMSE','MAE','R','Skill_Score','N','StdObs','StdMod'};
row_names  = arrayfun(@(x) sprintf('%dm', x), stat_depths, 'UniformOutput', false);
stats_T = array2table(NaN(nz, 8), 'VariableNames', stat_names, 'RowNames', row_names);
stats_S = array2table(NaN(nz, 8), 'VariableNames', stat_names, 'RowNames', row_names);

for iz = 1:nz
    obs = co_temp_argo(iz, :);
    mdl = co_temp_hycom(iz, :);
    [st, sk] = compute_stats(obs, mdl);
    ok = ~isnan(obs) & ~isnan(mdl);
    stats_T{iz, :} = [st.bias, st.rmse, st.mae, st.r, sk, st.n, ...
                      std(obs(ok), 'omitnan'), std(mdl(ok), 'omitnan')];

    obs = co_psal_argo(iz, :);
    mdl = co_psal_hycom(iz, :);
    [st, sk] = compute_stats(obs, mdl);
    ok = ~isnan(obs) & ~isnan(mdl);
    stats_S{iz, :} = [st.bias, st.rmse, st.mae, st.r, sk, st.n, ...
                      std(obs(ok), 'omitnan'), std(mdl(ok), 'omitnan')];
end

% ── MLD statistics ──
ok_mld = co_matched & ~isnan(co_mld_argo) & ~isnan(co_mld_hycom);
[mld_st, mld_ss] = compute_stats(co_mld_argo(ok_mld)', co_mld_hycom(ok_mld)');

% ── Surface (5 m proxy) statistics ──
[surf_T_st, ~] = compute_stats(co_temp_argo(1, :), co_temp_hycom(1, :));
[surf_S_st, ~] = compute_stats(co_psal_argo(1, :), co_psal_hycom(1, :));

%% ── Print summary to console ──
fprintf('\n================================================================\n');
fprintf('  TEMPERATURE VALIDATION  (HYCOM - Argo)  [degC]\n');
fprintf('================================================================\n');
fprintf('  %-7s %8s %8s %8s %8s %8s %6s\n', 'Depth','Bias','RMSE','MAE','R','Skill','N');
fprintf('  %s\n', repmat('-', 1, 63));
for iz = 1:nz
    v = stats_T{iz, :};
    fprintf('  %-7s %8.3f %8.3f %8.3f %8.3f %8.3f %6d\n', ...
        row_names{iz}, v(1), v(2), v(3), v(4), v(5), v(6));
end

fprintf('\n================================================================\n');
fprintf('  SALINITY VALIDATION  (HYCOM - Argo)  [psu]\n');
fprintf('================================================================\n');
fprintf('  %-7s %8s %8s %8s %8s %8s %6s\n', 'Depth','Bias','RMSE','MAE','R','Skill','N');
fprintf('  %s\n', repmat('-', 1, 63));
for iz = 1:nz
    v = stats_S{iz, :};
    fprintf('  %-7s %8.3f %8.3f %8.3f %8.3f %8.3f %6d\n', ...
        row_names{iz}, v(1), v(2), v(3), v(4), v(5), v(6));
end

fprintf('\n================================================================\n');
fprintf('  MLD  Bias=%.2fm  RMSE=%.2fm  MAE=%.2fm  R=%.3f  Skill=%.3f  N=%d\n', ...
    mld_st.bias, mld_st.rmse, mld_st.mae, mld_st.r, mld_ss, mld_st.n);
fprintf('================================================================\n\n');

%% ── Save CSVs ──
writetable(stats_T, fullfile(out_dir, 'stats_temperature.csv'), 'WriteRowNames', true);
writetable(stats_S, fullfile(out_dir, 'stats_salinity.csv'),    'WriteRowNames', true);

%% ══════════════════════════════════════════════════════════
%  STEP 5 – WRITE FULL TEXT REPORT
%% ══════════════════════════════════════════════════════════
fprintf('[5/5] Writing full report...\n');
fid = fopen(fullfile(out_dir, 'validation_full_report.txt'), 'w');

fprintf(fid,'================================================================\n');
fprintf(fid,'        HYCOM vs ARGO FULL VALIDATION REPORT\n');
fprintf(fid,'        Cyclone %s - %s\n', cfg.prefix, cfg.label);
fprintf(fid,'================================================================\n\n');
fprintf(fid,'Run date/time   : %s\n',   datestr(now));
fprintf(fid,'Argo source dir : %s\n',   argo_dir);
fprintf(fid,'HYCOM directory : %s\n',   hycom_dir);
fprintf(fid,'Output folder   : %s\n\n', out_dir);

fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'COLOCATION SETTINGS\n');
fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'  Max spatial distance : %d km\n',    max_dist_km);
fprintf(fid,'  Max time difference  : %d hours\n', max_dt_hours);
fprintf(fid,'  MLD method           : 0.2 deg C threshold from 5m reference (Argo and HYCOM)\n');
fprintf(fid,'  Argo QC filter       : flags 1 & 2 only (good / probably good)\n');
fprintf(fid,'  Adjusted fields      : yes (fallback to raw if unavailable)\n\n');

fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'DATA SUMMARY\n');
fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'  Argo daily files read        : %d\n',  numel(flist_a));
fprintf(fid,'  Total Argo profiles (post QC): %d\n',  nprof);
fprintf(fid,'  Profiles matched to HYCOM    : %d\n',  sum(co_matched));
fprintf(fid,'  Match rate                   : %.1f%%\n', 100*sum(co_matched)/nprof);
fprintf(fid,'  Argo date range  : %s  to  %s\n', ...
    datestr(min(argo_dnum),'yyyy-mmm-dd'), datestr(max(argo_dnum),'yyyy-mmm-dd'));
fprintf(fid,'  HYCOM date range : %s  to  %s\n', ...
    datestr(hycom_dnum(1),'yyyy-mmm-dd HH:MM'), datestr(hycom_dnum(end),'yyyy-mmm-dd HH:MM'));
fprintf(fid,'  HYCOM files used : %d\n', nf);
fprintf(fid,'  Argo lat range   : %.3f N  to  %.3f N\n', min(lat_a), max(lat_a));
fprintf(fid,'  Argo lon range   : %.3f E  to  %.3f E\n\n', min(lon_a), max(lon_a));

fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'COLOCATION QUALITY  (matched profiles only)\n');
fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'  Distance (km) : mean=%.2f  median=%.2f  std=%.2f  min=%.2f  max=%.2f\n', ...
    mean(co_dist_km(co_matched),'omitnan'), median(co_dist_km(co_matched),'omitnan'), ...
    std(co_dist_km(co_matched),'omitnan'),  min(co_dist_km(co_matched)), ...
    max(co_dist_km(co_matched)));
fprintf(fid,'  Time diff (h) : mean=%.2f  median=%.2f  std=%.2f  min=%.2f  max=%.2f\n\n', ...
    mean(co_dt_hours(co_matched),'omitnan'), median(co_dt_hours(co_matched),'omitnan'), ...
    std(co_dt_hours(co_matched),'omitnan'),  min(co_dt_hours(co_matched)), ...
    max(co_dt_hours(co_matched)));

fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'TEMPERATURE VALIDATION  (HYCOM - Argo)  [deg C]\n');
fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'%-8s %10s %10s %10s %10s %10s %8s %10s %10s\n', ...
    'Depth(m)','Bias','RMSE','MAE','R','Skill','N','StdArgo','StdHYCOM');
fprintf(fid,'%s\n', repmat('-', 1, 98));
for iz = 1:nz
    v = stats_T{iz, :};
    fprintf(fid,'%-8s %10.4f %10.4f %10.4f %10.4f %10.4f %8d %10.4f %10.4f\n', ...
        row_names{iz}, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8));
end

fprintf(fid,'\n--------------------------------------------------------------\n');
fprintf(fid,'SALINITY VALIDATION  (HYCOM - Argo)  [psu]\n');
fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'%-8s %10s %10s %10s %10s %10s %8s %10s %10s\n', ...
    'Depth(m)','Bias','RMSE','MAE','R','Skill','N','StdArgo','StdHYCOM');
fprintf(fid,'%s\n', repmat('-', 1, 98));
for iz = 1:nz
    v = stats_S{iz, :};
    fprintf(fid,'%-8s %10.4f %10.4f %10.4f %10.4f %10.4f %8d %10.4f %10.4f\n', ...
        row_names{iz}, v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8));
end

fprintf(fid,'\n--------------------------------------------------------------\n');
fprintf(fid,'MIXED LAYER DEPTH  (MLD)  [m]\n');
fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'  N matched          : %d\n',     mld_st.n);
fprintf(fid,'  Bias (HYCOM-Argo)  : %+.4f m\n', mld_st.bias);
fprintf(fid,'  RMSE               : %.4f m\n',  mld_st.rmse);
fprintf(fid,'  MAE                : %.4f m\n',  mld_st.mae);
fprintf(fid,'  Correlation  (R)   : %.4f\n',    mld_st.r);
fprintf(fid,'  Murphy Skill Score : %.4f\n',    mld_ss);
fprintf(fid,'  Argo  MLD - mean=%.2f  std=%.2f  min=%.2f  max=%.2f m\n', ...
    mean(co_mld_argo(ok_mld),'omitnan'), std(co_mld_argo(ok_mld),'omitnan'), ...
    min(co_mld_argo(ok_mld)), max(co_mld_argo(ok_mld)));
fprintf(fid,'  HYCOM MLD - mean=%.2f  std=%.2f  min=%.2f  max=%.2f m\n\n', ...
    mean(co_mld_hycom(ok_mld),'omitnan'), std(co_mld_hycom(ok_mld),'omitnan'), ...
    min(co_mld_hycom(ok_mld)), max(co_mld_hycom(ok_mld)));

fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'NEAR-SURFACE (5 m) STATISTICS\n');
fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'  Temperature : Bias=%+.4f C   RMSE=%.4f  MAE=%.4f  R=%.4f  N=%d\n', ...
    surf_T_st.bias, surf_T_st.rmse, surf_T_st.mae, surf_T_st.r, surf_T_st.n);
fprintf(fid,'  Salinity    : Bias=%+.4f psu  RMSE=%.4f  MAE=%.4f  R=%.4f  N=%d\n\n', ...
    surf_S_st.bias, surf_S_st.rmse, surf_S_st.mae, surf_S_st.r, surf_S_st.n);

fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'PER-PROFILE COLOCATION TABLE  (matched profiles only)\n');
fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'%-5s %-12s %8s %9s %9s %7s %8s %8s %8s %8s %10s %10s\n', ...
    'Idx','Date','Lat','Lon','Dist_km','dT_hrs', ...
    'T5_Argo','T5_HYCOM','S5_Argo','S5_HYCOM','MLD_Argo','MLD_HYCOM');
fprintf(fid,'%s\n', repmat('-', 1, 118));
matched_idx = find(co_matched);
for k = 1:length(matched_idx)
    ip = matched_idx(k);
    fprintf(fid,'%-5d %-12s %8.3f %9.3f %9.2f %7.2f %8.3f %8.3f %8.3f %8.3f %10.2f %10.2f\n', ...
        ip, datestr(argo_dnum(ip), 'yyyy-mm-dd'), ...
        lat_a(ip), lon_a(ip), co_dist_km(ip), co_dt_hours(ip), ...
        co_temp_argo(1, ip),  co_temp_hycom(1, ip), ...
        co_psal_argo(1, ip),  co_psal_hycom(1, ip), ...
        co_mld_argo(ip),      co_mld_hycom(ip));
end

fprintf(fid,'\n--------------------------------------------------------------\n');
fprintf(fid,'OUTPUT FILES\n');
fprintf(fid,'--------------------------------------------------------------\n');
fprintf(fid,'  %-35s : %s\n', 'validation_full_report.txt', 'This report');
fprintf(fid,'  %-35s : %s\n', 'stats_temperature.csv',      'Temperature stats per depth (CSV)');
fprintf(fid,'  %-35s : %s\n', 'stats_salinity.csv',         'Salinity stats per depth (CSV)');

fprintf(fid,'\n================================================================\n');
fprintf(fid,'                    END OF REPORT\n');
fprintf(fid,'================================================================\n');
fclose(fid);

fprintf('\n All outputs saved to: %s\n', out_dir);
fprintf('  Report : validation_full_report.txt\n');
fprintf('  CSV    : stats_temperature.csv, stats_salinity.csv\n\n');

%% ════════════════════════════════════════════
%                    HELPER FUNCTIONS
%% ═══════════════════════════════════════════

function [st, skill] = compute_stats(obs, mdl)
% Returns bias, RMSE, MAE, R, and Murphy Skill Score.
    obs = obs(:);  mdl = mdl(:);
    ok  = ~isnan(obs) & ~isnan(mdl);
    st.n = sum(ok);
    if st.n < 3
        st.bias = NaN;  st.rmse = NaN;  st.mae = NaN;
        st.r    = NaN;  skill   = NaN;
        return;
    end
    o = obs(ok);  m = mdl(ok);
    st.bias = mean(m - o);
    st.rmse = sqrt(mean((m - o).^2));
    st.mae  = mean(abs(m - o));
    cc      = corrcoef(o, m);
    st.r    = cc(1, 2);
    mse_clim = mean((o - mean(o)).^2);
    mse_mod  = mean((m - o).^2);
    skill    = 1 - mse_mod / (mse_clim + eps);
end

function mld = calc_mld(pres, temp)
% Mixed Layer Depth via 0.2 deg C threshold (de Boyer Montegut et al. 2004),
% referenced to 5 m (or the shallowest sample if deeper than 5 m).
    mld = NaN;
    if length(pres) < 3, return; end
    [pres, si] = sort(pres);
    temp = temp(si);
    [pres, ui] = unique(pres, 'stable');
    temp = temp(ui);
    if length(pres) < 3, return; end
    z_ref = max(5, pres(1));
    if z_ref > pres(end)
        mld = pres(end);
        return;
    end
    t_sfc = interp1(pres, temp, z_ref, 'linear');
    idx   = find(temp < t_sfc - 0.2, 1, 'first');
    if isempty(idx) || idx == 1
        mld = pres(end);   % well-mixed throughout
    else
        mld = interp1(temp(idx-1:idx), pres(idx-1:idx), t_sfc - 0.2, 'linear');
    end
end
