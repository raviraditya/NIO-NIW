% Sea-level anomaly context along each track. Reads the daily gridded SLA
% product and writes the per-cyclone summaries and Figure S9.

clc; clear; close all; warning('off','all');

DATA_DIR = 'path/to/data';        % folder with IMD, Final/All/Combined etc.
SLA_DIR  = 'path/to/sla';         % folder with cropped_ssha_YYYY.nc


%% STORM LIST
STORM_LIST = {'KYARR','AMPHAN','FANI','TAUKTAE'};   % runs all sequentially

for iStorm = 1:numel(STORM_LIST)
STORM_ID = STORM_LIST{iStorm};
close all;

%% CONFIG
cfg = get_sla_config(STORM_ID);

fprintf('\n============================================================\n');
fprintf(' SLA / EDDY CONTEXT v1  -  STORM: %s\n', cfg.name);
fprintf(' year_ref = %d   |   R_box = %d km   |   SLA file: %s\n', ...
    cfg.year_ref, cfg.R_box_km, cfg.sla_file);
fprintf('============================================================\n\n');

% Anticyclonic / cyclonic SLA thresholds (m).  +0.05 / -0.05 m follows
% the conservative thresholds used in mesoscale eddy literature for the
% NIO (Chacko 2021, Peng 2023).
SLA_ANTI =  0.05;
SLA_CYC  = -0.05;

% Box half-width (deg) used for track-centred mean SLA. Set so 2*R_box
% covers ~R_max for typical NIO TCs (R_max ~150-200 km ~ 1.5 deg).
R_box_deg = cfg.R_box_km / 111;

%% READ TRACK
T   = readtable(cfg.track_path);
tt  = datetime(T.year, T.month, T.day, T.hour, 0, 0);
latC = T.lat;
lonC = T.lon;
nt   = numel(tt);

%% READ SLA
lon_s  = ncread(cfg.sla_file, 'longitude');
lat_s  = ncread(cfg.sla_file, 'latitude');
t_days = ncread(cfg.sla_file, 'time');
sla    = ncread(cfg.sla_file, 'sla');   % lon x lat x time

t_ref  = datetime(cfg.year_ref, 1, 1);
t_sla  = t_ref + days(t_days);

% Make sure lon/lat are increasing
if lon_s(2) < lon_s(1)
    lon_s = flipud(lon_s(:));
    sla   = flip(sla, 1);
end
if lat_s(2) < lat_s(1)
    lat_s = flipud(lat_s(:));
    sla   = flip(sla, 2);
end

[LON_S, LAT_S] = meshgrid(lon_s, lat_s);   % [nlat x nlon]

%% PER-TIMESTEP DIAGNOSTICS
SLA_at_centre   = nan(nt,1);   % SLA at storm centre (m)
SLA_box_mean    = nan(nt,1);   % mean SLA inside R_box (m)
SLA_box_min     = nan(nt,1);
SLA_box_max     = nan(nt,1);
Frac_anti       = nan(nt,1);   % fraction of box area with SLA > +0.05
Frac_cyc        = nan(nt,1);   % fraction of box area with SLA < -0.05
EddyClass       = strings(nt,1);   % 'anticyclonic' | 'cyclonic' | 'neutral'

for it = 1:nt
    % nearest SLA day
    [~, kd] = min(abs(days(t_sla - tt(it))));
    F = sla(:,:,kd);     % nlon x nlat
    F = permute(F, [2 1]);   % nlat x nlon  to match LON_S/LAT_S

    % point sample at storm centre
    SLA_at_centre(it) = interp2(LON_S, LAT_S, F, lonC(it), latC(it), 'linear', NaN);

    % box mask
    inbox = (LON_S >= lonC(it)-R_box_deg) & (LON_S <= lonC(it)+R_box_deg) & ...
            (LAT_S >= latC(it)-R_box_deg) & (LAT_S <= latC(it)+R_box_deg);
    vals = F(inbox);
    vals = vals(~isnan(vals));
    if ~isempty(vals)
        SLA_box_mean(it) = mean(vals);
        SLA_box_min(it)  = min(vals);
        SLA_box_max(it)  = max(vals);
        Frac_anti(it)    = sum(vals >  SLA_ANTI) / numel(vals);
        Frac_cyc(it)     = sum(vals <  SLA_CYC ) / numel(vals);
    end

    % classification
    if SLA_box_mean(it) >  SLA_ANTI
        EddyClass(it) = "anticyclonic";
    elseif SLA_box_mean(it) <  SLA_CYC
        EddyClass(it) = "cyclonic";
    else
        EddyClass(it) = "neutral";
    end
end

%% STORM-PERIOD SUMMARY
S.storm              = cfg.name;
S.year               = cfg.year_ref;
S.R_box_km           = cfg.R_box_km;
S.SLA_anti_thresh    = SLA_ANTI;
S.SLA_cyc_thresh     = SLA_CYC;

S.SLA_centre_mean    = mean(SLA_at_centre, 'omitnan');
S.SLA_centre_median  = median(SLA_at_centre, 'omitnan');
S.SLA_centre_max     = max(SLA_at_centre, [], 'omitnan');
S.SLA_centre_min     = min(SLA_at_centre, [], 'omitnan');

S.SLA_box_mean       = mean(SLA_box_mean, 'omitnan');
S.SLA_box_median     = median(SLA_box_mean, 'omitnan');
S.SLA_box_max        = max(SLA_box_max, [], 'omitnan');
S.SLA_box_min        = min(SLA_box_min, [], 'omitnan');

S.Frac_anti_mean     = mean(Frac_anti, 'omitnan');
S.Frac_cyc_mean      = mean(Frac_cyc,  'omitnan');

S.frac_anticyclonic  = sum(EddyClass == "anticyclonic") / nt;
S.frac_cyclonic      = sum(EddyClass == "cyclonic")     / nt;
S.frac_neutral       = sum(EddyClass == "neutral")      / nt;

% Peak-NIKE snapshot (read from existing supercharged xlsx if present)
peak_idx = NaN;  peak_time = NaT;  SLA_at_peak = NaN;  EddyClass_at_peak = "n/a";
if exist(cfg.supercharged_xlsx, 'file')
    try
        Tnike = readtable(cfg.supercharged_xlsx, 'Sheet', 'Sheet1');
        % column names from Main_combined_new v6 (Sheet1)
        cand = {'NIKE_Coherent','NIKE_coherent_kJm2','NIKE_coherent_Jm2','NIKE_coh_col','NIKE_coh_kJm2'};
        col = '';
        for c = 1:numel(cand)
            if any(strcmpi(Tnike.Properties.VariableNames, cand{c}))
                col = cand{c}; break;
            end
        end
        if ~isempty(col)
            nike = Tnike.(col);
            [~, peak_idx] = max(nike);
            % find matching time in track
            if any(strcmpi(Tnike.Properties.VariableNames, 'time'))
                t_xlsx = Tnike.time;
                [~, peak_idx_track] = min(abs(tt - t_xlsx(peak_idx)));
            else
                peak_idx_track = peak_idx;
            end
            peak_time         = tt(peak_idx_track);
            SLA_at_peak       = SLA_box_mean(peak_idx_track);
            EddyClass_at_peak = EddyClass(peak_idx_track);
        end
    catch ME
        warning('Could not read peak-NIKE from %s: %s', cfg.supercharged_xlsx, ME.message);
    end
end
S.peak_NIKE_time        = peak_time;
S.SLA_at_peak_NIKE      = SLA_at_peak;
S.EddyClass_at_peak     = EddyClass_at_peak;

%% WRITE OUTPUTS
OUT_DIR = cfg.out_dir;
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

% (a) MAT file
save(fullfile(OUT_DIR, [cfg.name '_SLA.mat']), ...
    'tt','latC','lonC','SLA_at_centre','SLA_box_mean','SLA_box_min','SLA_box_max', ...
    'Frac_anti','Frac_cyc','EddyClass','S','cfg');

% (b) XLSX timeseries
Tout = table(tt, latC, lonC, SLA_at_centre, SLA_box_mean, SLA_box_min, SLA_box_max, ...
    Frac_anti, Frac_cyc, EddyClass, ...
    'VariableNames', {'time','lat','lon','SLA_centre_m','SLA_box_mean_m', ...
                      'SLA_box_min_m','SLA_box_max_m','Frac_anticyc','Frac_cyc','EddyClass'});
writetable(Tout, fullfile(OUT_DIR, [cfg.name '_SLA.xlsx']), 'Sheet','timeseries');

Tsum = struct2table(rmfield(S, intersect({'cfg'}, fieldnames(S))), 'AsArray', true);
writetable(Tsum, fullfile(OUT_DIR, [cfg.name '_SLA.xlsx']), 'Sheet','summary');

%% STORE PER-STORM DATA FOR COMBINED FIG
ALL(iStorm).name        = cfg.name;
ALL(iStorm).tt          = tt;
ALL(iStorm).latC        = latC;
ALL(iStorm).lonC        = lonC;
ALL(iStorm).SLA_centre  = SLA_at_centre;
ALL(iStorm).SLA_box     = SLA_box_mean;
ALL(iStorm).peak_time   = peak_time;
ALL(iStorm).SLA_at_peak = SLA_at_peak;
ALL(iStorm).S           = S;

%% CONSOLE REPORT
fprintf('\n--- %s summary ---\n', cfg.name);
fprintf('   track-mean SLA at centre  : %+.3f m\n', S.SLA_centre_mean);
fprintf('   track-mean box SLA        : %+.3f m\n', S.SLA_box_mean);
fprintf('   anti / cyc / neutral frac : %5.1f / %5.1f / %5.1f %%\n', ...
    100*S.frac_anticyclonic, 100*S.frac_cyclonic, 100*S.frac_neutral);
if ~isnat(peak_time)
    fprintf('   SLA at peak NIKE          : %+.3f m  (%s)\n', SLA_at_peak, EddyClass_at_peak);
end
fprintf('\nWrote outputs to: %s\n', OUT_DIR);

end  % end storm loop
%  COMBINED FIGURE — all four storms, 2 rows x 4 columns
%  Row (a): SLA box-mean time series with peak-NIKE marker
%  Row (b): SLA at storm centre vs box-mean
%  Each subplot has its own x and y range (matches main figures style).
COMBINED_OUT = fullfile(pwd, 'SLA_EDDY_OUTPUT');
if ~exist(COMBINED_OUT,'dir'), mkdir(COMBINED_OUT); end

storm_order = {'KYARR','AMPHAN','FANI','TAUKTAE'};
storm_basin = {'AS','BoB','BoB','AS'};

% reorder ALL to match canonical paper order
ord = zeros(1,4);
for k = 1:4
    for j = 1:numel(ALL)
        if strcmp(ALL(j).name, storm_order{k}), ord(k) = j; break; end
    end
end
ALL = ALL(ord);

SLA_ANTI =  0.05;
SLA_CYC  = -0.05;

fC = figure('Color','w','Position',[40 40 1500 700]);
tl = tiledlayout(fC, 2, 4, 'Padding','compact', 'TileSpacing','compact');

for k = 1:4
    D = ALL(k);

    % --- Row (a): box-mean SLA with anti/cyc shading ---
    nexttile(k);
    ax = gca; hold(ax,'on');
    % shaded bands for anti / cyc thresholds
    yl_pad = max(abs([D.SLA_box; D.SLA_centre]),[],'omitnan') * 1.15;
    if isempty(yl_pad) || isnan(yl_pad), yl_pad = 0.25; end
    fill([D.tt(1) D.tt(end) D.tt(end) D.tt(1)], ...
         [SLA_ANTI SLA_ANTI yl_pad yl_pad], [1.0 0.85 0.85], ...
         'EdgeColor','none','FaceAlpha',0.5);
    fill([D.tt(1) D.tt(end) D.tt(end) D.tt(1)], ...
         [-yl_pad -yl_pad SLA_CYC SLA_CYC], [0.85 0.85 1.0], ...
         'EdgeColor','none','FaceAlpha',0.5);

    plot(D.tt, D.SLA_box, 'k-', 'LineWidth',1.4);
    yline(SLA_ANTI,'--','Color',[0.5 0 0]);
    yline(SLA_CYC ,'--','Color',[0 0 0.5]);
    yline(0,':','Color',[0.4 0.4 0.4]);

    if ~isnat(D.peak_time)
        xline(D.peak_time, 'k-', 'LineWidth',1.2);
        % annotate SLA at peak NIKE
        text(0.02, 0.92, sprintf('SLA@peak = %+.3f m', D.SLA_at_peak), ...
            'Units','normalized','FontSize',9,'FontWeight','bold', ...
            'BackgroundColor',[1 1 1 0.7]);
    end

    ylim([-yl_pad yl_pad]);
    xlim([D.tt(1) D.tt(end)]);
    grid on; box on;
    title(sprintf('%s (%s)', storm_order{k}, storm_basin{k}), 'FontSize',11);
    if k == 1, ylabel('SLA box-mean (m)'); end
    set(ax,'FontSize',9);

    % --- Row (b): SLA at centre vs box mean ---
    nexttile(4 + k);
    ax = gca; hold(ax,'on');
    plot(D.tt, D.SLA_centre, 'b-', 'LineWidth',1.2);
    plot(D.tt, D.SLA_box,    'r-', 'LineWidth',1.2);
    yline(0,':','Color',[0.4 0.4 0.4]);
    if ~isnat(D.peak_time), xline(D.peak_time, 'k-', 'LineWidth',1.0); end

    yl_pad2 = max(abs([D.SLA_box; D.SLA_centre]),[],'omitnan')*1.15;
    if isempty(yl_pad2) || isnan(yl_pad2), yl_pad2 = 0.25; end
    ylim([-yl_pad2 yl_pad2]);
    xlim([D.tt(1) D.tt(end)]);
    grid on; box on;
    if k == 1
        ylabel('SLA (m)');
        legend({'at centre','box mean (R=150 km)'}, 'Location','best','FontSize',8);
    end
    xlabel('Date (UTC)');
    set(ax,'FontSize',9);

    % print storm-period stats
    frac_anti = D.S.frac_anticyclonic * 100;
    frac_cyc  = D.S.frac_cyclonic     * 100;
    text(0.02, 0.10, sprintf('anti %.0f%% / cyc %.0f%%', frac_anti, frac_cyc), ...
        'Units','normalized','FontSize',8, ...
        'BackgroundColor',[1 1 1 0.7]);
end

% row labels using annotation
annotation(fC,'textbox',[0.005 0.74 0.02 0.04], 'String','(a)', ...
    'EdgeColor','none','FontWeight','bold','FontSize',12);
annotation(fC,'textbox',[0.005 0.30 0.02 0.04], 'String','(b)', ...
    'EdgeColor','none','FontWeight','bold','FontSize',12);

%sgtitle('Mesoscale eddy context: SLA along storm tracks', 'FontWeight','bold','FontSize',12);

exportgraphics(fC, fullfile(COMBINED_OUT,'SLA_Combined_AllStorms.png'), 'Resolution',300);
fprintf('\nCombined figure written: %s\n', fullfile(COMBINED_OUT,'SLA_Combined_AllStorms.png'));

fprintf('\n============================================================\n');
fprintf(' ALL STORMS COMPLETE.\n');
fprintf('============================================================\n');
%  CONFIG FUNCTION
function cfg = get_sla_config(storm_id)
    base_track   = [fullfile(DATA_DIR,'IMD') filesep];
    base_sla     = [SLA_DIR filesep];
    base_super   = [fullfile(DATA_DIR,'Final','All','Combined') filesep];
    cfg.out_dir  = fullfile(pwd, 'SLA_EDDY_OUTPUT');
    cfg.R_box_km = 150;   % match R_NI in Main_combined_new.m

    switch upper(storm_id)
        case 'KYARR'
            cfg.name              = 'KYARR';
            cfg.year_ref          = 2019;
            cfg.track_path        = [base_track 'Kyarr.csv'];
            cfg.sla_file          = [base_sla   'cropped_ssha_2019.nc'];
            cfg.supercharged_xlsx = [base_super 'KYARR_V16_SUPERCHARGED.xlsx'];

        case 'AMPHAN'
            cfg.name              = 'AMPHAN';
            cfg.year_ref          = 2020;
            cfg.track_path        = [base_track 'Amphan.csv'];
            cfg.sla_file          = [base_sla   'cropped_ssha_2020.nc'];
            cfg.supercharged_xlsx = [base_super 'AMPHAN_V16_SUPERCHARGED.xlsx'];

        case 'FANI'
            cfg.name              = 'FANI';
            cfg.year_ref          = 2019;
            cfg.track_path        = [base_track 'Fani.csv'];
            cfg.sla_file          = [base_sla   'cropped_ssha_2019.nc'];
            cfg.supercharged_xlsx = [base_super 'FANI_V16_SUPERCHARGED.xlsx'];

        case 'TAUKTAE'
            cfg.name              = 'TAUKTAE';
            cfg.year_ref          = 2021;
            cfg.track_path        = [base_track 'Tauktae.csv'];
            cfg.sla_file          = [base_sla   'cropped_ssha_2021.nc'];
            cfg.supercharged_xlsx = [base_super 'TAUKTAE_V16_SUPERCHARGED.xlsx'];

        otherwise
            error('Unknown STORM_ID: %s. Valid: KYARR | AMPHAN | FANI | TAUKTAE', storm_id);
    end
end
%  Local red-white-blue colormap (avoids needing colormap toolbox)
function cmap = redblue_local(n)
    if nargin < 1, n = 256; end
    half = floor(n/2);
    r = [linspace(0.03,1,half) linspace(1,0.7,n-half)]';
    g = [linspace(0.19,1,half) linspace(1,0.0,n-half)]';
    b = [linspace(0.42,1,half) linspace(1,0.15,n-half)]';
    cmap = [r g b];
end
