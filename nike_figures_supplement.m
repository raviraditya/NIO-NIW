% Supporting Information figures S1-S8. Reads *_V17.mat.

clc; clear; close all;

DATA_DIR = 'path/to/data';        % folder with HYCOM, IMD, Stress, ARGO etc.
GSW_PATH = 'path/to/GSW-Matlab';  % Gibbs SeaWater toolbox

%% USER CONFIG
BASE_DIR = fullfile(DATA_DIR,'Final','All','Combined');
OUT_DIR  = fullfile(pwd, 'PAPER1_FIGURES_SUPP');

% HYCOM paths per storm.
% Fig S7 (R_NI sensitivity) needs velocity (UV3Z).
% Fig S1 (MLD criteria)     needs temperature/salinity (TS3Z).
HYCOM_UV = struct( ...
    'KYARR',   fullfile(DATA_DIR,'HYCOM','Kyarr_UV3Z'),   ...
    'AMPHAN',  fullfile(DATA_DIR,'HYCOM','Amphan_UV3Z'),  ...
    'FANI',    fullfile(DATA_DIR,'HYCOM','Fani_UV3Z'),    ...
    'TAUKTAE', fullfile(DATA_DIR,'HYCOM','Tauktae_UV3Z'));

HYCOM_TS = struct( ...
    'KYARR',   fullfile(DATA_DIR,'HYCOM','Kyarr_TS3Z'),   ...
    'AMPHAN',  fullfile(DATA_DIR,'HYCOM','Amphan_TS3Z'),  ...
    'FANI',    fullfile(DATA_DIR,'HYCOM','Fani_TS3Z'),    ...
    'TAUKTAE', fullfile(DATA_DIR,'HYCOM','Tauktae_TS3Z'));

% Merged HYCOM-vs-Argo validation report used by Fig S2
MERGED_VALIDATION_REPORT = fullfile(DATA_DIR,'HYCOM','Validation','argo_hycom_report.txt');


% GSW toolbox (for the Fig S1 density calculation)
GSW_PATH = GSW_PATH;
if exist(GSW_PATH,'dir'), addpath(genpath(GSW_PATH)); end

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

STORM_NAMES = {'KYARR','AMPHAN','FANI','TAUKTAE'};
STORM_BASIN = {'AS',   'BoB',   'BoB', 'AS'};

% PAPER figure numbers. Figs 1 and 7 are slow (HYCOM re-reads).
FIGURES_TO_MAKE = [1 2 3 4 5 6 7 8];   % S9 comes from the separate SLA script

n_cyc = numel(STORM_NAMES);

%% GLOBAL STYLE
% [S9] Stop the interactive toolbar being captured in exported PNGs.
set(groot, ...
    'defaultFigureToolBar','none', ...
    'defaultFigureMenuBar','none', ...
    'defaultFigureColor','w', ...
    'defaultAxesFontName','Helvetica', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesLineWidth',1.2, ...
    'defaultLineLineWidth',1.7, ...
    'defaultAxesTitleFontWeight','bold');

c_grid = [0.85 0.85 0.85];

beautify = @() set(gca, ...
    'TickDir','out','Box','on','Layer','top', ...
    'FontSize',10,'GridAlpha',0.20, ...
    'XGrid','on','YGrid','on','GridColor',c_grid);

Omega  = 7.2921e-5;
dt_hr  = 3;
Fs_cpd = 24/dt_hr;
f_M2   = 1.9323;
f_K1   = 1.0027;
rho0   = 1025;
Hmax   = 700;

fprintf('\n================================================================\n');
fprintf(' PAPER 1 - SUPPORTING INFORMATION FIGURES S1..S8\n');
fprintf(' Output : %s\n', OUT_DIR);
fprintf('================================================================\n');
fprintf('\n Figure map (paper number -> file -> was called in draft)\n');
fprintf(' ---------------------------------------------------------------\n');
map_note = {'S1','draft S3'; 'S2','draft S2'; 'S3','draft S5'; 'S4','draft S6'; ...
            'S5','draft S7'; 'S6','draft S8'; 'S7','draft S9'; 'S8','draft S10'};
for ii = 1:size(map_note,1)
    fprintf('  %-4s  %-42s  (%s)\n', map_note{ii,1}, figFileName(map_note{ii,1}), map_note{ii,2});
end
fprintf(' ---------------------------------------------------------------\n');

%% LOAD COMMON DATA ONCE (v7)
%  Reads *_V17.mat directly. Builds a T table with the field names the
D = struct();

for i = 1:n_cyc
    storm = STORM_NAMES{i};
    fprintf('Loading %s ...\n', storm);

    mat_main = fullfile(BASE_DIR, [storm, '_V17.mat']);
    if ~isfile(mat_main), error('Missing %s', mat_main); end
    M = load(mat_main);

    time = M.time_tc(:);
    nT   = numel(time);

    T = table();
    T.time_tc        = time;
    T.NIKE_Broadband = getcol(M,'NIKE_broadband_Jm2',nT);
    T.NIKE_Coherent  = getcol(M,'NIKE_coherent_Jm2',nT);
    T.MLD            = getcol(M,'MLD_m',nT);
    T.lat            = getcol(M,'lat_tc',nT);
    T.lon            = getcol(M,'lon_tc',nT);

    D(i).name  = storm;
    D(i).basin = STORM_BASIN{i};
    D(i).time  = time;
    D(i).T     = T;
    D(i).M     = M;
    D(i).MLD   = T.MLD;
    D(i).lat   = T.lat;
    D(i).lon   = T.lon;
    D(i).lat_mean = mean(D(i).lat,'omitnan');
    if ~isfinite(D(i).lat_mean), D(i).lat_mean = 16; end
    D(i).f_inertial = abs(2*Omega*sind(D(i).lat_mean)) * 86400 / (2*pi);

    if isfield(M,'z'), D(i).z = M.z(:); else, D(i).z = []; end

    D(i).NIKE_coh_z   = fldz(M,'NIKE_coherent_z_save');
    D(i).NIKE_brd_z   = fldz(M,'NIKE_broadband_z_save');
    D(i).NIKE_coh_z_R = fldz(M,'NIKE_coherent_z_right');
    D(i).NIKE_coh_z_L = fldz(M,'NIKE_coherent_z_left');
    D(i).N2           = fld(M,'N2_save');
    D(i).S2           = fld(M,'S2_save');
    D(i).MLD_main     = colv(fld(M,'MLD_m'));

    % surface velocity time series (window centre column)
    if isfield(M,'u_surf_save') && isfield(M,'v_surf_save')
        i0c = ceil(size(M.u_surf_save,2)/2);
        D(i).u_surf = M.u_surf_save(:, i0c);
        D(i).v_surf = M.v_surf_save(:, i0c);
    else
        D(i).u_surf = []; D(i).v_surf = [];
    end

    D(i).nike_col = T.NIKE_Coherent;
end

%% START QC FILE
fprintf('================================================================\n');
fprintf(' PAPER 1 - SUPPORTING INFORMATION FIGURES QC SUMMARY\n');
fprintf(' Figure numbers follow the Supporting Information\n');
fprintf(' Generated : %s\n', datestr(now)); %#ok<TNOW1,DATST>
fprintf('================================================================\n\n');

%% FIGURES
fig_ok = false(1,8);

if ismember(1, FIGURES_TO_MAKE)      % was make_figS3
    try
        fig_ok(1) = make_figS1_MLDMethods(D, n_cyc, OUT_DIR, beautify, HYCOM_TS, Hmax);
    catch ME
        report_fig_error('S1 (MLD criterion comparison)', ME);
    end
end

if ismember(2, FIGURES_TO_MAKE)      % was make_figS2
    try
        fig_ok(2) = make_figS2_HycomArgo(STORM_NAMES, MERGED_VALIDATION_REPORT, OUT_DIR, beautify);
    catch ME
        report_fig_error('S2 (HYCOM-Argo validation)', ME);
    end
end

if ismember(3, FIGURES_TO_MAKE)      % was make_figS5
    try
        fig_ok(3) = make_figS3_RotarySpectra(D, n_cyc, OUT_DIR, beautify, Fs_cpd, f_M2, f_K1);
    catch ME
        report_fig_error('S3 (rotary spectra)', ME);
    end
end

if ismember(4, FIGURES_TO_MAKE)      % was make_figS6
    try
        fig_ok(4) = make_figS4_Spectrograms(D, n_cyc, OUT_DIR, beautify, Fs_cpd, f_M2);
    catch ME
        report_fig_error('S4 (inertial spectrograms)', ME);
    end
end

if ismember(5, FIGURES_TO_MAKE)      % was make_figS7
    try
        fig_ok(5) = make_figS5_VertProfiles(D, n_cyc, OUT_DIR, beautify, rho0);
    catch ME
        report_fig_error('S5 (vertical profile snapshots)', ME);
    end
end

if ismember(6, FIGURES_TO_MAKE)      % was make_figS8
    try
        fig_ok(6) = make_figS6_RLDepth(D, n_cyc, OUT_DIR, beautify);
    catch ME
        report_fig_error('S6 (R/L depth-resolved)', ME);
    end
end

if ismember(7, FIGURES_TO_MAKE)      % was make_figS9
    try
        fig_ok(7) = make_figS7_RNISensitivity(D, n_cyc, OUT_DIR, beautify, HYCOM_UV, Hmax, rho0);
    catch ME
        report_fig_error('S7 (R_NI sensitivity)', ME);
    end
end

if ismember(8, FIGURES_TO_MAKE)      % was make_figS10
    try
        fig_ok(8) = make_figS8_PhaseSpace(D, n_cyc, OUT_DIR, beautify);
    catch ME
        report_fig_error('S8 (N2-S2 phase space)', ME);
    end
end


%% Save check
% Names come from the SAME helper the writers use, so this check can
% never look for a file the writers never intended to create.
fprintf('\n---------------------------------------------------------------\n');
fprintf(' Save check - output file names\n');
fprintf('---------------------------------------------------------------\n');
n_missing = 0;
for k = FIGURES_TO_MAKE
    tag   = sprintf('S%d', k);
    fname = figFileName(tag);
    fpath = fullfile(OUT_DIR, fname);
    if isfile(fpath)
        info = dir(fpath);
        fprintf('  [OK]   %-42s %8.1f kB\n', fname, info.bytes/1024);
    else
        fprintf(2,'  [MISS] %-42s\n', fname);
        n_missing = n_missing + 1;
    end
end

% [S5] S10 and S11 are written by Combined_figures_main_v7.m (the two figures
% moved out of the main text). Report on them so the SI package is complete.
% FigS9 is built by the separate SLA script; S10 and S11 come from
% Combined_figures_main_v7.m (the two figures moved out of the main text).
for ext = {'S9','S10','S11'}
    en = figFileName(ext{1});
    if isfile(fullfile(OUT_DIR, en))
        info = dir(fullfile(OUT_DIR, en));
        fprintf('  [OK]   %-42s %8.1f kB  (from main figure script)\n', en, info.bytes/1024);
    else
        fprintf('  [--]   %-42s  copy from MAIN_FIGURES\n', en);
    end
end
fprintf('---------------------------------------------------------------\n');

% Ground truth: every PNG actually sitting in the folder right now.
fprintf('\n Full contents of %s :\n', OUT_DIR);
dd = dir(fullfile(OUT_DIR,'*.png'));
if isempty(dd)
    fprintf(2,'   (no PNG files found)\n');
else
    [~,ord] = sort({dd.name});
    for ii = ord
        fprintf('   %-44s %8.1f kB   %s\n', dd(ii).name, dd(ii).bytes/1024, dd(ii).date);
    end
end

if n_missing == 0
    fprintf('\n All eight script-built SI figures are present and correctly named.\n');
else
    fprintf(2,'\n %d figure(s) missing - scroll up for the reported error.\n', n_missing);
end
fprintf('\n=== DONE ===\n');
%  Output file names
%  Single source of truth for the whole script.
function name = figFileName(tag)
    switch upper(tag)
        case 'S1', name = 'FigS1_MLD_Method_Comparison.png';
        case 'S2', name = 'FigS2_HYCOM_Argo_Validation.png';
        case 'S3', name = 'FigS3_Rotary_Spectra_Detail.png';
        case 'S4', name = 'FigS4_Inertial_Spectrograms.png';
        case 'S5', name = 'FigS5_Vertical_Profile_Snapshots.png';
        case 'S6', name = 'FigS6_RL_Depth_Resolved.png';
        case 'S7', name = 'FigS7_RNI_Sensitivity.png';
        case 'S8', name = 'FigS8_PhaseSpace_N2_S2.png';
        % [S4] S9 is now built here. S10 and S11 are the two figures moved out
        % of the main text; they are produced by Combined_figures_main_v7.m,
        % which writes them under these names.
        case 'S9',  name = 'FigS9_SLA_Environment.png';
        case 'S10', name = 'FigS10_Cross_Storm_Summary.png';
        case 'S11', name = 'FigS11_Hysteresis.png';
        otherwise, name = sprintf('Fig%s.png', upper(tag));
    end
end

function report_fig_error(label, ME)
    fprintf(2, '\n*** ERROR building Fig %s: %s ***\n', label, ME.message);
    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks','off'));
end

function report_tile_error(label, ME)
    fprintf(2, '   [tile warning] %s: %s\n', label, ME.message);
end
%  ROBUST FIGURE SAVER
function ok = saveFigureRobust(fig, OUT_DIR, filename)

    ok = false;
    if nargin < 3 || isempty(filename)
        warning('saveFigureRobust called without a filename.'); return;
    end
    if isempty(fig) || ~isgraphics(fig) || ~isvalid(fig)
        warning('Figure handle for %s is invalid (deleted before save).', filename); return;
    end
    if ~exist(OUT_DIR,'dir')
        try, mkdir(OUT_DIR); catch, end
    end

    try
        set(fig,'PaperPositionMode','auto');
        figure(fig);
    catch
    end
    drawnow; pause(0.4); drawnow;

    outfile = fullfile(OUT_DIR, filename);
    if isfile(outfile)
        try
            delete(outfile);
        catch ME
            warning('Could not delete old figure file: %s\n%s', outfile, ME.message);
        end
    end

    saved_ok = local_try_all(fig, outfile, filename);

    if ~saved_ok
        alt = fullfile(pwd, filename);
        warning('Retrying save of %s into the current folder: %s', filename, pwd);
        saved_ok = local_try_all(fig, alt, filename);
        if saved_ok, outfile = alt; end
    end

    if ~saved_ok
        warning('FIGURE SAVE FAILED after all methods: %s', filename); return;
    end

    info = dir(outfile);
    fprintf(' [FIGURE SAVED] %s  (%.1f kB)\n', outfile, info.bytes/1024);
    ok = true;
end

function saved_ok = local_try_all(fig, outfile, filename)
    saved_ok = false;

    % Raster first: complex pcolor/shading-interp meshes and dense text
    % overlays can defeat vector export.
    try
        exportgraphics(fig, outfile, 'Resolution', 300, ...
            'ContentType','image', 'BackgroundColor','white');
        drawnow; saved_ok = file_nonempty(outfile);
    catch ME
        warning('exportgraphics (image) failed for %s:\n%s', filename, ME.message);
    end

    if ~saved_ok
        try
            exportgraphics(fig, outfile, 'Resolution', 300);
            drawnow; saved_ok = file_nonempty(outfile);
        catch ME
            warning('exportgraphics (auto) failed for %s:\n%s', filename, ME.message);
        end
    end

    if ~saved_ok
        try
            print(fig, outfile, '-dpng', '-r300', '-opengl');
            drawnow; saved_ok = file_nonempty(outfile);
        catch ME
            warning('print (-opengl) failed for %s:\n%s', filename, ME.message);
        end
    end

    if ~saved_ok
        try
            print(fig, outfile, '-dpng', '-r200', '-painters');
            drawnow; saved_ok = file_nonempty(outfile);
        catch ME
            warning('print (-painters) failed for %s:\n%s', filename, ME.message);
        end
    end

    if ~saved_ok
        try
            saveas(fig, outfile);
            drawnow; saved_ok = file_nonempty(outfile);
        catch ME
            warning('saveas failed for %s:\n%s', filename, ME.message);
        end
    end
end

function tf = file_nonempty(f)
    tf = false;
    if isfile(f)
        info = dir(f);
        tf = info.bytes > 0;
    end
end
%  SAFE-PLOTTING HELPERS
function ok = safe_xlim(a, b)
    ok = false;
    if isempty(a) || isempty(b), return; end
    a = a(1); b = b(1);
    if isdatetime(a) || isdatetime(b)
        if ~isdatetime(a) || ~isdatetime(b), return; end
        if isnat(a) || isnat(b) || a >= b, return; end
    else
        if ~isfinite(a) || ~isfinite(b) || a >= b, return; end
    end
    try, xlim([a b]); ok = true; catch, end
end

function ok = safe_ylim(lo, hi)
    ok = false;
    if isempty(lo) || isempty(hi), return; end
    lo = double(lo(1)); hi = double(hi(1));
    if ~isfinite(lo) || ~isfinite(hi) || hi <= lo, return; end
    try, ylim([lo hi]); ok = true; catch, end
end

function safe_clim(cl)
    if numel(cl) < 2 || any(~isfinite(cl)), cl = [0 1]; end
    cl = double(cl(:))'; cl = cl(1:2);
    if cl(2) <= cl(1), cl(2) = cl(1) + max(abs(cl(1))*0.1, 1); end
    try
        clim(cl);
    catch
        try, caxis(cl); catch, end %#ok<CAXIS>
    end
end

function w = win_hann(n)
% hann() needs the Signal Processing Toolbox; identical window here.
    n = max(1, round(n));
    if n == 1, w = 1; return; end
    w = 0.5*(1 - cos(2*pi*(0:n-1)'/(n-1)));
end

function dkm = dist_km(lat0, lon0, LAT, LON)
% Haversine great-circle distance in km (no Mapping Toolbox needed).
    Re = 6371.0;
    dlat = deg2rad(LAT - lat0);
    dlon = deg2rad(LON - lon0);
    a = sin(dlat/2).^2 + cosd(lat0).*cosd(LAT).*sin(dlon/2).^2;
    dkm = 2*Re*asin(min(1, sqrt(a)));
end
%  FIG S1 (draft S3) - MLD CRITERION COMPARISON
%  writes FigS1_MLD_Method_Comparison.png
%  Caption: HYCOM native, Delta-rho = 0.03, Delta-sigma_theta = 0.125
%           (paper criterion), Delta-T = 0.5, vs Argo climatology.
function ok = make_figS1_MLDMethods(D, n_cyc, OUT_DIR, beautify, HYCOM_PATHS, Hmax)
    fprintf('\n--- Fig S1 : MLD criterion comparison (slow, re-reads HYCOM) ---\n');
    ok = false;

    c_hycom = [0.40 0.40 0.40];
    c_m1    = [0.85 0.20 0.20];
    c_m2    = [0.10 0.45 0.75];
    c_m3    = [0.55 0.15 0.55];
    c_argo  = [0.10 0.55 0.30];

    THRESH_DRHO  = 0.03;
    THRESH_DT_M2 = 0.2;
    THRESH_DT_M3 = 0.5;
    MLD_MIN = 5; MLD_MAX = 100;

    ARGO_CLIM = struct('KYARR',35,'AMPHAN',25,'FANI',25,'TAUKTAE',30);
    ARGO_CITE = struct( ...
        'KYARR','de Boyer Montegut 2004 AS Oct', ...
        'AMPHAN','Li 2017; Girishkumar 2011 BoB May', ...
        'FANI','Li 2017; Girishkumar 2011 BoB May', ...
        'TAUKTAE','de Boyer Montegut 2004 AS May');

    tmpl = struct('valid',false,'name','','time',datetime.empty(0,1),'keep',[], ...
        'MLD_HYCOM',[],'MLD_M1',[],'MLD_M2',[],'MLD_M3',[], ...
        'MLD_HYCOM_raw',[],'MLD_M1_raw',[],'MLD_M2_raw',[],'MLD_M3_raw',[], ...
        'xmin',NaT,'xmax',NaT,'argo_clim',NaN,'argo_cite','');
    R = repmat(tmpl, n_cyc, 1);
    for i = 1:n_cyc, R(i).name = D(i).name; end

    for i = 1:n_cyc
        storm = D(i).name;
        fprintf('  %s ... ', storm);
        try
            if ~isfield(HYCOM_PATHS, storm)
                warning('No HYCOM path for %s', storm); fprintf('skipped\n'); continue;
            end
            hycom_path = HYCOM_PATHS.(storm);
            if ~exist(hycom_path,'dir')
                warning('HYCOM dir missing: %s', hycom_path); fprintf('skipped\n'); continue;
            end

            T = D(i).T; time = D(i).time;
            lat_tc = D(i).lat; lon_tc = D(i).lon;
            nt = height(T);

            [tsfiles, tstimes, tsfi, tsti] = index_hycom_dir(hycom_path);
            if isempty(tsfiles)
                warning('No .nc in %s', hycom_path); fprintf('skipped\n'); continue;
            end
            fn_ref = tsfiles{1};
            try
                lon_g = ncread(fn_ref,'lon'); lat_g = ncread(fn_ref,'lat');
            catch
                lon_g = ncread(fn_ref,'longitude'); lat_g = ncread(fn_ref,'latitude');
            end
            depth = ncread(fn_ref,'depth'); idz = find(depth <= Hmax);
            z = depth(idz); nz = length(z);

            MLD_HYCOM = nan(nt,1); MLD_M1 = nan(nt,1);
            MLD_M2    = nan(nt,1); MLD_M3 = nan(nt,1);

            for it = 1:nt
                if isnan(lat_tc(it)) || isnan(lon_tc(it)), continue; end
                [dd,ki] = min(abs(tstimes - time(it)));
                if isempty(ki) || dd > hours(2), continue; end
                fn = tsfiles{tsfi(ki)}; tt = tsti(ki);
                ix = find(lon_g >= lon_tc(it)-2 & lon_g <= lon_tc(it)+2);
                iy = find(lat_g >= lat_tc(it)-2 & lat_g <= lat_tc(it)+2);
                if isempty(ix) || isempty(iy), continue; end
                try
                    T_box = ncread(fn,'water_temp',[ix(1) iy(1) 1 tt],[length(ix) length(iy) nz 1]);
                    S_box = ncread(fn,'salinity',  [ix(1) iy(1) 1 tt],[length(ix) length(iy) nz 1]);
                    T_box = squeeze(T_box); S_box = squeeze(S_box);
                catch
                    continue;
                end
                [LonSub,LatSub] = meshgrid(lon_g(ix),lat_g(iy));
                % [S6] Main_combined_v7_FINAL.m forms its T/S profile as
                % mean(tb,[1 2]) over the WHOLE +/-2 deg sub-box, not over the
                % 150 km near-inertial disk. Masking to 150 km here made this
                % figure disagree with the ILD in Table 2 by 3-11 m. Use the
                % same box so the two are comparable.
                mask_NI = ones(size(LatSub));

                T_pf = nan(nz,1); S_pf = nan(nz,1);
                for kk = 1:nz
                    sT = T_box(:,:,kk); sS = S_box(:,:,kk);
                    if ~isequal(size(sT), size(mask_NI)), sT = sT'; sS = sS'; end
                    T_pf(kk) = mean(sT.*mask_NI, 'all','omitnan');
                    S_pf(kk) = mean(sS.*mask_NI, 'all','omitnan');
                end
                valid = isfinite(T_pf) & isfinite(S_pf) & z >= 0;
                if sum(valid) < 3, continue; end
                T_v = T_pf(valid); S_v = S_pf(valid); z_v = z(valid);

                try
                    if exist('gsw_rho','file') == 2
                        SA  = gsw_SA_from_SP(S_v, z_v, lon_tc(it), lat_tc(it));
                        CT  = gsw_CT_from_t(SA, T_v, z_v);
                        rho = gsw_rho(SA, CT, 10);
                    else
                        rho = 1000 + (28 - 0.20*T_v + 0.78*(S_v-35));
                    end
                catch
                    rho = 1000 + (28 - 0.20*T_v + 0.78*(S_v-35));
                end

                [~,iz10] = min(abs(z_v - 10));
                rho_ref = rho(iz10); T_ref = T_v(iz10);

                ix1 = find(rho - rho_ref >= THRESH_DRHO, 1);
                if ~isempty(ix1) && ix1 > 1
                    m = (rho(ix1)-rho(ix1-1))/(z_v(ix1)-z_v(ix1-1)+eps);
                    MLD_M1(it) = z_v(ix1-1) + (rho_ref+THRESH_DRHO - rho(ix1-1))/(m+eps);
                elseif ~isempty(ix1)
                    MLD_M1(it) = z_v(ix1);
                else
                    MLD_M1(it) = max(z_v);
                end

                ix2 = find(T_ref - T_v >= THRESH_DT_M2, 1);
                if ~isempty(ix2) && ix2 > 1
                    m = (T_v(ix2)-T_v(ix2-1))/(z_v(ix2)-z_v(ix2-1)+eps);
                    MLD_M2(it) = z_v(ix2-1) + (T_ref - THRESH_DT_M2 - T_v(ix2-1))/(m+eps);
                elseif ~isempty(ix2)
                    MLD_M2(it) = z_v(ix2);
                else
                    MLD_M2(it) = max(z_v);
                end

                ix3 = find(T_ref - T_v >= THRESH_DT_M3, 1);
                if ~isempty(ix3) && ix3 > 1
                    m = (T_v(ix3)-T_v(ix3-1))/(z_v(ix3)-z_v(ix3-1)+eps);
                    MLD_M3(it) = z_v(ix3-1) + (T_ref - THRESH_DT_M3 - T_v(ix3-1))/(m+eps);
                elseif ~isempty(ix3)
                    MLD_M3(it) = z_v(ix3);
                else
                    MLD_M3(it) = max(z_v);
                end

                MLD_HYCOM(it) = MLD_M2(it);
            end

            MLD_M1    = max(min(MLD_M1,    MLD_MAX), MLD_MIN);
            MLD_M2    = max(min(MLD_M2,    MLD_MAX), MLD_MIN);
            MLD_M3    = max(min(MLD_M3,    MLD_MAX), MLD_MIN);
            MLD_HYCOM = max(min(MLD_HYCOM, MLD_MAX), MLD_MIN);

            valid_any = isfinite(MLD_M1) | isfinite(MLD_M2) | isfinite(MLD_M3);
            if ~any(valid_any)
                warning('No valid MLD retrieved for %s', storm); fprintf('no data\n'); continue;
            end
            i_first = find(valid_any,1,'first'); i_last = find(valid_any,1,'last');
            keep = false(size(valid_any)); keep(i_first:i_last) = true;

            MLD_HYCOM_p = MLD_HYCOM; MLD_M1_p = MLD_M1;
            MLD_M2_p    = MLD_M2;    MLD_M3_p = MLD_M3;
            MLD_HYCOM_p(keep) = fillmissing(MLD_HYCOM(keep),'linear','EndValues','none');
            MLD_M1_p(keep)    = fillmissing(MLD_M1(keep),   'linear','EndValues','none');
            MLD_M2_p(keep)    = fillmissing(MLD_M2(keep),   'linear','EndValues','none');
            MLD_M3_p(keep)    = fillmissing(MLD_M3(keep),   'linear','EndValues','none');

            R(i).MLD_HYCOM     = MLD_HYCOM_p;
            R(i).MLD_HYCOM_raw = MLD_HYCOM;
            R(i).MLD_M1        = MLD_M1_p;
            R(i).MLD_M1_raw    = MLD_M1;

            % The M2 series shows the density criterion actually used in
            % the paper (Delta-sigma_theta = 0.125) when available.
            if ~isempty(D(i).MLD_main) && numel(D(i).MLD_main) == numel(MLD_M2)
                MLD_main_p = D(i).MLD_main;
                MLD_main_p(keep) = fillmissing(D(i).MLD_main(keep),'linear','EndValues','none');
                R(i).MLD_M2     = MLD_main_p;
                R(i).MLD_M2_raw = D(i).MLD_main;
            else
                R(i).MLD_M2     = MLD_M2_p;
                R(i).MLD_M2_raw = MLD_M2;
            end

            R(i).MLD_M3     = MLD_M3_p;
            R(i).MLD_M3_raw = MLD_M3;

            % [S8] Restrict every criterion to the timesteps the MAIN analysis
            % actually kept. Main_combined_v7_FINAL.m drops whole timesteps at
            % the shallow-water guard and on missing reads, so its ILD mean is
            % over 39/60/36 steps for Amphan/Fani/Tauktae, not the full track.
            % Averaging this figure over its own larger set made the ILD look
            % 2-5 m deeper than Table 2 for those three (Kyarr, with 77/77
            % valid, matched exactly). Use one common sample for all columns.
            if ~isempty(D(i).MLD_main) && numel(D(i).MLD_main) == numel(MLD_M1)
                bad_main = ~isfinite(D(i).MLD_main);
                R(i).MLD_HYCOM_raw(bad_main) = NaN;
                R(i).MLD_M1_raw(bad_main)    = NaN;
                R(i).MLD_M3_raw(bad_main)    = NaN;
            end
            R(i).time       = time;
            R(i).keep       = keep;
            R(i).xmin       = min(time(keep));
            R(i).xmax       = max(time(keep));
            R(i).argo_clim  = ARGO_CLIM.(storm);
            R(i).argo_cite  = ARGO_CITE.(storm);
            R(i).valid      = true;
            fprintf('done\n');
        catch ME
            fprintf(2,'failed (%s)\n', ME.message);
        end
    end

    all_mld = [];
    for i = 1:n_cyc
        if ~R(i).valid, continue; end
        k = R(i).keep;
        all_mld = [all_mld; R(i).MLD_HYCOM(k); R(i).MLD_M1(k); ...
                            R(i).MLD_M2(k);    R(i).MLD_M3(k)]; %#ok<AGROW>
    end
    all_mld = all_mld(isfinite(all_mld));
    if isempty(all_mld), ylim_a = [0 80];
    else, ylim_a = [0, max(prctile(all_mld,99)*1.15, 80)]; end
    ylim_b = ylim_a;

    fig = figure('Name','Fig S1 - MLD Criterion Comparison', ...
        'Units','normalized','Position',[0.05 0.10 0.92 0.78]);
    tiledlayout(2, n_cyc, 'TileSpacing','compact','Padding','compact');

    % ---- Row (a) : time series ----
    for i = 1:n_cyc
        try
            nexttile(i);
            if ~R(i).valid, title(R(i).name,'FontSize',12); axis off; continue; end
            s = R(i); k = s.keep; tk = s.time(k);

            % [S10] Values sitting exactly on the MLD_MAX clamp are not mixed
            % layers, they are the limit of the search range, and they appeared
            % as vertical spikes to 100 m at the end of three records. Blank
            % them for display. The TC-period means in panel (b) are unaffected,
            % since those timesteps lie outside the active window.
            uc = @(x) local_unclamp(x, MLD_MIN, MLD_MAX);
            h_hy = plot(tk, uc(s.MLD_HYCOM(k)),'Color',c_hycom,'LineWidth',1.5,'LineStyle',':'); hold on;
            h_m1 = plot(tk, uc(s.MLD_M1(k)),   'Color',c_m1,   'LineWidth',1.8,'LineStyle','--');
            h_m2 = plot(tk, uc(s.MLD_M2(k)),   'Color',c_m2,   'LineWidth',2.5,'LineStyle','-');
            h_m3 = plot(tk, uc(s.MLD_M3(k)),   'Color',c_m3,   'LineWidth',1.8,'LineStyle','-.');
            yline(s.argo_clim,'-',sprintf('Argo %d m', s.argo_clim), ...
                'Color',c_argo,'LineWidth',2.0,'FontSize',8, ...
                'LabelHorizontalAlignment','left','HandleVisibility','off');

            safe_xlim(s.xmin, s.xmax); safe_ylim(ylim_a(1), ylim_a(2));
            set(gca,'YDir','reverse');
            title(s.name,'FontSize',12);
            set_xax_date(gca);
            xlabel('Date (UTC)','FontSize',10);
            if i == 1
                ylabel({'\bf(a)';'MLD (m)'},'FontSize',11,'FontWeight','bold');
                legend([h_hy h_m1 h_m2 h_m3], ...
                    {'\DeltaT = 0.2 ^\circC (paper ILD)','\Delta\rho = 0.03 kg m^{-3}', ...
                     '\Delta\sigma_\theta = 0.125 (paper MLD)','\DeltaT = 0.5 ^\circC'}, ...
                    'Box','off','Location','southwest','FontSize',8);
            else
                ylabel('MLD (m)','FontSize',11,'FontWeight','bold');
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S1 row a tile %d', i), ME);
        end
    end

    % ---- Row (b) : means +/- sd against Argo climatology ----
    for i = 1:n_cyc
        try
            nexttile(n_cyc + i);
            if ~R(i).valid, axis off; continue; end
            s = R(i); k = s.keep;

            means = [mean(s.MLD_HYCOM_raw(k),'omitnan'), mean(s.MLD_M1_raw(k),'omitnan'), ...
                     mean(s.MLD_M2_raw(k),'omitnan'),    mean(s.MLD_M3_raw(k),'omitnan')];
            stds  = [std(s.MLD_HYCOM_raw(k),'omitnan'),  std(s.MLD_M1_raw(k),'omitnan'), ...
                     std(s.MLD_M2_raw(k),'omitnan'),     std(s.MLD_M3_raw(k),'omitnan')];
            means(~isfinite(means)) = NaN;
            stds(~isfinite(stds))   = 0;

            b = bar(1:4, means, 0.6, 'FaceColor','flat','EdgeColor','k','LineWidth',1.2);
            b.CData = [c_hycom; c_m1; c_m2; c_m3]; hold on;
            errorbar(1:4, means, stds, 'k','LineStyle','none','LineWidth',1.4,'CapSize',8);
            % [S15] label moved to the right-hand end of the line; on the left it
            % overlapped the first bar's value annotation.
            yline(s.argo_clim, '-', sprintf('Argo %d m', s.argo_clim), ...
                'Color',c_argo,'LineWidth',2.0,'FontSize',9, ...
                'LabelHorizontalAlignment','right', ...
                'LabelVerticalAlignment','bottom','HandleVisibility','off');
            for j = 1:4
                if isfinite(means(j))
                    text(j, means(j) + max(stds(j),2) + 1, sprintf('%.0f', means(j)), ...
                        'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
                end
            end

            bias = means - s.argo_clim;
            lbls = {'\DeltaT_{0.2}','\Delta\rho','\Delta\sigma_\theta','\DeltaT_{0.5}'};
            text(0.6, ylim_b(2)*0.85, 'Argo bias (m):','FontSize',8,'FontWeight','bold');
            for j = 1:4
                text(0.6, ylim_b(2)*(0.78 - 0.06*j), ...
                    sprintf('  %s: %+.1f', lbls{j}, bias(j)),'FontSize',8);
            end

            set(gca,'XTick',1:4,'XTickLabel',lbls,'TickLabelInterpreter','tex');
            safe_ylim(ylim_b(1), ylim_b(2));
            if i == 1
                ylabel({'\bf(b)';'Mean MLD (m)'},'FontSize',11,'FontWeight','bold');
                h_argo = plot(NaN,NaN,'-','Color',c_argo,'LineWidth',2.0);
                legend(h_argo,{'Argo climatology'},'Box','off','Location','northeast','FontSize',8);
            else
                ylabel('Mean MLD (m)','FontSize',11,'FontWeight','bold');
            end
            xlabel(sprintf('%s (%s)', s.name, s.argo_cite),'FontSize',9);
            beautify();
        catch ME
            report_tile_error(sprintf('S1 row b tile %d', i), ME);
        end
    end

    % Titles omitted - the AGU caption sits below the figure
    ok = saveFigureRobust(fig, OUT_DIR, figFileName('S1'));
    if ok, fprintf(' Saved Fig S1\n'); else, fprintf(2,' Fig S1 was NOT saved\n'); end

    try
        fprintf('\n[FIG S1 - MLD CRITERION COMPARISON]\n');
        fprintf('%-8s | %-9s | %-9s | %-9s | %-9s | %-9s\n', ...
            'Storm','dT0.2(ILD)','d-rho0.03','d-sig0.125','dT0.5','Argo');
        fprintf('%s\n', repmat('-',1,70));
        for i = 1:n_cyc
            if ~R(i).valid, continue; end
            s = R(i); k = s.keep;
            fprintf('%-8s | %-9.1f | %-9.1f | %-9.1f | %-9.1f | %-9d\n', ...
                s.name, mean(s.MLD_HYCOM_raw(k),'omitnan'), mean(s.MLD_M1_raw(k),'omitnan'), ...
                mean(s.MLD_M2_raw(k),'omitnan'),  mean(s.MLD_M3_raw(k),'omitnan'), s.argo_clim);
        end
    catch ME
        report_tile_error('S1 QC block', ME);
    end
end
%  FIG S2 (draft S2) - HYCOM vs ARGO VALIDATION SCATTER
%  writes FigS2_HYCOM_Argo_Validation.png
%  Reads the merged validation report; no re-computation.
function ok = make_figS2_HycomArgo(STORM_NAMES, MERGED_REPORT, OUT_DIR, beautify)
    fprintf('\n--- Fig S2 : HYCOM-Argo validation scatter ---\n');
    ok = false;
    n_cyc = numel(STORM_NAMES);

    storm_color = containers.Map( ...
        {'KYARR','AMPHAN','FANI','TAUKTAE'}, ...
        {[0.85 0.40 0.10], [0.10 0.40 0.75], [0.20 0.65 0.45], [0.80 0.30 0.55]});

    if ~isfile(MERGED_REPORT)
        warning('Merged validation report not found: %s -- Fig S2 skipped.', MERGED_REPORT);
        return;
    end

    txt_all = fileread(MERGED_REPORT);
    blocks  = strsplit(txt_all, 'HYCOM vs ARGO FULL VALIDATION REPORT');

    parsed = containers.Map('KeyType','char','ValueType','any');
    for b = 2:numel(blocks)
        blk = blocks{b};
        nm  = regexp(blk, 'Cyclone\s+(\w+)', 'tokens', 'once');
        if isempty(nm), continue; end
        storm = upper(nm{1});

        i_start = strfind(blk, 'PER-PROFILE COLOCATION TABLE');
        if isempty(i_start), continue; end
        sub = blk(i_start:end);
        i_dash = strfind(sub, '----');
        if isempty(i_dash), continue; end
        data_block = sub(i_dash(1):end);
        i_end = strfind(data_block, 'OUTPUT FILES');
        if ~isempty(i_end), data_block = data_block(1:i_end(1)-1); end

        lines = strsplit(data_block, newline);
        T5a = nan(500,1); T5h = T5a; S5a = T5a; S5h = T5a; MLDa = T5a; MLDh = T5a;
        cnt = 0;
        for L = 1:numel(lines)
            ln = strtrim(lines{L});
            if isempty(ln) || ln(1)=='-' || ln(1)=='=', continue; end
            tok = strsplit(ln);
            if numel(tok) < 12, continue; end
            if ~isfinite(str2double(tok{1})), continue; end
            cnt = cnt + 1;
            T5a(cnt)  = str2double(tok{7});  T5h(cnt)  = str2double(tok{8});
            S5a(cnt)  = str2double(tok{9});  S5h(cnt)  = str2double(tok{10});
            MLDa(cnt) = str2double(tok{11}); MLDh(cnt) = str2double(tok{12});
        end
        if cnt == 0, continue; end
        parsed(storm) = struct('T5a',T5a(1:cnt),'T5h',T5h(1:cnt), ...
                               'S5a',S5a(1:cnt),'S5h',S5h(1:cnt), ...
                               'MLDa',MLDa(1:cnt),'MLDh',MLDh(1:cnt));
        fprintf('  Parsed %s (n=%d)\n', storm, cnt);
    end

    tmpl = struct('valid',false,'name','','T5a',[],'T5h',[],'S5a',[],'S5h',[], ...
        'MLDa',[],'MLDh',[],'color',[0 0 0]);
    V = repmat(tmpl, n_cyc, 1);
    for i = 1:n_cyc
        storm = STORM_NAMES{i};
        V(i).name = storm;
        if ~isKey(parsed, upper(storm))
            warning('No block for %s in merged report', storm); continue;
        end
        s = parsed(upper(storm));
        V(i).T5a = s.T5a;  V(i).T5h = s.T5h;
        V(i).S5a = s.S5a;  V(i).S5h = s.S5h;
        V(i).MLDa = s.MLDa; V(i).MLDh = s.MLDh;
        if isKey(storm_color, storm), V(i).color = storm_color(storm);
        else, V(i).color = [0.3 0.3 0.3]; end
        V(i).valid = true;
    end

    Tall = []; Sall = []; Mall = [];
    for i = 1:n_cyc
        if ~V(i).valid, continue; end
        Tall = [Tall; V(i).T5a;  V(i).T5h];  %#ok<AGROW>
        Sall = [Sall; V(i).S5a;  V(i).S5h];  %#ok<AGROW>
        Mall = [Mall; V(i).MLDa; V(i).MLDh]; %#ok<AGROW>
    end
    Tall = Tall(isfinite(Tall)); Sall = Sall(isfinite(Sall)); Mall = Mall(isfinite(Mall));
    if isempty(Tall), Tlim = [25 32]; else, Tlim = [floor(min(Tall))-0.5, ceil(max(Tall))+0.5]; end
    if isempty(Sall), Slim = [30 37]; else, Slim = [floor(min(Sall))-0.5, ceil(max(Sall))+0.5]; end
    if isempty(Mall), Mlim = [0 100]; else, Mlim = [0, ceil(max(Mall)/10)*10+5]; end

    fig = figure('Name','Fig S2 - HYCOM-Argo Validation', ...
        'Units','normalized','Position',[0.04 0.05 0.92 0.88]);
    tiledlayout(3, n_cyc, 'TileSpacing','compact','Padding','compact');

    % ---- Row (a) : temperature at 5 m ----
    for i = 1:n_cyc
        try
            nexttile(i);
            if ~V(i).valid, title(STORM_NAMES{i},'FontSize',12); axis off; continue; end
            v = V(i); okm = isfinite(v.T5a) & isfinite(v.T5h);
            plot(Tlim,Tlim,'k--','LineWidth',1.0); hold on;
            scatter(v.T5a(okm), v.T5h(okm), 36, v.color,'filled','MarkerEdgeColor','k','LineWidth',0.5);
            [bias,rmse,r] = scatter_stats(v.T5a(okm), v.T5h(okm));
            text(0.04,0.96,sprintf('Bias = %+.2f\nRMSE = %.2f\nR = %.2f\nN = %d',bias,rmse,r,sum(okm)), ...
                'Units','normalized','VerticalAlignment','top','FontSize',8, ...
                'BackgroundColor','w','EdgeColor',[0.5 0.5 0.5],'Margin',3);
            safe_xlim(Tlim(1),Tlim(2)); safe_ylim(Tlim(1),Tlim(2)); axis square;
            title(v.name,'FontSize',12);
            xlabel('Argo T at 5 m (^\circC)','FontSize',9);
            if i == 1
                ylabel({'\bf(a)';'HYCOM T at 5 m (^\circC)'},'FontSize',10,'FontWeight','bold');
            else
                ylabel('HYCOM T (^\circC)','FontSize',10,'FontWeight','bold');
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S2 row a tile %d', i), ME);
        end
    end

    % ---- Row (b) : salinity at 5 m ----
    for i = 1:n_cyc
        try
            nexttile(n_cyc + i);
            if ~V(i).valid, axis off; continue; end
            v = V(i); okm = isfinite(v.S5a) & isfinite(v.S5h);
            plot(Slim,Slim,'k--','LineWidth',1.0); hold on;
            scatter(v.S5a(okm), v.S5h(okm), 36, v.color,'filled','MarkerEdgeColor','k','LineWidth',0.5);
            [bias,rmse,r] = scatter_stats(v.S5a(okm), v.S5h(okm));
            text(0.04,0.96,sprintf('Bias = %+.2f\nRMSE = %.2f\nR = %.2f\nN = %d',bias,rmse,r,sum(okm)), ...
                'Units','normalized','VerticalAlignment','top','FontSize',8, ...
                'BackgroundColor','w','EdgeColor',[0.5 0.5 0.5],'Margin',3);
            safe_xlim(Slim(1),Slim(2)); safe_ylim(Slim(1),Slim(2)); axis square;
            xlabel('Argo S at 5 m (psu)','FontSize',9);
            if i == 1
                ylabel({'\bf(b)';'HYCOM S at 5 m (psu)'},'FontSize',10,'FontWeight','bold');
            else
                ylabel('HYCOM S (psu)','FontSize',10,'FontWeight','bold');
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S2 row b tile %d', i), ME);
        end
    end

    % ---- Row (c) : MLD (0.2 degC threshold ref 5 m, both datasets) ----
    for i = 1:n_cyc
        try
            nexttile(2*n_cyc + i);
            if ~V(i).valid, axis off; continue; end
            v = V(i); okm = isfinite(v.MLDa) & isfinite(v.MLDh);
            plot(Mlim,Mlim,'k--','LineWidth',1.0); hold on;
            scatter(v.MLDa(okm), v.MLDh(okm), 36, v.color,'filled','MarkerEdgeColor','k','LineWidth',0.5);
            [bias,rmse,r] = scatter_stats(v.MLDa(okm), v.MLDh(okm));
            text(0.04,0.96,sprintf('Bias = %+.1f m\nRMSE = %.1f\nR = %.2f\nN = %d',bias,rmse,r,sum(okm)), ...
                'Units','normalized','VerticalAlignment','top','FontSize',8, ...
                'BackgroundColor','w','EdgeColor',[0.5 0.5 0.5],'Margin',3);
            safe_xlim(Mlim(1),Mlim(2)); safe_ylim(Mlim(1),Mlim(2)); axis square;
            xlabel('Argo MLD (m)','FontSize',9);
            if i == 1
                ylabel({'\bf(c)';'HYCOM MLD (m)'},'FontSize',10,'FontWeight','bold');
            else
                ylabel('HYCOM MLD (m)','FontSize',10,'FontWeight','bold');
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S2 row c tile %d', i), ME);
        end
    end

    ok = saveFigureRobust(fig, OUT_DIR, figFileName('S2'));
    if ok, fprintf(' Saved Fig S2\n'); else, fprintf(2,' Fig S2 was NOT saved\n'); end

    try
        fprintf('\n[FIG S2 - HYCOM-ARGO VALIDATION]\n');
        fprintf('%-8s | %-5s | %-9s | %-9s | %-9s\n','Storm','N','T_bias','S_bias','MLD_bias');
        fprintf('%s\n', repmat('-',1,55));
        for i = 1:n_cyc
            if ~V(i).valid, continue; end
            v = V(i);
            okT = isfinite(v.T5a)  & isfinite(v.T5h);
            okS = isfinite(v.S5a)  & isfinite(v.S5h);
            okM = isfinite(v.MLDa) & isfinite(v.MLDh);
            fprintf('%-8s | %-5d | %+9.3f | %+9.3f | %+9.2f\n', v.name, sum(okT), ...
                mean(v.T5h(okT)-v.T5a(okT)), mean(v.S5h(okS)-v.S5a(okS)), ...
                mean(v.MLDh(okM)-v.MLDa(okM)));
        end
    catch ME
        report_tile_error('S2 QC block', ME);
    end
end

function [bias, rmse, r] = scatter_stats(x, y)
    bias = NaN; rmse = NaN; r = NaN;
    if isempty(x) || isempty(y), return; end
    bias = mean(y - x);
    rmse = sqrt(mean((y - x).^2));
    if numel(x) > 1 && std(x) > 0 && std(y) > 0
        cc = corrcoef(x,y); r = cc(1,2);
    end
end
%  FIG S3 (draft S5) - DETAILED ROTARY SPECTRA
%  writes FigS3_Rotary_Spectra_Detail.png
function ok = make_figS3_RotarySpectra(D, n_cyc, OUT_DIR, beautify, Fs_cpd, f_M2, f_K1)
    fprintf('\n--- Fig S3 : Detailed rotary spectra ---\n');
    ok = false;

    c_cw       = [0.05 0.32 0.58];
    c_cw_fill  = [0.55 0.75 0.95];
    c_ccw      = [0.85 0.45 0.10];
    c_ccw_fill = [0.96 0.75 0.55];
    c_f        = [0.10 0.10 0.10];
    c_2f       = [0.55 0.55 0.55];
    c_M2c      = [0.50 0.30 0.65];
    c_K1c      = [0.30 0.55 0.25];
    c_strong_cw  = [0.85 0.92 1.00];
    c_strong_ccw = [1.00 0.88 0.85];
    c_weak       = [0.95 0.95 0.95];
    c_RC_line    = [0.30 0.10 0.50];

    tmpl = struct('valid',false,'name','','f_cw',[],'P_cw',[],'f_ccw',[],'P_ccw',[], ...
        'P_cw_lo',[],'P_cw_hi',[],'P_ccw_lo',[],'P_ccw_hi',[],'f_inert',NaN, ...
        'f_common',[],'RC_for_plot',[],'RC_at_f',NaN,'peak_cw_f',NaN,'pk_cw_pwr',NaN);
    SPEC = repmat(tmpl, n_cyc, 1);
    for i = 1:n_cyc, SPEC(i).name = D(i).name; end

    has_sig = exist('butter','file') == 2 && exist('filtfilt','file') == 2;

    for i = 1:n_cyc
        try
            d = D(i);
            if isempty(d.u_surf) || isempty(d.v_surf), continue; end

            u = fillmissing(d.u_surf,'linear');
            v = fillmissing(d.v_surf,'linear');
            if length(u) < 12 || all(~isfinite(u)) || all(~isfinite(v)), continue; end
            u = detrend(u,'linear'); v = detrend(v,'linear');

            f_inert = d.f_inertial;
            fnyq = Fs_cpd/2;
            f_lo = 0.5*f_inert/fnyq; f_hi = min(1.5*f_inert/fnyq, 0.95);
            if has_sig && length(u) >= 12 && f_hi > f_lo && f_lo > 0
                [b_bp,a_bp] = butter(2,[f_lo,f_hi],'bandpass');
                u = filtfilt(b_bp,a_bp,u); v = filtfilt(b_bp,a_bp,v);
            else
                u = u - movmean(u,8,'omitnan'); v = v - movmean(v,8,'omitnan');
            end

            N = length(u); w = win_hann(N);
            W = (u(:)+1i*v(:)).*w(:);
            Y = fftshift(fft(W));
            df = Fs_cpd/N; f_vec = (-N/2:N/2-1)*df;
            P = abs(Y).^2 / (sum(w.^2)*df); P(P<=0) = NaN;

            idx_cw = f_vec < 0; idx_ccw = f_vec > 0;
            f_cw  = abs(f_vec(idx_cw));  P_cw  = P(idx_cw);
            f_ccw = f_vec(idx_ccw);      P_ccw = P(idx_ccw);
            [f_cw,jc]   = sort(f_cw);  P_cw  = P_cw(jc);
            [f_ccw,jcc] = sort(f_ccw); P_ccw = P_ccw(jcc);

            SPEC(i).valid    = true;
            SPEC(i).f_cw     = f_cw(:);  SPEC(i).P_cw  = P_cw(:);
            SPEC(i).f_ccw    = f_ccw(:); SPEC(i).P_ccw = P_ccw(:);
            SPEC(i).P_cw_lo  = SPEC(i).P_cw  * 0.6;
            SPEC(i).P_cw_hi  = SPEC(i).P_cw  * 1.66;
            SPEC(i).P_ccw_lo = SPEC(i).P_ccw * 0.6;
            SPEC(i).P_ccw_hi = SPEC(i).P_ccw * 1.66;
            SPEC(i).f_inert  = f_inert;

            f_common = linspace(0.05, 4, 250);
            P_cw_i   = interp1(f_cw, P_cw,  f_common,'linear',0);
            P_ccw_i  = interp1(f_ccw,P_ccw, f_common,'linear',0);
            SPEC(i).f_common    = f_common;
            SPEC(i).RC_for_plot = (P_cw_i - P_ccw_i)./(P_cw_i + P_ccw_i + eps);
            [~,idx_f] = min(abs(f_common - f_inert));
            SPEC(i).RC_at_f = SPEC(i).RC_for_plot(idx_f);

            in_band = f_cw >= 0.5*f_inert & f_cw <= 1.5*f_inert;
            if any(in_band)
                [pkv,pki] = max(P_cw(in_band));
                f_in = f_cw(in_band);
                SPEC(i).peak_cw_f = f_in(pki); SPEC(i).pk_cw_pwr = pkv;
            end
        catch ME
            report_tile_error(sprintf('S3 prep %s', D(i).name), ME);
        end
    end

    all_p = [];
    for i = 1:n_cyc
        if SPEC(i).valid, all_p = [all_p; SPEC(i).P_cw(:); SPEC(i).P_ccw(:)]; end %#ok<AGROW>
    end
    all_p(~isfinite(all_p) | all_p <= 0) = [];
    if ~isempty(all_p)
        ylim_a_log = [10^floor(log10(min(all_p)*0.5)), 10^ceil(log10(max(all_p)*1.3))];
    else
        ylim_a_log = [1e-4, 1e2];
    end

    fig = figure('Name','Fig S3 - Detailed Rotary Spectra', ...
        'Units','normalized','Position',[0.05 0.05 0.92 0.88]);
    tiledlayout(2, n_cyc, 'TileSpacing','compact','Padding','compact');

    % ---- Row (a) : CW / CCW spectra ----
    for i = 1:n_cyc
        try
            nexttile(i);
            if ~SPEC(i).valid, title(SPEC(i).name,'FontSize',12); axis off; continue; end
            s = SPEC(i);

            f_band_cw = [s.f_cw(:); flipud(s.f_cw(:))];
            P_band_cw = [s.P_cw_lo(:); flipud(s.P_cw_hi(:))];
            nm = min(length(f_band_cw),length(P_band_cw));
            f_band_cw = f_band_cw(1:nm); P_band_cw = P_band_cw(1:nm);
            vcw = isfinite(f_band_cw) & isfinite(P_band_cw) & P_band_cw > 0;
            if any(vcw)
                fill(f_band_cw(vcw), P_band_cw(vcw), c_cw_fill, ...
                    'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
            end
            hold on;

            f_band_ccw = [s.f_ccw(:); flipud(s.f_ccw(:))];
            P_band_ccw = [s.P_ccw_lo(:); flipud(s.P_ccw_hi(:))];
            nm = min(length(f_band_ccw),length(P_band_ccw));
            f_band_ccw = f_band_ccw(1:nm); P_band_ccw = P_band_ccw(1:nm);
            vcc = isfinite(f_band_ccw) & isfinite(P_band_ccw) & P_band_ccw > 0;
            if any(vcc)
                fill(f_band_ccw(vcc), P_band_ccw(vcc), c_ccw_fill, ...
                    'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
            end

            h_cw  = loglog(s.f_cw,  s.P_cw,  'Color',c_cw, 'LineWidth',2.2);
            h_ccw = loglog(s.f_ccw, s.P_ccw, '--','Color',c_ccw,'LineWidth',1.8);

            xline(s.f_inert,  '-','Color',c_f,  'LineWidth',1.4,'HandleVisibility','off');
            xline(2*s.f_inert,':','Color',c_2f, 'LineWidth',1.0,'HandleVisibility','off');
            xline(f_M2,      '--','Color',c_M2c,'LineWidth',1.0,'HandleVisibility','off');
            xline(f_K1,       ':','Color',c_K1c,'LineWidth',1.0,'HandleVisibility','off');

            ytxt = ylim_a_log(2)*0.4;
            text(s.f_inert*1.04,   ytxt,      'f','FontSize',9,'FontWeight','bold','Color',c_f);
            text(2*s.f_inert*1.04, ytxt,     '2f','FontSize',8,'Color',c_2f);
            text(f_M2*1.04,        ytxt,    'M_2','FontSize',8,'Color',c_M2c);
            text(f_K1*1.04,        ytxt*0.5,'K_1','FontSize',8,'Color',c_K1c);

            if isfinite(s.peak_cw_f)
                plot(s.peak_cw_f, s.pk_cw_pwr,'o', ...
                    'MarkerEdgeColor','k','MarkerFaceColor',c_cw,'MarkerSize',8);
                text(s.peak_cw_f*1.05, s.pk_cw_pwr*1.4, ...
                    sprintf('peak = %.3f cpd', s.peak_cw_f), ...
                    'FontSize',8,'FontWeight','bold','BackgroundColor','w');
                text(0.09, ylim_a_log(1)*3, sprintf('f_{pk}/f = %.2f', s.peak_cw_f/s.f_inert), ...
                    'FontSize',9,'FontWeight','bold','BackgroundColor','w', ...
                    'EdgeColor',[0.5 0.5 0.5],'Margin',2);
            end

            safe_xlim(0.08, 4); safe_ylim(ylim_a_log(1), ylim_a_log(2));
            set(gca,'XScale','log','YScale','log');
            title(s.name,'FontSize',12);
            if i == 1
                ylabel({'\bf(a)';'Spectral density';'(m^2 s^{-2} cpd^{-1})'}, ...
                    'FontSize',11,'FontWeight','bold');
                legend([h_cw h_ccw],{'CW (NIW signature)','CCW'}, ...
                    'Box','off','Location','southwest','FontSize',9);
            else
                ylabel('Spectral density (m^2 s^{-2} cpd^{-1})','FontSize',11,'FontWeight','bold');
            end
            xlabel('Frequency (cpd)','FontSize',10);
            beautify();
        catch ME
            report_tile_error(sprintf('S3 row a tile %d', i), ME);
        end
    end

    % ---- Row (b) : rotary coefficient ----
    for i = 1:n_cyc
        try
            nexttile(n_cyc + i);
            if ~SPEC(i).valid, axis off; continue; end
            s = SPEC(i);

            x_band = [0 4 4 0];
            patch(x_band,[0.5 0.5 1.0 1.0],    c_strong_cw, 'EdgeColor','none','FaceAlpha',0.6,'HandleVisibility','off');
            hold on;
            patch(x_band,[-1.0 -1.0 -0.5 -0.5],c_strong_ccw,'EdgeColor','none','FaceAlpha',0.6,'HandleVisibility','off');
            patch(x_band,[-0.5 -0.5 0.5 0.5],  c_weak,      'EdgeColor','none','FaceAlpha',0.6,'HandleVisibility','off');
            yline( 0.5,'k--','LineWidth',1.0,'HandleVisibility','off');
            yline( 0.0,'k-', 'LineWidth',0.8,'HandleVisibility','off');
            yline(-0.5,'k--','LineWidth',1.0,'HandleVisibility','off');

            nm = min(length(s.P_cw),length(s.P_ccw));
            P_total = s.P_cw(1:nm) + s.P_ccw(1:nm);
            mask_RC = (interp1(s.f_cw, s.P_cw, s.f_common,'linear',0) + ...
                       interp1(s.f_ccw,s.P_ccw,s.f_common,'linear',0)) > prctile(P_total,25);
            RC_plot = s.RC_for_plot; RC_plot(~mask_RC) = NaN;
            plot(s.f_common, RC_plot, 'Color',c_RC_line,'LineWidth',2.2);

            xline(s.f_inert,'k-','LineWidth',1.4,'HandleVisibility','off');
            xline(f_M2,'--','Color',c_M2c,'LineWidth',1.0,'HandleVisibility','off');

            if isfinite(s.RC_at_f)
                plot(s.f_inert, s.RC_at_f,'o', ...
                    'MarkerEdgeColor','k','MarkerFaceColor',c_RC_line,'MarkerSize',9);
                text(s.f_inert*1.05, s.RC_at_f, sprintf('|RC| = %.3f', abs(s.RC_at_f)), ...
                    'FontSize',9,'FontWeight','bold','BackgroundColor','w', ...
                    'EdgeColor',[0.5 0.5 0.5],'Margin',2);
            end
            text(s.f_inert*1.04, 0.88,'f','FontSize',9,'FontWeight','bold');

            safe_xlim(0, 4); safe_ylim(-1.05, 1.05);
            xlabel('Frequency (cpd)','FontSize',10);
            if i == 1
                ylabel({'\bf(b)';'Rotary coefficient C_R';'(positive = CW dom.)'}, ...
                    'FontSize',11,'FontWeight','bold');
                h_scw  = patch(NaN,NaN,c_strong_cw, 'EdgeColor','none','FaceAlpha',0.6);
                h_w    = patch(NaN,NaN,c_weak,      'EdgeColor','none','FaceAlpha',0.6);
                h_sccw = patch(NaN,NaN,c_strong_ccw,'EdgeColor','none','FaceAlpha',0.6);
                legend([h_scw h_w h_sccw], ...
                    {'Strong CW (>+0.5)','Weak (|C_R|<0.5)','Strong CCW (<-0.5)'}, ...
                    'Box','off','Location','southeast','FontSize',8);
            else
                ylabel('Rotary coefficient C_R','FontSize',11,'FontWeight','bold');
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S3 row b tile %d', i), ME);
        end
    end

    ok = saveFigureRobust(fig, OUT_DIR, figFileName('S3'));
    if ok, fprintf(' Saved Fig S3\n'); else, fprintf(2,' Fig S3 was NOT saved\n'); end

    try
        fprintf('\n[FIG S3 - DETAILED ROTARY SPECTRA]\n');
        fprintf('%-8s | %-10s | %-10s | %-10s\n','Storm','f_pk(cpd)','f_pk/f','|RC|@f');
        fprintf('%s\n', repmat('-',1,50));
        for i = 1:n_cyc
            if ~SPEC(i).valid, continue; end
            s = SPEC(i);
            fprintf('%-8s | %-10.3f | %-10.3f | %-10.3f\n', ...
                s.name, s.peak_cw_f, s.peak_cw_f/s.f_inert, abs(s.RC_at_f));
        end
    catch ME
        report_tile_error('S3 QC block', ME);
    end
end
%  FIG S4 (draft S6) - INERTIAL SPECTROGRAMS
%  writes FigS4_Inertial_Spectrograms.png
%  Caption: Hann windows, 85% overlap, NFFT = 256 (df ~ 0.03 cpd)
function ok = make_figS4_Spectrograms(D, n_cyc, OUT_DIR, beautify, Fs_cpd, f_M2)
    fprintf('\n--- Fig S4 : Inertial spectrograms ---\n');
    ok = false;

    c_finert = [0.10 0.10 0.10];
    c_2f     = [0.55 0.55 0.55];
    c_M2c    = [0.50 0.30 0.65];
    c_peak   = [0.85 0.20 0.20];
    use_cmo  = exist('cmocean.m','file') == 2;

    NFFT_PAD     = 256;    % df = Fs/NFFT = 0.03125 cpd for Fs = 8 cpd
    OVERLAP_FRAC = 0.85;   % matches the SI caption
    WIN_MIN      = 16;
    WIN_MAX      = 48;

    tmpl = struct('valid',false,'name','','f_plot',[],'time_sp',datetime.empty(0,1), ...
        'P_plot',[],'f_inert',NaN,'time_peak',NaT,'xmin',NaT,'xmax',NaT);
    SP = repmat(tmpl, n_cyc, 1);
    for i = 1:n_cyc, SP(i).name = D(i).name; end

    for i = 1:n_cyc
        try
            d = D(i);
            if isempty(d.u_surf) || isempty(d.v_surf), continue; end
            u = fillmissing(d.u_surf,'linear','EndValues','nearest');
            v = fillmissing(d.v_surf,'linear','EndValues','nearest');
            time = d.time;

            valid_uv = isfinite(u) & isfinite(v);
            if ~any(valid_uv), continue; end
            i_first = find(valid_uv,1,'first'); i_last = find(valid_uv,1,'last');
            keep = false(size(valid_uv)); keep(i_first:i_last) = true;
            u = u(keep); v = v(keep); time = time(keep);
            N = length(u);
            if N < WIN_MIN*2, continue; end

            WIN_LEN = max(WIN_MIN, min(WIN_MAX, floor(N/3)));
            OVERLAP = floor(WIN_LEN * OVERLAP_FRAC);
            win_inc = max(1, WIN_LEN - OVERLAP);
            n_win   = floor((N - WIN_LEN)/win_inc) + 1;
            if n_win < 2, continue; end

            df     = Fs_cpd/NFFT_PAD;
            f_vec  = (-NFFT_PAD/2 : NFFT_PAD/2-1) * df;
            idx_cw = f_vec < 0;
            f_pos_full = abs(f_vec(idx_cw));
            [f_pos_full, ic_sort] = sort(f_pos_full);

            time_sp = NaT(n_win,1);
            P_cw_sp = nan(length(f_pos_full), n_win);
            w = win_hann(WIN_LEN);
            wnorm = sum(w.^2) * df;

            for j = 1:n_win
                i1 = (j-1)*win_inc + 1; i2 = i1 + WIN_LEN - 1;
                if i2 > N, break; end
                uu = detrend(u(i1:i2),'linear'); uu = uu(:);
                vv = detrend(v(i1:i2),'linear'); vv = vv(:);
                W  = (uu + 1i*vv) .* w;
                Y  = fftshift(fft(W, NFFT_PAD));
                P  = abs(Y).^2 / wnorm;
                P_cw_w = P(idx_cw);
                P_cw_sp(:,j) = P_cw_w(ic_sort);
                time_sp(j)   = time(round((i1+i2)/2));
            end

            col_ok = any(isfinite(P_cw_sp),1) & ~isnat(time_sp(:)');
            if ~any(col_ok), continue; end
            j_first = find(col_ok,1,'first'); j_last = find(col_ok,1,'last');
            kp = false(size(col_ok)); kp(j_first:j_last) = true;
            P_cw_sp = P_cw_sp(:,kp);
            time_sp = time_sp(kp);
            if numel(time_sp) < 2 || numel(f_pos_full) < 2, continue; end

            SP(i).valid   = true;
            SP(i).f_plot  = f_pos_full(:);
            SP(i).time_sp = time_sp;
            SP(i).P_plot  = 10*log10(P_cw_sp + eps);
            SP(i).f_inert = d.f_inertial;
            if ~isempty(d.nike_col) && any(isfinite(d.nike_col))
                [~, ip] = max(d.nike_col);
                SP(i).time_peak = d.time(ip);
            end
            SP(i).xmin = min(time_sp); SP(i).xmax = max(time_sp);
        catch ME
            report_tile_error(sprintf('S4 prep %s', D(i).name), ME);
        end
    end

    all_P = [];
    for i = 1:n_cyc
        if SP(i).valid, all_P = [all_P; SP(i).P_plot(:)]; end %#ok<AGROW>
    end
    all_P(~isfinite(all_P)) = [];
    if ~isempty(all_P), clim_set = prctile(all_P,[5 99]); else, clim_set = [-50 10]; end
    if ~all(isfinite(clim_set)) || clim_set(2) <= clim_set(1), clim_set = [-50 10]; end

    fig = figure('Name','Fig S4 - Inertial Spectrograms', ...
        'Units','normalized','Position',[0.05 0.20 0.95 0.55]);
    tiledlayout(1, n_cyc, 'TileSpacing','compact','Padding','compact');

    for i = 1:n_cyc
        try
            nexttile(i);
            if ~SP(i).valid, title(SP(i).name,'FontSize',12); axis off; continue; end
            s = SP(i);

            pcolor(s.time_sp, s.f_plot, s.P_plot); shading interp;
            if use_cmo, colormap(gca, cmocean('thermal')); else, colormap(gca, hot); end
            safe_clim(clim_set); hold on;

            % f marked white-under-black, per the SI caption
            yline(s.f_inert,  '-','Color','w',     'LineWidth',2.6,'HandleVisibility','off');
            yline(s.f_inert,  '-','Color',c_finert,'LineWidth',1.4,'HandleVisibility','off');
            yline(2*s.f_inert,':','Color',c_2f,    'LineWidth',1.0,'HandleVisibility','off');
            yline(f_M2,      '--','Color',c_M2c,   'LineWidth',1.0,'HandleVisibility','off');

            % BackgroundColor must be an RGB triplet (RGBA errors)
            xoff = s.xmin + (s.xmax - s.xmin)*0.02;
            text(xoff, s.f_inert*1.08, 'f', ...
                'FontSize',10,'FontWeight','bold','Color','w', ...
                'BackgroundColor',[0.15 0.15 0.15],'Margin',1);
            text(xoff, 2*s.f_inert*1.05, '2f', ...
                'FontSize',9,'Color','w','BackgroundColor',[0.25 0.25 0.25],'Margin',1);
            text(xoff, f_M2*1.05, 'M_2', ...
                'FontSize',9,'Color','w','BackgroundColor',[0.25 0.25 0.25],'Margin',1);

            if ~isnat(s.time_peak)
                plot(s.time_peak, 0.07, '^', ...
                    'MarkerEdgeColor','k','MarkerFaceColor',c_peak, ...
                    'MarkerSize',10,'LineWidth',1.2);
                text(s.time_peak, 0.085, 'peak NIKE','HorizontalAlignment','center', ...
                    'FontSize',8,'FontWeight','bold','Color','w');
            end

            safe_xlim(s.xmin, s.xmax);
            set(gca,'YScale','log');
            safe_ylim(0.05, 4);
            set(gca,'YLimMode','manual');
            title(s.name,'FontSize',12);
            set_xax_date(gca);
            xlabel('Date (UTC)','FontSize',10);
            if i == 1
                ylabel({'\bf(a)';'Frequency (cpd)'},'FontSize',11,'FontWeight','bold');
            else
                ylabel('Frequency (cpd)','FontSize',11,'FontWeight','bold');
            end
            if i == n_cyc
                cb = colorbar;
                cb.Label.String     = 'Spectral Power (dB)';
                cb.Label.FontWeight = 'bold'; cb.Label.FontSize = 10;
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S4 tile %d', i), ME);
        end
    end

    ok = saveFigureRobust(fig, OUT_DIR, figFileName('S4'));
    if ok, fprintf(' Saved Fig S4\n'); else, fprintf(2,' Fig S4 was NOT saved\n'); end

    try
        fprintf('\n[FIG S4 - INERTIAL SPECTROGRAMS]\n');
        fprintf('%-8s | %-12s\n','Storm','f_inert(cpd)');
        fprintf('%s\n', repmat('-',1,30));
        for i = 1:n_cyc
            if SP(i).valid
                fprintf('%-8s | %-12.3f\n', SP(i).name, SP(i).f_inert);
            end
        end
    catch ME
        report_tile_error('S4 QC block', ME);
    end
end
%  FIG S5 (draft S7) - VERTICAL PROFILE SNAPSHOTS AT PEAK NIKE
%  writes FigS5_Vertical_Profile_Snapshots.png
function ok = make_figS5_VertProfiles(D, n_cyc, OUT_DIR, beautify, rho0)
    fprintf('\n--- Fig S5 : Vertical profile snapshots ---\n');
    ok = false;

    c_u   = [0.10 0.45 0.75]; c_u_fill   = [0.55 0.75 0.95];
    c_coh = [0.30 0.10 0.50]; c_coh_fill = [0.60 0.40 0.75];
    c_brd = [0.40 0.55 0.20];
    c_mld = [0.95 0.50 0.10];
    c_N2  = [0.05 0.40 0.65]; c_S2 = [0.85 0.30 0.10];
    ZMAX_PLOT = 300;

    tmpl = struct('valid',false,'name','','time_pk',NaT,'z',[],'z_mid',[], ...
        'NIKE_coh_pk',[],'NIKE_brd_pk',[],'U_NI',[],'MLD_pk',NaN, ...
        'N2',[],'S2',[],'Ri',[],'pk_kJm2',NaN);
    S = repmat(tmpl, n_cyc, 1);
    for i = 1:n_cyc, S(i).name = D(i).name; end

    for i = 1:n_cyc
        try
            d = D(i);
            if isempty(d.NIKE_coh_z) || isempty(d.z), continue; end
            iz_p = d.z <= ZMAX_PLOT;
            [pk_val, i_pk] = max(d.nike_col);
            if isempty(pk_val) || ~isfinite(pk_val), continue; end

            NIKE_coh_pk = d.NIKE_coh_z(i_pk, iz_p);
            if ~isempty(d.NIKE_brd_z)
                NIKE_brd_pk = d.NIKE_brd_z(i_pk, iz_p);
            else
                NIKE_brd_pk = nan(size(NIKE_coh_pk));
            end
            MLD_pk = d.MLD(i_pk);
            U_NI   = sqrt(2 * NIKE_coh_pk*1000 / rho0);

            z_mid    = 0.5*(d.z(1:end-1) + d.z(2:end));
            iz_mid_p = z_mid <= ZMAX_PLOT;
            z_mid_p  = z_mid(iz_mid_p);
            if ~isempty(d.N2), N2_pk = d.N2(i_pk,:); N2_pk_p = N2_pk(iz_mid_p); else, N2_pk_p = nan(size(z_mid_p)); end
            if ~isempty(d.S2), S2_pk = d.S2(i_pk,:); S2_pk_p = S2_pk(iz_mid_p); else, S2_pk_p = nan(size(z_mid_p)); end
            Ri_pk = N2_pk_p ./ (S2_pk_p + eps);
            Ri_pk(~isfinite(Ri_pk) | Ri_pk <= 0) = NaN;

            S(i).valid       = true;
            S(i).time_pk     = d.time(i_pk);
            S(i).z           = d.z(iz_p);
            S(i).z_mid       = z_mid_p;
            S(i).NIKE_coh_pk = NIKE_coh_pk(:);
            S(i).NIKE_brd_pk = NIKE_brd_pk(:);
            S(i).U_NI        = U_NI(:);
            S(i).MLD_pk      = MLD_pk;
            S(i).N2          = N2_pk_p(:);
            S(i).S2          = S2_pk_p(:);
            S(i).Ri          = Ri_pk(:);
            S(i).pk_kJm2     = pk_val/1000;
        catch ME
            report_tile_error(sprintf('S5 prep %s', D(i).name), ME);
        end
    end

    all_U = []; all_NIKE = []; all_N2 = []; all_S2 = [];
    for i = 1:n_cyc
        if ~S(i).valid, continue; end
        all_U    = [all_U;    S(i).U_NI];                             %#ok<AGROW>
        all_NIKE = [all_NIKE; S(i).NIKE_coh_pk; S(i).NIKE_brd_pk];    %#ok<AGROW>
        all_N2   = [all_N2;   S(i).N2];                               %#ok<AGROW>
        all_S2   = [all_S2;   S(i).S2];                               %#ok<AGROW>
    end
    all_U    = all_U(isfinite(all_U));
    all_NIKE = all_NIKE(isfinite(all_NIKE));
    if isempty(all_U),    xlim_a = [0 1]; else, xlim_a = [0, max(prctile(all_U,99)*1.15, eps)]; end
    if isempty(all_NIKE), xlim_b = [0 1]; else, xlim_b = [0, max(prctile(all_NIKE,99)*1.15, eps)]; end
    all_N2(all_N2 <= 0 | ~isfinite(all_N2)) = [];
    all_S2(all_S2 <= 0 | ~isfinite(all_S2)) = [];
    if ~isempty(all_N2) && ~isempty(all_S2)
        xlim_c = [10^floor(log10(min([all_N2;all_S2]))), 10^ceil(log10(max([all_N2;all_S2])))];
    else
        xlim_c = [1e-7, 1e-3];
    end

    fig = figure('Name','Fig S5 - Vertical Profile Snapshots', ...
        'Units','normalized','Position',[0.05 0.05 0.92 0.90]);
    tiledlayout(3, n_cyc, 'TileSpacing','compact','Padding','compact');

    % ---- Row (a) : |U_NI| ----
    for i = 1:n_cyc
        try
            nexttile(i);
            if ~S(i).valid, title(S(i).name,'FontSize',11); axis off; continue; end
            s = S(i);
            z_v = s.z; U_v = s.U_NI;
            valid = isfinite(U_v) & isfinite(z_v);
            if any(valid)
                xf = [U_v(valid); zeros(sum(valid),1)];
                yf = [z_v(valid); flipud(z_v(valid))];
                fill(xf, yf, c_u_fill, 'FaceAlpha',0.30,'EdgeColor','none'); hold on;
                plot(U_v(valid), z_v(valid),'Color',c_u,'LineWidth',2.2);
                if isfinite(s.MLD_pk)
                    yline(s.MLD_pk,'-', 'Color','w',  'LineWidth',2.4,'HandleVisibility','off');
                    yline(s.MLD_pk,'--','Color',c_mld,'LineWidth',1.8,'HandleVisibility','off');
                    text(xlim_a(2)*0.95, s.MLD_pk - 8, sprintf('MLD = %.0f m', s.MLD_pk), ...
                        'HorizontalAlignment','right','FontSize',8,'FontWeight','bold', ...
                        'Color',c_mld,'BackgroundColor','w');
                end
            end
            set(gca,'YDir','reverse');
            safe_xlim(xlim_a(1), xlim_a(2)); safe_ylim(0, ZMAX_PLOT);
            title(sprintf('%s\n(peak NIKE: %s)', s.name, datestr(s.time_pk,'mm/dd HH:MM')), ...
                'FontSize',11); %#ok<DATST>
            xlabel('|U_{NI}| (m s^{-1})','FontSize',10);
            if i == 1
                ylabel({'\bf(a)';'Depth (m)'},'FontSize',11,'FontWeight','bold');
            else
                ylabel('Depth (m)','FontSize',11,'FontWeight','bold');
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S5 row a tile %d', i), ME);
        end
    end

    % ---- Row (b) : coherent / broadband NIKE ----
    for i = 1:n_cyc
        try
            nexttile(n_cyc + i);
            if ~S(i).valid, axis off; continue; end
            s = S(i);
            z_v = s.z; coh = s.NIKE_coh_pk; brd = s.NIKE_brd_pk;
            vc = isfinite(coh) & isfinite(z_v);
            vb = isfinite(brd) & isfinite(z_v);
            h_coh = []; h_brd = [];
            if any(vc)
                xf = [coh(vc); zeros(sum(vc),1)];
                yf = [z_v(vc); flipud(z_v(vc))];
                fill(xf,yf,c_coh_fill,'FaceAlpha',0.45,'EdgeColor','none'); hold on;
                h_coh = plot(coh(vc), z_v(vc),'Color',c_coh,'LineWidth',2.5);
            end
            hold on;
            if any(vb)
                h_brd = plot(brd(vb), z_v(vb),'--','Color',c_brd,'LineWidth',1.6);
            end
            if isfinite(s.MLD_pk)
                yline(s.MLD_pk,'-', 'Color','w',  'LineWidth',2.4,'HandleVisibility','off');
                yline(s.MLD_pk,'--','Color',c_mld,'LineWidth',1.8,'HandleVisibility','off');
            end
            % [S13] moved clear of the MLD line, which sits near the top of the
            % panel; the profiles decay to zero below ~100 m so mid-depth is free.
            text(0.96, 0.45, ...
                sprintf('Col. NIKE_{coh} = %.1f kJ/m^2', s.pk_kJm2), ...
                'Units','normalized','HorizontalAlignment','right', ...
                'FontSize',9,'FontWeight','bold','BackgroundColor','w', ...
                'EdgeColor',[0.5 0.5 0.5],'Margin',2);
            set(gca,'YDir','reverse');
            safe_xlim(xlim_b(1), xlim_b(2)); safe_ylim(0, ZMAX_PLOT);
            xlabel('NIKE (kJ m^{-3})','FontSize',10);
            if i == 1
                ylabel({'\bf(b)';'Depth (m)'},'FontSize',11,'FontWeight','bold');
                if ~isempty(h_coh) && ~isempty(h_brd)
                    % [S1] renamed to match the manuscript: Eq. (6) is the
                    % detrended windowed variance, not a bandpassed quantity.
                    legend([h_coh h_brd],{'Coherent NIKE','Windowed KE'}, ...
                        'Box','off','Location','southeast','FontSize',8);
                end
            else
                ylabel('Depth (m)','FontSize',11,'FontWeight','bold');
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S5 row b tile %d', i), ME);
        end
    end

    % ---- Row (c) : N^2 and S^2, with Ri < 0.25 bands ----
    for i = 1:n_cyc
        try
            nexttile(2*n_cyc + i);
            if ~S(i).valid, axis off; continue; end
            s = S(i);
            N2_v = s.N2; S2_v = s.S2;
            N2_v(N2_v <= 0 | ~isfinite(N2_v)) = NaN;
            S2_v(S2_v <= 0 | ~isfinite(S2_v)) = NaN;
            h_N2 = semilogx(N2_v, s.z_mid,'Color',c_N2,'LineWidth',2.2); hold on;
            h_S2 = semilogx(S2_v, s.z_mid,'--','Color',c_S2,'LineWidth',1.8);
            if isfinite(s.MLD_pk)
                yline(s.MLD_pk,'-', 'Color','w',  'LineWidth',2.4,'HandleVisibility','off');
                yline(s.MLD_pk,'--','Color',c_mld,'LineWidth',1.8,'HandleVisibility','off');
            end
            if any(isfinite(s.Ri))
                sub_idx = s.Ri < 0.25 & isfinite(s.Ri);
                if any(sub_idx)
                    sub_z = s.z_mid(sub_idx);
                    for kk = 1:length(sub_z)
                        % RGB only - RGBA errors on several releases
                        yline(sub_z(kk),'-','Color',[0.75 0.75 0.75], ...
                            'LineWidth',0.8,'HandleVisibility','off');
                    end
                    % [S14] was pinned in data units and ran off the right edge.
                    text(0.97, 0.06, sprintf('Ri<0.25: %d levels', sum(sub_idx)), ...
                        'Units','normalized','HorizontalAlignment','right', ...
                        'VerticalAlignment','top', ...
                        'FontSize',8,'Color',[0.4 0.4 0.4],'FontWeight','bold');
                end
            end
            set(gca,'YDir','reverse','XScale','log');
            safe_xlim(xlim_c(1), xlim_c(2)); safe_ylim(0, ZMAX_PLOT);
            xlabel('N^2, S^2 (s^{-2})','FontSize',10);
            if i == 1
                ylabel({'\bf(c)';'Depth (m)'},'FontSize',11,'FontWeight','bold');
                legend([h_N2 h_S2],{'N^2 stratification','S^2 shear'}, ...
                    'Box','off','Location','southeast','FontSize',8);
            else
                ylabel('Depth (m)','FontSize',11,'FontWeight','bold');
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S5 row c tile %d', i), ME);
        end
    end

    ok = saveFigureRobust(fig, OUT_DIR, figFileName('S5'));
    if ok, fprintf(' Saved Fig S5\n'); else, fprintf(2,' Fig S5 was NOT saved\n'); end

    try
        fprintf('\n[FIG S5 - VERTICAL PROFILES AT PEAK NIKE]\n');
        fprintf('%-8s | %-15s | %-9s | %-12s\n','Storm','PeakTime','MLD(m)','PkColNIKE(kJ/m2)');
        fprintf('%s\n', repmat('-',1,60));
        for i = 1:n_cyc
            if ~S(i).valid, continue; end
            s = S(i);
            fprintf('%-8s | %-15s | %-9.1f | %-12.2f\n', ...
                s.name, datestr(s.time_pk,'mm/dd HH:MM'), s.MLD_pk, s.pk_kJm2); %#ok<DATST>
        end
    catch ME
        report_tile_error('S5 QC block', ME);
    end
end
%  FIG S6 (draft S8) - R/L DEPTH-RESOLVED COHERENT NIKE
%  writes FigS6_RL_Depth_Resolved.png
function ok = make_figS6_RLDepth(D, n_cyc, OUT_DIR, beautify)
    fprintf('\n--- Fig S6 : R/L depth-resolved ---\n');
    ok = false;
    use_cmo = exist('cmocean.m','file') == 2;
    ZMAX_PLOT = 300;

    tmpl = struct('valid',false,'name','','time',datetime.empty(0,1),'keep',[], ...
        'z',[],'NIKE_R',[],'NIKE_L',[],'MLD',[],'xmin',NaT,'xmax',NaT);
    S = repmat(tmpl, n_cyc, 1);
    for i = 1:n_cyc, S(i).name = D(i).name; end

    for i = 1:n_cyc
        try
            d = D(i);
            if isempty(d.NIKE_coh_z_R) || isempty(d.NIKE_coh_z_L) || isempty(d.z), continue; end
            iz_p   = d.z <= ZMAX_PLOT;
            R_full = d.NIKE_coh_z_R(:, iz_p);
            L_full = d.NIKE_coh_z_L(:, iz_p);

            valid_t = any(isfinite(R_full),2) & any(isfinite(L_full),2);
            if ~any(valid_t), continue; end
            i_first = find(valid_t,1,'first'); i_last = find(valid_t,1,'last');
            keep = false(size(valid_t)); keep(i_first:i_last) = true;

            R_full(keep,:) = fillmissing(R_full(keep,:),'nearest',1);
            L_full(keep,:) = fillmissing(L_full(keep,:),'nearest',1);

            z_p    = d.z(iz_p);
            col_ok = any(isfinite(R_full(keep,:)),1) & any(isfinite(L_full(keep,:)),1);
            if sum(col_ok) < 2 || sum(keep) < 2, continue; end

            S(i).valid  = true;
            S(i).time   = d.time;
            S(i).keep   = keep;
            S(i).z      = z_p(col_ok);
            S(i).NIKE_R = R_full(:, col_ok);
            S(i).NIKE_L = L_full(:, col_ok);
            S(i).MLD    = d.MLD;
            S(i).xmin   = min(d.time(keep));
            S(i).xmax   = max(d.time(keep));
        catch ME
            report_tile_error(sprintf('S6 prep %s', D(i).name), ME);
        end
    end

    all_data = [];
    for i = 1:n_cyc
        if ~S(i).valid, continue; end
        k = S(i).keep;
        all_data = [all_data; S(i).NIKE_R(k,:); S(i).NIKE_L(k,:)]; %#ok<AGROW>
    end
    all_data = all_data(:);
    all_data(~isfinite(all_data) | all_data <= 0) = [];
    if ~isempty(all_data), clim_set = [0, prctile(all_data,75)]; else, clim_set = [0 1]; end
    if ~isfinite(clim_set(2)) || clim_set(2) <= 0, clim_set(2) = 1; end

    fig = figure('Name','Fig S6 - R/L Depth-Resolved NIKE', ...
        'Units','normalized','Position',[0.02 0.05 0.98 0.88]);
    tiledlayout(2, n_cyc, 'TileSpacing','compact','Padding','compact');

    % ---- Row (a) : right of track ----
    for i = 1:n_cyc
        try
            nexttile(i);
            if ~S(i).valid, title(S(i).name,'FontSize',11); axis off; continue; end
            s = S(i); k = s.keep; time_k = s.time(k);

            pcolor(time_k, s.z, s.NIKE_R(k,:)'); shading interp;
            set(gca,'YDir','reverse');
            if use_cmo, colormap(gca, cmocean('thermal')); else, colormap(gca, parula); end
            safe_clim(clim_set); hold on;
            plot(time_k, s.MLD(k),'w-','LineWidth',2.4);
            plot(time_k, s.MLD(k),'k--','LineWidth',1.2);

            safe_xlim(s.xmin, s.xmax); safe_ylim(0, ZMAX_PLOT);
            title(sprintf('%s : Right of track', s.name),'FontSize',11);
            set_xax_date(gca);
            xlabel('Date (UTC)','FontSize',10);
            if i == 1
                ylabel({'\bf(a)';'Depth (m)';'NIKE_{coh,right}'},'FontSize',11,'FontWeight','bold');
            else
                ylabel('Depth (m)','FontSize',11,'FontWeight','bold');
            end
            if i == n_cyc
                cb = colorbar;
                cb.Label.String     = 'NIKE_{coh} (kJ m^{-3})';
                cb.Label.FontWeight = 'bold'; cb.Label.FontSize = 10;
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S6 row a tile %d', i), ME);
        end
    end

    % ---- Row (b) : left of track ----
    for i = 1:n_cyc
        try
            nexttile(n_cyc + i);
            if ~S(i).valid, axis off; continue; end
            s = S(i); k = s.keep; time_k = s.time(k);

            pcolor(time_k, s.z, s.NIKE_L(k,:)'); shading interp;
            set(gca,'YDir','reverse');
            if use_cmo, colormap(gca, cmocean('thermal')); else, colormap(gca, hot); end
            safe_clim(clim_set); hold on;
            plot(time_k, s.MLD(k),'w-','LineWidth',2.4);
            plot(time_k, s.MLD(k),'k--','LineWidth',1.2);

            safe_xlim(s.xmin, s.xmax); safe_ylim(0, ZMAX_PLOT);
            title(sprintf('%s : Left of track', s.name),'FontSize',11);
            set_xax_date(gca);
            xlabel('Date (UTC)','FontSize',10);
            if i == 1
                ylabel({'\bf(b)';'Depth (m)';'NIKE_{coh,left}'},'FontSize',11,'FontWeight','bold');
                h_dummy = plot(NaN,NaN,'k--','LineWidth',1.4);
                legend(h_dummy,{'MLD'},'Box','off','Location','southwest','FontSize',8);
            else
                ylabel('Depth (m)','FontSize',11,'FontWeight','bold');
            end
            if i == n_cyc
                cb = colorbar;
                cb.Label.String     = 'NIKE_{coh} (kJ m^{-3})';
                cb.Label.FontWeight = 'bold'; cb.Label.FontSize = 10;
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S6 row b tile %d', i), ME);
        end
    end

    ok = saveFigureRobust(fig, OUT_DIR, figFileName('S6'));
    if ok, fprintf(' Saved Fig S6\n'); else, fprintf(2,' Fig S6 was NOT saved\n'); end

    try
        fprintf('\n[FIG S6 - R/L DEPTH-RESOLVED]\n');
        fprintf('%-8s | %-10s | %-10s | %-10s | %-10s | %-9s | %-9s\n', ...
            'Storm','MeanR_ML','MeanR_Dp','MeanL_ML','MeanL_Dp','R/L_ML','R/L_Dp');
        fprintf('%s\n', repmat('-',1,80));
        for i = 1:n_cyc
            if ~S(i).valid, continue; end
            s = S(i); k = s.keep;
            z_v = s.z; R_mat = s.NIKE_R(k,:); L_mat = s.NIKE_L(k,:); mld_k = s.MLD(k);
            nrows = size(R_mat,1);
            R_ml = nan(nrows,1); R_dp = nan(nrows,1);
            L_ml = nan(nrows,1); L_dp = nan(nrows,1);
            for it = 1:nrows
                if isnan(mld_k(it)), continue; end
                is_ml = z_v <= mld_k(it);
                is_dp = z_v > mld_k(it) & z_v <= 200;
                if any(is_ml)
                    R_ml(it) = mean(R_mat(it,is_ml),'omitnan');
                    L_ml(it) = mean(L_mat(it,is_ml),'omitnan');
                end
                if any(is_dp)
                    R_dp(it) = mean(R_mat(it,is_dp),'omitnan');
                    L_dp(it) = mean(L_mat(it,is_dp),'omitnan');
                end
            end
            mR_ml = mean(R_ml,'omitnan'); mR_dp = mean(R_dp,'omitnan');
            mL_ml = mean(L_ml,'omitnan'); mL_dp = mean(L_dp,'omitnan');
            fprintf('%-8s | %-10.3f | %-10.3f | %-10.3f | %-10.3f | %-9.2f | %-9.2f\n', ...
                s.name, mR_ml, mR_dp, mL_ml, mL_dp, mR_ml/(mL_ml+eps), mR_dp/(mL_dp+eps));
        end
    catch ME
        report_tile_error('S6 QC block', ME);
    end
end
%  FIG S7 (draft S9) - R_NI SENSITIVITY (100 / 150 / 200 km)
%  writes FigS7_RNI_Sensitivity.png
function ok = make_figS7_RNISensitivity(D, n_cyc, OUT_DIR, beautify, HYCOM_PATHS, Hmax, rho0)
    fprintf('\n--- Fig S7 : R_NI sensitivity (slow, re-reads HYCOM) ---\n');
    ok = false;

    c_r100 = [0.20 0.45 0.70]; c_r150 = [0.30 0.55 0.25]; c_r200 = [0.85 0.20 0.20];
    c_radii = [c_r100; c_r150; c_r200];
    R_NI_TEST = [100 150 200]; n_RNI = length(R_NI_TEST);

    tmpl = struct('valid',false,'name','','time',datetime.empty(0,1),'keep',[], ...
        'NIKE_brd',[],'NIKE_brd_raw',[],'AR',[],'xmin',NaT,'xmax',NaT);
    S = repmat(tmpl, n_cyc, 1);
    for i = 1:n_cyc, S(i).name = D(i).name; end

    for i = 1:n_cyc
        storm = D(i).name;
        try
            if ~isfield(HYCOM_PATHS,storm), warning('No HYCOM path %s',storm); continue; end
            hycom_path = HYCOM_PATHS.(storm);
            if ~exist(hycom_path,'dir'), warning('HYCOM dir missing: %s',hycom_path); continue; end
            fprintf('  %s ... ', storm);

            T = D(i).T; time = D(i).time; lat_tc = D(i).lat; lon_tc = D(i).lon;
            nt = height(T);

            heading_rad = nan(nt,1);
            for it = 1:nt
                if it == 1
                    dl = lat_tc(2)-lat_tc(1);       dlo = lon_tc(2)-lon_tc(1);
                elseif it == nt
                    dl = lat_tc(nt)-lat_tc(nt-1);   dlo = lon_tc(nt)-lon_tc(nt-1);
                else
                    dl = lat_tc(it+1)-lat_tc(it-1); dlo = lon_tc(it+1)-lon_tc(it-1);
                end
                dy = dl*111; dx = dlo*111*cosd(lat_tc(it));
                heading_rad(it) = atan2(dy,dx);
            end
            heading_rad = unwrap(heading_rad);
            heading_rad = movmean(heading_rad,3,'omitnan');

            [uvfiles, uvtimes, uvfi, uvti] = index_hycom_dir(hycom_path);
            if isempty(uvfiles), warning('No .nc in %s', hycom_path); fprintf('skipped\n'); continue; end
            fn_ref = uvfiles{1};
            try
                lon_g = ncread(fn_ref,'lon'); lat_g = ncread(fn_ref,'lat');
            catch
                lon_g = ncread(fn_ref,'longitude'); lat_g = ncread(fn_ref,'latitude');
            end
            depth = ncread(fn_ref,'depth'); idz = find(depth <= Hmax);
            z = depth(idz); nz = length(z);

            NIKE_brd_R = nan(nt,n_RNI);
            NIKE_R_R   = nan(nt,n_RNI);
            NIKE_L_R   = nan(nt,n_RNI);

            for it = 1:nt
                if isnan(lat_tc(it)) || isnan(lon_tc(it)), continue; end
                [dd,ki] = min(abs(uvtimes - time(it)));
                if isempty(ki) || dd > hours(2), continue; end
                fn = uvfiles{uvfi(ki)}; tt = uvti(ki);
                ix = find(lon_g >= lon_tc(it)-2 & lon_g <= lon_tc(it)+2);
                iy = find(lat_g >= lat_tc(it)-2 & lat_g <= lat_tc(it)+2);
                if isempty(ix) || isempty(iy), continue; end
                try
                    u_box = ncread(fn,'water_u',[ix(1) iy(1) 1 tt],[length(ix) length(iy) nz 1]);
                    v_box = ncread(fn,'water_v',[ix(1) iy(1) 1 tt],[length(ix) length(iy) nz 1]);
                    u_box = squeeze(u_box); v_box = squeeze(v_box);
                catch
                    continue;
                end
                [LonSub,LatSub] = meshgrid(lon_g(ix),lat_g(iy));
                dist  = dist_km(lat_tc(it), lon_tc(it), LatSub, LonSub);
                h_rad = heading_rad(it);
                dy_pt = (LatSub - lat_tc(it))*111;
                dx_pt = (LonSub - lon_tc(it))*111*cosd(lat_tc(it));
                cross_z = dx_pt*sin(h_rad) - dy_pt*cos(h_rad);

                for ir = 1:n_RNI
                    R_NI = R_NI_TEST(ir);
                    mask = double(dist <= R_NI); mask(mask==0) = NaN;
                    mask_R = mask; mask_R(cross_z <  0) = NaN;
                    mask_L = mask; mask_L(cross_z >= 0) = NaN;
                    u_pf = nan(nz,1); v_pf = nan(nz,1);
                    u_pR = nan(nz,1); v_pR = nan(nz,1);
                    u_pL = nan(nz,1); v_pL = nan(nz,1);
                    for kk = 1:nz
                        sU = u_box(:,:,kk); sV = v_box(:,:,kk);
                        if ~isequal(size(sU),size(mask)), sU = sU'; sV = sV'; end
                        u_pf(kk) = mean(sU.*mask,  'all','omitnan'); v_pf(kk) = mean(sV.*mask,  'all','omitnan');
                        u_pR(kk) = mean(sU.*mask_R,'all','omitnan'); v_pR(kk) = mean(sV.*mask_R,'all','omitnan');
                        u_pL(kk) = mean(sU.*mask_L,'all','omitnan'); v_pL(kk) = mean(sV.*mask_L,'all','omitnan');
                    end
                    Nf = 0.5*rho0*(u_pf.^2+v_pf.^2);
                    NR = 0.5*rho0*(u_pR.^2+v_pR.^2);
                    NL = 0.5*rho0*(u_pL.^2+v_pL.^2);
                    vf = isfinite(Nf); vR = isfinite(NR); vL = isfinite(NL);
                    if sum(vf) >= 2, NIKE_brd_R(it,ir) = trapz(z(vf),Nf(vf)); end
                    if sum(vR) >= 2, NIKE_R_R(it,ir)   = trapz(z(vR),NR(vR)); end
                    if sum(vL) >= 2, NIKE_L_R(it,ir)   = trapz(z(vL),NL(vL)); end
                end
            end
            NIKE_brd_R = NIKE_brd_R/1000; NIKE_R_R = NIKE_R_R/1000; NIKE_L_R = NIKE_L_R/1000;

            AR_R = nan(nt,n_RNI);
            for ir = 1:n_RNI
                for it = 1:nt
                    Rv = NIKE_R_R(it,ir); Lv = NIKE_L_R(it,ir); tot = NIKE_brd_R(it,ir);
                    if isfinite(Rv) && isfinite(Lv) && isfinite(tot) && tot > 0 ...
                            && Rv > 0.10*tot && Lv > 0.10*tot
                        AR_R(it,ir) = min(Rv/Lv, 10);
                    end
                end
            end

            valid_any = any(isfinite(NIKE_brd_R),2);
            if ~any(valid_any), fprintf('no data\n'); continue; end
            i_first = find(valid_any,1,'first'); i_last = find(valid_any,1,'last');
            keep = false(size(valid_any)); keep(i_first:i_last) = true;

            NIKE_brd_R_p = NIKE_brd_R;
            for ir = 1:n_RNI
                NIKE_brd_R_p(keep,ir) = fillmissing(NIKE_brd_R(keep,ir),'linear','EndValues','none');
            end

            S(i).valid        = true;
            S(i).time         = time;
            S(i).keep         = keep;
            S(i).NIKE_brd     = NIKE_brd_R_p;
            S(i).NIKE_brd_raw = NIKE_brd_R;
            S(i).AR           = AR_R;
            S(i).xmin         = min(time(keep));
            S(i).xmax         = max(time(keep));
            fprintf('done\n');
        catch ME
            fprintf(2,'failed (%s)\n', ME.message);
        end
    end

    all_brd = [];
    for i = 1:n_cyc
        if S(i).valid
            k = S(i).keep;
            all_brd = [all_brd; S(i).NIKE_brd(k,:)]; %#ok<AGROW>
        end
    end
    all_brd = all_brd(isfinite(all_brd));
    if isempty(all_brd), ylim_a = [0 30];
    else, ylim_a = [0, max(prctile(all_brd(:),99)*1.15, 1)]; end

    AR_means = nan(n_cyc,n_RNI); AR_stds = nan(n_cyc,n_RNI);
    for i = 1:n_cyc
        if ~S(i).valid, continue; end
        k = S(i).keep;
        for ir = 1:n_RNI
            ar_v = S(i).AR(k,ir); ar_v = ar_v(isfinite(ar_v));
            if ~isempty(ar_v)
                AR_means(i,ir) = mean(ar_v);
                AR_stds(i,ir)  = std(ar_v);
            end
        end
    end
    ylim_b = [0, max(AR_means(:)+AR_stds(:),[],'omitnan')*1.15];
    if ~isfinite(ylim_b(2)) || ylim_b(2) < 4, ylim_b(2) = 4; end

    fig = figure('Name','Fig S7 - R_NI Sensitivity', ...
        'Units','normalized','Position',[0.05 0.10 0.92 0.78]);
    tiledlayout(2, n_cyc, 'TileSpacing','compact','Padding','compact');

    % ---- Row (a) : broadband NIKE per radius ----
    for i = 1:n_cyc
        try
            nexttile(i);
            if ~S(i).valid, title(S(i).name,'FontSize',12); axis off; continue; end
            s = S(i); k = s.keep; tk = s.time(k);

            order   = [3 1 2];
            h_lines = gobjects(n_RNI,1);
            for jj = 1:length(order)
                ir = order(jj);
                y_k = s.NIKE_brd(k,ir);
                valid = isfinite(y_k);
                if ~any(valid), continue; end
                if ir == 2
                    h_lines(ir) = plot(tk(valid), y_k(valid),     'Color',c_radii(ir,:),'LineWidth',2.6);
                elseif ir == 1
                    h_lines(ir) = plot(tk(valid), y_k(valid),'--','Color',c_radii(ir,:),'LineWidth',1.8);
                else
                    h_lines(ir) = plot(tk(valid), y_k(valid),':', 'Color',c_radii(ir,:),'LineWidth',1.8);
                end
                hold on;
            end

            peak_vals = max(s.NIKE_brd_raw(k,:),[],1,'omitnan');
            if all(isfinite(peak_vals)) && mean(peak_vals) > 0
                sens_pct = 100*(max(peak_vals)-min(peak_vals))/mean(peak_vals);
                text(s.xmin + 0.02*(s.xmax-s.xmin), ylim_a(2)*0.88, ...
                    sprintf('Peak spread: %.1f%%', sens_pct), ...
                    'FontSize',9,'FontWeight','bold','BackgroundColor','w', ...
                    'EdgeColor',[0.5 0.5 0.5],'Margin',2);
            end

            safe_xlim(s.xmin, s.xmax); safe_ylim(ylim_a(1), ylim_a(2));
            title(s.name,'FontSize',12);
            set_xax_date(gca);
            xlabel('Date (UTC)','FontSize',10);
            if i == 1
                % [S7] This is the column-integrated TOTAL kinetic energy of the
                % disk-mean velocity profile: no detrending, no harmonic fit. It
                % is NOT the coherent NIKE of Table 2 and must not be labelled as
                % such. The figure tests sensitivity to R_NI, for which a simpler
                % metric applied identically at all three radii is sufficient.
                ylabel({'\bf(a)';'Column KE, disk-mean';'(kJ m^{-2})'},'FontSize',11,'FontWeight','bold');
                lbl_all = {'R_{NI}=100 km','R_{NI}=150 km (paper)','R_{NI}=200 km'};
                sel = isgraphics(h_lines);
                if any(sel)
                    legend(h_lines(sel), lbl_all(sel),'Box','off','Location','northeast','FontSize',8);
                end
            else
                ylabel('Column KE, disk-mean (kJ m^{-2})','FontSize',11,'FontWeight','bold');
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S7 row a tile %d', i), ME);
        end
    end

    % ---- Row (b) : mean R/L per radius ----
    for i = 1:n_cyc
        try
            nexttile(n_cyc + i);
            if ~S(i).valid, axis off; continue; end
            means_i = AR_means(i,:); stds_i = AR_stds(i,:);
            stds_i(~isfinite(stds_i)) = 0;

            b = bar(1:n_RNI, means_i, 0.55,'FaceColor','flat','EdgeColor','k','LineWidth',1.2);
            b.CData = c_radii; hold on;
            errorbar(1:n_RNI, means_i, stds_i,'k','LineStyle','none','LineWidth',1.4,'CapSize',8);
            rectangle('Position',[1.65, 0, 0.7, max(ylim_b(2)*0.95, eps)], ...
                'EdgeColor',[0.30 0.55 0.25],'LineWidth',2.5,'LineStyle','--');
            for ir = 1:n_RNI
                if isfinite(means_i(ir))
                    text(ir, means_i(ir)+max(stds_i(ir),0.1)+0.15, sprintf('%.2f', means_i(ir)), ...
                        'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
                end
            end
            if all(isfinite(means_i)) && mean(means_i) > 0
                sens_RL = 100*(max(means_i)-min(means_i))/mean(means_i);
                text(0.6, ylim_b(2)*0.92, sprintf('R/L spread: %.1f%%', sens_RL), ...
                    'FontSize',8,'FontWeight','bold','BackgroundColor','w', ...
                    'EdgeColor',[0.5 0.5 0.5],'Margin',2);
            end
            yline(1.5,'--','lit. lo (1.5)','Color',[0.4 0.4 0.4],'LineWidth',1.0, ...
                'FontSize',7,'HandleVisibility','off');
            yline(3.0,'--','lit. hi (3.0)','Color',[0.4 0.4 0.4],'LineWidth',1.0, ...
                'FontSize',7,'HandleVisibility','off');
            set(gca,'XTick',1:n_RNI,'XTickLabel',{'100 km','150 km','200 km'});
            safe_ylim(ylim_b(1), ylim_b(2));
            if i == 1
                ylabel({'\bf(b)';'Mean R/L (column KE)'},'FontSize',11,'FontWeight','bold');
            else
                ylabel('Mean R/L (column KE)','FontSize',11,'FontWeight','bold');
            end
            xlabel('R_{NI} (km)','FontSize',10);
            beautify();
        catch ME
            report_tile_error(sprintf('S7 row b tile %d', i), ME);
        end
    end

    ok = saveFigureRobust(fig, OUT_DIR, figFileName('S7'));
    if ok, fprintf(' Saved Fig S7\n'); else, fprintf(2,' Fig S7 was NOT saved\n'); end

    try
        fprintf('\n[FIG S7 - R_NI SENSITIVITY]\n');
        fprintf('  NOTE: mean R/L of column KE on the disk-mean profile.\n');
        fprintf('  Not the median coherent R/L of Table 2; use for radius sensitivity only.\n');
        fprintf('%-8s | %-12s | %-12s | %-12s\n','Storm','R/L@100','R/L@150','R/L@200');
        fprintf('%s\n', repmat('-',1,60));
        for i = 1:n_cyc
            if ~S(i).valid, continue; end
            fprintf('%-8s | %-12.2f | %-12.2f | %-12.2f\n', ...
                S(i).name, AR_means(i,1), AR_means(i,2), AR_means(i,3));
        end
    catch ME
        report_tile_error('S7 QC block', ME);
    end
end
%  FIG S8 (draft S10) - N^2-S^2 PHASE SPACE
%  writes FigS8_PhaseSpace_N2_S2.png
function ok = make_figS8_PhaseSpace(D, n_cyc, OUT_DIR, beautify)
    fprintf('\n--- Fig S8 : N^2-S^2 phase space ---\n');
    ok = false;

    c_Ri025 = [0.85 0.20 0.20];
    c_Ri05  = [0.55 0.30 0.65];
    c_Ri10  = [0.10 0.55 0.30];
    LOG_N2_RANGE = [-7,-3]; LOG_S2_RANGE = [-9,-3]; N_BINS = 50;
    use_cmo = exist('cmocean.m','file') == 2;

    tmpl = struct('valid',false,'name','','logN2',[],'logS2',[],'Ri_v',[], ...
        'n_total',0,'n_subKH',0,'n_diff',0,'n_marg',0,'n_stable',0);
    S = repmat(tmpl, n_cyc, 1);
    for i = 1:n_cyc, S(i).name = D(i).name; end

    for i = 1:n_cyc
        try
            d = D(i);
            if isempty(d.N2) || isempty(d.S2), continue; end
            valid_t = any(isfinite(d.N2),2) | any(isfinite(d.S2),2);
            if ~any(valid_t), continue; end
            i_first = find(valid_t,1,'first'); i_last = find(valid_t,1,'last');
            keep = false(size(valid_t)); keep(i_first:i_last) = true;
            N2_use = d.N2(keep,:); S2_use = d.S2(keep,:);
            N2_v = N2_use(:); S2_v = S2_use(:);
            % exclude the numerical S^2 floor (unresolved, not weak, shear)
            mask_ok = isfinite(N2_v) & isfinite(S2_v) & N2_v > 1e-8 & S2_v > 1.1e-7;
            N2_v = N2_v(mask_ok); S2_v = S2_v(mask_ok);
            if isempty(N2_v), continue; end

            Ri_v = N2_v ./ S2_v;
            S(i).valid    = true;
            S(i).logN2    = log10(N2_v);
            S(i).logS2    = log10(S2_v);
            S(i).Ri_v     = Ri_v;
            S(i).n_total  = length(Ri_v);
            S(i).n_subKH  = sum(Ri_v < 0.25);
            S(i).n_diff   = sum(Ri_v >= 0.25 & Ri_v < 0.5);
            S(i).n_marg   = sum(Ri_v >= 0.5  & Ri_v < 1.0);
            S(i).n_stable = sum(Ri_v >= 1.0);
        catch ME
            report_tile_error(sprintf('S8 prep %s', D(i).name), ME);
        end
    end

    edges_N2 = linspace(LOG_N2_RANGE(1),LOG_N2_RANGE(2),N_BINS+1);
    edges_S2 = linspace(LOG_S2_RANGE(1),LOG_S2_RANGE(2),N_BINS+1);

    hist_data = cell(n_cyc,1); all_max = 0;
    for i = 1:n_cyc
        if ~S(i).valid, continue; end
        counts = histcounts2(S(i).logN2, S(i).logS2, edges_N2, edges_S2);
        hist_data{i} = counts; all_max = max(all_max, max(counts(:)));
    end
    nonzero_counts = [];
    for i = 1:n_cyc
        if isempty(hist_data{i}), continue; end
        c = hist_data{i}(:); nonzero_counts = [nonzero_counts; c(c>0)]; %#ok<AGROW>
    end
    if ~isempty(nonzero_counts), clim_set = [0, prctile(nonzero_counts,95)];
    else, clim_set = [0, max(all_max,1)]; end
    if ~isfinite(clim_set(2)) || clim_set(2) <= 0, clim_set(2) = 1; end

    fig = figure('Name','Fig S8 - N^2 S^2 Phase Space', ...
        'Units','normalized','Position',[0.05 0.20 0.95 0.55]);
    tiledlayout(1, n_cyc, 'TileSpacing','compact','Padding','compact');

    for i = 1:n_cyc
        try
            nexttile(i);
            if ~S(i).valid, title(S(i).name,'FontSize',12); axis off; continue; end
            s = S(i);

            histogram2(s.logN2, s.logS2, edges_N2, edges_S2, ...
                'DisplayStyle','tile','ShowEmptyBins','off','EdgeColor','none');
            hold on;
            if use_cmo, colormap(gca, cmocean('thermal')); else, colormap(gca, hot); end
            safe_clim(clim_set);

            x_line  = linspace(LOG_N2_RANGE(1),LOG_N2_RANGE(2),200);
            h_Ri025 = plot(x_line, x_line - log10(0.25),'--','Color',c_Ri025,'LineWidth',2.0);
            h_Ri05  = plot(x_line, x_line - log10(0.5), '-.','Color',c_Ri05, 'LineWidth',1.6);
            h_Ri10  = plot(x_line, x_line - log10(1.0), ':', 'Color',c_Ri10, 'LineWidth',2.0);

            % RGB only - RGBA background errors on several releases
            text(LOG_N2_RANGE(1)+0.2, LOG_S2_RANGE(2)-0.4, 'Ri < 0.25 (KH-active)', ...
                'FontSize',9,'FontWeight','bold','Color',c_Ri025,'BackgroundColor','w','Margin',1);
            text(LOG_N2_RANGE(2)-1.5, LOG_S2_RANGE(1)+0.4, 'Ri > 1 (stable)', ...
                'FontSize',9,'FontWeight','bold','Color',c_Ri10,'BackgroundColor','w','Margin',1);
            if s.n_total > 0
                pct_subKH = 100 * s.n_subKH / s.n_total;
                text(LOG_N2_RANGE(2)-0.6, LOG_S2_RANGE(2)-0.4, ...
                    sprintf('Ri<0.25: %.1f%%', pct_subKH), ...
                    'FontSize',9,'FontWeight','bold','BackgroundColor','w', ...
                    'EdgeColor',[0.5 0.5 0.5],'Margin',2,'HorizontalAlignment','right');
            end

            safe_xlim(LOG_N2_RANGE(1), LOG_N2_RANGE(2));
            safe_ylim(LOG_S2_RANGE(1), LOG_S2_RANGE(2));
            title(s.name,'FontSize',12);
            xlabel('log_{10} N^2 (s^{-2})','FontSize',10);
            if i == 1
                ylabel({'\bf(a)';'log_{10} S^2 (s^{-2})'},'FontSize',11,'FontWeight','bold');
                legend([h_Ri025 h_Ri05 h_Ri10], ...
                    {'Ri = 0.25 (KH instability)','Ri = 0.50 (enhanced mixing)', ...
                     'Ri = 1.00 (marginal stability)'}, ...
                    'Box','off','Location','southwest','FontSize',7);
                % [S12] southwest, then nudged down-left into the empty corner
                % below log10(S^2) = -7.5. northwest put it on top of both the
                % 'Ri < 0.25 (KH-active)' label and the percentage box.
                lg = findobj(gcf,'Type','Legend'); lg = lg(1);
                lg.Units = 'normalized';
                lg.Position(2) = lg.Position(2) - 0.04;
            else
                ylabel('log_{10} S^2 (s^{-2})','FontSize',11,'FontWeight','bold');
            end
            if i == n_cyc
                cb = colorbar;
                cb.Label.String     = 'Bin count (data density)';
                cb.Label.FontWeight = 'bold'; cb.Label.FontSize = 10;
            end
            beautify();
        catch ME
            report_tile_error(sprintf('S8 tile %d', i), ME);
        end
    end

    ok = saveFigureRobust(fig, OUT_DIR, figFileName('S8'));
    if ok, fprintf(' Saved Fig S8\n'); else, fprintf(2,' Fig S8 was NOT saved\n'); end

    try
        fprintf('\n[FIG S8 - PHASE SPACE N^2 vs S^2]\n');
        fprintf('%-8s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s\n', ...
            'Storm','Total#','Ri<0.25','Ri.25-.5','Ri.5-1','Ri>1','%KH');
        fprintf('%s\n', repmat('-',1,80));
        for i = 1:n_cyc
            if ~S(i).valid, continue; end
            s = S(i);
            if s.n_total > 0, pct = 100*s.n_subKH/s.n_total; else, pct = NaN; end
            fprintf('%-8s | %-10d | %-10d | %-10d | %-10d | %-10d | %-10.1f\n', ...
                s.name, s.n_total, s.n_subKH, s.n_diff, s.n_marg, s.n_stable, pct);
        end
    catch ME
        report_tile_error('S8 QC block', ME);
    end
end
%  COMPACT DATE AXIS
%  Avoids linspace() on datetimes and is fully wrapped, so a tick
%  formatting quirk can never abort a figure.
function set_xax_date(ax)
    try
        xl = ax.XLim;
        if numel(xl) < 2, return; end

        if isdatetime(xl)
            span = xl(2) - xl(1);
            tv = xl(1) + (0:3)/3 * span;      % duration arithmetic
            try
                ax.XAxis.TickValues      = tv;
                ax.XAxis.TickLabelFormat = 'dd/MM';
            catch
            end
        else
            tv = linspace(xl(1), xl(2), 4);
            ax.XTick = tv;
            try, datetick(ax,'x','dd/mm','keepticks','keeplimits'); catch, end
        end

        yr_str = datestr(tv(end),'yyyy'); %#ok<DATST>
        ax.XTickLabelRotation = 0;
        ax.FontSize = 8;
        text(ax, 0.98, -0.18, yr_str, ...
            'Units','normalized','HorizontalAlignment','right', ...
            'VerticalAlignment','top','FontSize',8, ...
            'FontAngle','italic','Color',[0.35 0.35 0.35]);
    catch
        % axis cosmetics must never kill a figure
    end
end
%  LOADER HELPERS (v7)
function v = getcol(M,f,n)
    if isfield(M,f) && numel(M.(f))==n, v = M.(f)(:); else, v = nan(n,1); end
end

function A = fld(M,f)
    if isfield(M,f), A = M.(f); else, A = []; end
end

function A = fldz(M,f)
    if isfield(M,f), A = M.(f)/1000; else, A = []; end
end

function v = colv(x)
    if isempty(x), v = []; else, v = x(:); end
end
%  HYCOM directory time-indexer (v7)
%  Builds a lookup of every timestep across all .nc files in a folder,
%  regardless of naming. 'files' is per-row, so fidx is just the row
%  index (the old placeholder line appended a wrong-length vector).
function [files, times, fidx, tidx] = index_hycom_dir(hycom_path)
    files = {}; times = datetime.empty(0,1); fidx = []; tidx = [];
    d = dir(fullfile(hycom_path,'*.nc'));
    if isempty(d), return; end

    flist = cell(numel(d),1);
    for i = 1:numel(d), flist{i} = fullfile(d(i).folder, d(i).name); end

    for i = 1:numel(flist)
        fn = flist{i};
        try
            tv = ncread(fn,'time');
        catch
            continue;
        end
        try
            info = ncinfo(fn,'time');
        catch
            continue;
        end
        un = '';
        for a = 1:numel(info.Attributes)
            if strcmpi(info.Attributes(a).Name,'units'), un = info.Attributes(a).Value; end
        end
        tok = regexp(un,'hours since\s+([\d\-]+)[ T]([\d:]+)','tokens','once');
        if isempty(tok)
            tok = regexp(un,'hours since\s+([\d\-]+)','tokens','once');
            if isempty(tok), continue; end
            t0 = datetime(tok{1},'InputFormat','yyyy-MM-dd');
        else
            t0 = datetime([tok{1} ' ' tok{2}],'InputFormat','yyyy-MM-dd HH:mm:ss');
        end
        tt = t0 + hours(double(tv(:)));
        n  = numel(tt);
        times = [times; tt];                 %#ok<AGROW>
        files = [files; repmat({fn}, n, 1)]; %#ok<AGROW>
        tidx  = [tidx;  (1:n)'];             %#ok<AGROW>
    end

    if isempty(files), return; end
    [times, o] = sort(times);
    files = files(o);
    tidx  = tidx(o);
    fidx  = (1:numel(files))';
end

function y = local_unclamp(x, floor_val, cap)
    y = x;
    y(y >= cap - 1e-6)       = NaN;
    y(y <= floor_val + 1e-6) = NaN;
end
