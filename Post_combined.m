clc; clear; close all;

% Unified cyclone post-processing. Reads the core output produced by
% Main_combined.m and computes penetration, conversion efficiency,
% Ri statistics, R/L asymmetry classification, and Langmuir flags.
% Output: <STORM>_POST_CONSOLIDATED_V16.{xlsx,mat} + *_POST_SUMMARY_v6.txt

STORM_ID = 'KYARR';   % KYARR | AMPHAN | FANI | TAUKTAE

cfg = get_post_config(STORM_ID);

fprintf('\n=== %s POST-PROCESSING ===\n', cfg.name);

mat_file  = cfg.mat_file;
xlsx_file = cfg.xlsx_file;
if ~isfile(mat_file) || ~isfile(xlsx_file)
    error('Main core output not found. Run Main_combined.m first.\n  Expected: %s and %s', ...
        mat_file, xlsx_file);
end

load(mat_file);
Tmain = readtable(xlsx_file);
nt    = height(Tmain);

z  = z(:);
nz = length(z);

rho0  = 1025;
g     = 9.81;
Omega = 7.2921e-5;

% Cap coherence > 1.0 (harmonic-fit overshoot in low-coherence regimes)
if ismember('Coherence_Ratio', Tmain.Properties.VariableNames)
    over_one = sum(Tmain.Coherence_Ratio > 1, 'omitnan');
    if over_one > 0
        fprintf('Capping %d coherence values > 1.0\n', over_one);
        Tmain.Coherence_Ratio(Tmain.Coherence_Ratio > 1) = 1.0;
    end
end

% Storage
Frac_NIKE_below_MLD       = nan(nt,1);
NIW_Penetration_Frac      = nan(nt,1);
Centroid_MLD_diff_m       = nan(nt,1);
Penetration_Efficiency    = nan(nt,1);

isStrongLC                = false(nt,1);
isWeakLC                  = false(nt,1);

Ri_actual_min_100m        = nan(nt,1);
Ri_actual_median_100m     = nan(nt,1);
Ri_wave_proxy_min         = nan(nt,1);

WindPower_SpinUpTime_s    = nan(nt,1);
WindPower_negative_flag   = false(nt,1);
WindPower_ConvEff_dimless = nan(nt,1);   % capped at 1.0
WindPower_ConvEff_raw     = nan(nt,1);   % uncapped
WindPower_ConvEff_capped  = false(nt,1);

NIW_Penetration_right     = nan(nt,1);
NIW_Penetration_left      = nan(nt,1);
Asymmetry_Class           = strings(nt,1);

% Main loop
for it = 1:nt

    MLD_it = Tmain.MLD(it);
    if isnan(MLD_it), continue; end

    Etot_coh = Tmain.NIKE_Coherent(it);
    Etot_brd = Tmain.NIKE_Broadband(it);
    if Etot_coh <= 0 || isnan(Etot_coh), continue; end

    Frac_NIKE_below_MLD(it) = Tmain.NIKE_Deep(it) / (Etot_coh + eps);

    if isfinite(Etot_brd) && Etot_brd > 0
        NIW_Penetration_Frac(it) = Tmain.NIKE_Deep(it) / (Etot_brd + eps);
    end

    Centroid_MLD_diff_m(it)    = Tmain.Centroid(it) - MLD_it;
    Penetration_Efficiency(it) = max(-2, min(2, Centroid_MLD_diff_m(it) / (MLD_it + eps)));

    if ismember('La_t', Tmain.Properties.VariableNames)
        if ~isnan(Tmain.La_t(it))
            isStrongLC(it) = Tmain.La_t(it) < 0.35;
            isWeakLC(it)   = Tmain.La_t(it) > 0.7;
        end
    end

    % Ri: upper-100m scan, 5th percentile as robust min (Jampana 2018);
    % also a wave-band proxy in the MLD..250m range
    if it <= size(Ri_save,1)
        Ri_prof = Ri_save(it,:)';
        z_Ri    = z(1:length(Ri_prof));

        idx_100 = z_Ri <= 100;
        if any(idx_100)
            Ri_upper = Ri_prof(idx_100);
            Ri_upper = Ri_upper(isfinite(Ri_upper));
            if numel(Ri_upper) >= 5
                Ri_actual_min_100m(it)    = prctile(Ri_upper, 5);
                Ri_actual_median_100m(it) = median(Ri_upper);
            elseif ~isempty(Ri_upper)
                Ri_actual_min_100m(it)    = min(Ri_upper);
                Ri_actual_median_100m(it) = median(Ri_upper);
            end
        end

        idx_wave = z_Ri > MLD_it & z_Ri < 250;
        if any(idx_wave)
            Ri_wave = Ri_prof(idx_wave);
            Ri_wave = Ri_wave(isfinite(Ri_wave));
            if numel(Ri_wave) >= 5
                Ri_wave_proxy_min(it) = prctile(Ri_wave, 5);
            elseif ~isempty(Ri_wave)
                Ri_wave_proxy_min(it) = min(Ri_wave);
            end
        end
    end

    % Wind power spin-up time
    if ismember('WindPower_Wm2', Tmain.Properties.VariableNames)
        wp = Tmain.WindPower_Wm2(it);
        if isfinite(wp)
            if wp < 0, WindPower_negative_flag(it) = true; end
            if wp > 0 && isfinite(Etot_coh)
                WindPower_SpinUpTime_s(it) = Etot_coh / wp;
            end
        end
    end

    % Conversion efficiency: capped value for the headline figure,
    % raw value for the wind-deficit diagnostic
    if ismember('WindWork_cum_Jm2', Tmain.Properties.VariableNames)
        ww = Tmain.WindWork_cum_Jm2(it);
        if isfinite(ww) && ww > 0 && isfinite(Etot_coh)
            ce_raw = Etot_coh / ww;
            WindPower_ConvEff_dimless(it) = min(ce_raw, 1.0);
            WindPower_ConvEff_raw(it)     = ce_raw;
            if ce_raw > 1.0
                WindPower_ConvEff_capped(it) = true;
            end
        end
    end

    % Right/left penetration fractions (depth-resolved)
    if exist('NIKE_coherent_z_right','var') && exist('NIKE_coherent_z_left','var')
        prof_R = NIKE_coherent_z_right(it,:)';
        prof_L = NIKE_coherent_z_left(it,:)';
        idx_deep = (z > MLD_it);

        v_R_deep = isfinite(prof_R) & idx_deep;
        v_R_all  = isfinite(prof_R);
        if sum(v_R_deep) >= 2 && sum(v_R_all) >= 2
            E_R_deep  = trapz(z(v_R_deep), prof_R(v_R_deep));
            E_R_total = trapz(z(v_R_all),  prof_R(v_R_all));
            if E_R_total > 0
                NIW_Penetration_right(it) = E_R_deep / E_R_total;
            end
        end
        v_L_deep = isfinite(prof_L) & idx_deep;
        v_L_all  = isfinite(prof_L);
        if sum(v_L_deep) >= 2 && sum(v_L_all) >= 2
            E_L_deep  = trapz(z(v_L_deep), prof_L(v_L_deep));
            E_L_total = trapz(z(v_L_all),  prof_L(v_L_all));
            if E_L_total > 0
                NIW_Penetration_left(it) = E_L_deep / E_L_total;
            end
        end
    end

    % Asymmetry class (four bins around symmetric +/-15%)
    if ismember('Asymmetry_Ratio', Tmain.Properties.VariableNames)
        ar = Tmain.Asymmetry_Ratio(it);
        if isfinite(ar)
            if     ar >= 1.5,                Asymmetry_Class(it) = "right-bias";
            elseif ar >= 1.18 && ar < 1.5,   Asymmetry_Class(it) = "transitional";
            elseif ar >= 0.85 && ar < 1.18,  Asymmetry_Class(it) = "symmetric";
            elseif ar > 0.67  && ar < 0.85,  Asymmetry_Class(it) = "transitional";
            else,                             Asymmetry_Class(it) = "left-bias";
            end
        else
            Asymmetry_Class(it) = "undefined";
        end
    end

end

% Pull through optional columns
get_col = @(name) ifelse(ismember(name, Tmain.Properties.VariableNames), ...
    Tmain.(char(name)), nan(nt,1));

La_t_col      = get_col('La_t');
WindPower_col = get_col('WindPower_Wm2');
WindWork_col  = get_col('WindWork_cum_Jm2');
Stress_col    = get_col('Stress');
dPE_dt_col    = get_col('dPE_dt');

NIKE_Coh_R   = get_col('NIKE_Coherent_right');
NIKE_Coh_L   = get_col('NIKE_Coherent_left');
NIKE_Brd_R   = get_col('NIKE_Broadband_right');
NIKE_Brd_L   = get_col('NIKE_Broadband_left');
Asym_Ratio   = get_col('Asymmetry_Ratio');
Stress_R     = get_col('WindStress_right');
Stress_L     = get_col('WindStress_left');
WP_R         = get_col('WindPower_right');
WP_L         = get_col('WindPower_left');
Heading_deg  = get_col('TC_heading_deg');
R_NI_km      = get_col('R_NI_used_km');

% Final table
Tpost = table( ...
    Tmain.time, Tmain.lat, Tmain.lon, Tmain.MLD, ...
    Tmain.NIKE_Coherent, Tmain.NIKE_Broadband, ...
    Tmain.NIKE_ML, Tmain.NIKE_Deep, ...
    Frac_NIKE_below_MLD, NIW_Penetration_Frac, ...
    Tmain.Centroid, Centroid_MLD_diff_m, Penetration_Efficiency, ...
    La_t_col, isStrongLC, isWeakLC, Tmain.Coherence_Ratio, ...
    Tmain.Ri_Bulk, Ri_actual_min_100m, Ri_actual_median_100m, Ri_wave_proxy_min, ...
    Tmain.NIKE_WKB, Tmain.NIW_PE, Tmain.PE_KE_Ratio, ...
    Tmain.Cgz_Proxy, Tmain.EnergyLoss, ...
    WindPower_col, WindWork_col, ...
    WindPower_SpinUpTime_s, WindPower_ConvEff_dimless, WindPower_ConvEff_raw, ...
    WindPower_ConvEff_capped, WindPower_negative_flag, ...
    Stress_col, dPE_dt_col, ...
    NIKE_Coh_R, NIKE_Coh_L, NIKE_Brd_R, NIKE_Brd_L, ...
    Asym_Ratio, Asymmetry_Class, ...
    NIW_Penetration_right, NIW_Penetration_left, ...
    Stress_R, Stress_L, WP_R, WP_L, ...
    Heading_deg, R_NI_km);

Tpost.Properties.VariableNames = { ...
    'time','lat','lon','MLD_m', ...
    'NIKE_Coherent_Jm2','NIKE_Broadband_Jm2', ...
    'NIKE_above_MLD_Jm2','NIKE_below_MLD_Jm2', ...
    'Frac_NIKE_below_MLD','NIW_Penetration_Frac', ...
    'NIKE_Centroid_m','Centroid_MLD_diff_m','Penetration_Efficiency', ...
    'La_t','isStrongLC','isWeakLC','Coherence_Ratio', ...
    'Ri_Bulk_MLD','Ri_actual_min_100m','Ri_actual_median_100m','Ri_wave_proxy_min', ...
    'NIKE_WKB_Jm2','NIW_PE_Jm2','PE_KE_Ratio', ...
    'Inferred_Cgz_m_day','Eff_Energy_Loss_Wm2', ...
    'WindPower_Wm2','WindWork_cum_Jm2', ...
    'WindPower_SpinUpTime_s','WindPower_ConvEff_dimless','WindPower_ConvEff_raw', ...
    'WindPower_ConvEff_capped','WindPower_negative_flag', ...
    'WindStress_Nm2','dPE_dt_Wm2', ...
    'NIKE_Coherent_right_Jm2','NIKE_Coherent_left_Jm2', ...
    'NIKE_Broadband_right_Jm2','NIKE_Broadband_left_Jm2', ...
    'Asymmetry_Ratio','Asymmetry_Class', ...
    'NIW_Penetration_right','NIW_Penetration_left', ...
    'WindStress_right_Nm2','WindStress_left_Nm2', ...
    'WindPower_right_Wm2','WindPower_left_Wm2', ...
    'TC_heading_deg','R_NI_used_km'};

writetable(Tpost, cfg.out_xlsx);
save(cfg.out_mat, 'Tpost', 'z');

% Summary report
storm_tag = upper(extractBefore(xlsx_file, '_V16'));
summary_file = sprintf('%s_POST_SUMMARY_v6.txt', storm_tag);
fid = fopen(summary_file, 'w');
wl = @(varargin) [fprintf(varargin{:}), fprintf(fid, varargin{:})];

wl('\n=== %s POST SUMMARY ===\n', storm_tag);
wl('Total timesteps          : %d\n', nt);

wl('\n--- PENETRATION (broadband denominator) ---\n');
mean_MLD = mean(Tmain.MLD,'omitnan');
NPF = NIW_Penetration_Frac(isfinite(NIW_Penetration_Frac));
if ~isempty(NPF)
    wl('Mean MLD                 : %.1f m\n', mean_MLD);
    wl('Mean NIW penetration     : %.3f\n', mean(NPF));
    wl('Median NIW penetration   : %.3f\n', median(NPF));
    wl('Max NIW penetration      : %.3f\n', max(NPF));
    wl('Steps in 0.10-0.45       : %.1f %%\n', 100*sum(NPF>=0.10 & NPF<=0.45)/length(NPF));
end

wl('\n--- WIND CONVERSION ---\n');
SUT = WindPower_SpinUpTime_s(isfinite(WindPower_SpinUpTime_s));
CE  = WindPower_ConvEff_dimless(isfinite(WindPower_ConvEff_dimless));
n_capped = sum(WindPower_ConvEff_capped);

if ~isempty(SUT)
    wl('Median spin-up time      : %.0f s (= %.1f h)\n', median(SUT), median(SUT)/3600);
    wl('Mean spin-up time        : %.0f s (= %.1f h)\n', mean(SUT),   mean(SUT)/3600);
end
if ~isempty(CE)
    wl('Median conversion eff.   : %.3f\n', median(CE));
    wl('Mean conversion eff.     : %.3f\n', mean(CE));
    wl('Steps in 0.20-0.50       : %.1f %%\n', 100*sum(CE>=0.20 & CE<=0.50)/length(CE));
    wl('Steps at cap (=1.0)      : %d / %d\n', n_capped, length(CE));
end
wl('Negative-WP timesteps    : %d\n', sum(WindPower_negative_flag));

wl('\n--- MIXING (Ri, upper 100 m) ---\n');
Ri_min = Ri_actual_min_100m(isfinite(Ri_actual_min_100m));
Ri_med = Ri_actual_median_100m(isfinite(Ri_actual_median_100m));
if ~isempty(Ri_min)
    wl('Min Ri (upper 100m)      : %.3f\n', min(Ri_min));
    wl('Median of per-step 5%%ile : %.3f\n', median(Ri_min));
    wl('Steps with Ri < 0.25     : %d / %d\n', sum(Ri_min < 0.25), length(Ri_min));
end
if ~isempty(Ri_med)
    wl('Mean of per-step median  : %.3f\n', mean(Ri_med));
end

wl('\n--- ASYMMETRY ---\n');
AR = Asym_Ratio(isfinite(Asym_Ratio));
if ~isempty(AR)
    wl('Mean R/L ratio           : %.2f\n', mean(AR));
    wl('Median R/L ratio         : %.2f\n', median(AR));
    wl('Right-bias fraction      : %.1f %%\n', 100*sum(Asymmetry_Class == "right-bias")/nt);
    wl('Symmetric fraction       : %.1f %%\n', 100*sum(Asymmetry_Class == "symmetric") /nt);
    wl('Left-bias fraction       : %.1f %%\n', 100*sum(Asymmetry_Class == "left-bias") /nt);
    wl('Transitional fraction    : %.1f %%\n', 100*sum(Asymmetry_Class == "transitional")/nt);
end
NPR = NIW_Penetration_right(isfinite(NIW_Penetration_right));
NPL = NIW_Penetration_left(isfinite(NIW_Penetration_left));
if ~isempty(NPR), wl('Mean penetration (right) : %.3f\n', mean(NPR)); end
if ~isempty(NPL), wl('Mean penetration (left)  : %.3f\n', mean(NPL)); end

wl('\n--- LANGMUIR ---\n');
wl('Strong LC fraction       : %.1f %% (La_t < 0.35)\n', 100*sum(isStrongLC)/nt);
wl('Weak LC fraction         : %.1f %% (La_t > 0.70)\n', 100*sum(isWeakLC)/nt);

fclose(fid);
fprintf('Summary saved: %s\n', summary_file);
fprintf('=== %s POST: DONE ===\n', cfg.name);


% --- Local helpers ---
function out = ifelse(cond, a, b)
    if cond, out = a; else, out = b; end
end

function cfg = get_post_config(storm_id)
% Per-storm output filenames; should match those produced by Main_combined.m.
    switch upper(storm_id)
        case 'KYARR'
            cfg.name      = 'KYARR';
            cfg.mat_file  = 'KYARR_V16_SUPERCHARGED.mat';
            cfg.xlsx_file = 'KYARR_V16_SUPERCHARGED.xlsx';
            cfg.out_xlsx  = 'KYARR_POST_CONSOLIDATED_V16.xlsx';
            cfg.out_mat   = 'KYARR_POST_CONSOLIDATED_V16.mat';

        case 'AMPHAN'
            cfg.name      = 'AMPHAN';
            cfg.mat_file  = 'AMPHAN_V16_SUPERCHARGED.mat';
            cfg.xlsx_file = 'AMPHAN_V16_SUPERCHARGED.xlsx';
            cfg.out_xlsx  = 'AMPHAN_POST_CONSOLIDATED_V16.xlsx';
            cfg.out_mat   = 'AMPHAN_POST_CONSOLIDATED_V16.mat';

        case 'FANI'
            cfg.name      = 'FANI';
            cfg.mat_file  = 'FANI_V16_SUPERCHARGED.mat';
            cfg.xlsx_file = 'FANI_V16_SUPERCHARGED.xlsx';
            cfg.out_xlsx  = 'FANI_POST_CONSOLIDATED_V16.xlsx';
            cfg.out_mat   = 'FANI_POST_CONSOLIDATED_V16.mat';

        case 'TAUKTAE'
            cfg.name      = 'TAUKTAE';
            cfg.mat_file  = 'TAUKTAE_V16_SUPERCHARGED.mat';
            cfg.xlsx_file = 'TAUKTAE_V16_SUPERCHARGED.xlsx';
            cfg.out_xlsx  = 'TAUKTAE_POST_CONSOLIDATED_V16.xlsx';
            cfg.out_mat   = 'TAUKTAE_POST_CONSOLIDATED_V16.mat';

        otherwise
            error('Unknown STORM_ID: %s. Valid: KYARR | AMPHAN | FANI | TAUKTAE', storm_id);
    end
end
