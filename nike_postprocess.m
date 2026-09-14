% Post-processing: reads the *.mat files and writes the per-cyclone
% summary statistics used in Tables 2 and 3. Run after the main analysis.

clc; clear; close all; warning('off','all');


%% RUN ALL STORMS
% Processes every storm in one run. To do a single storm instead,
% set STORM_LIST = {'AMPHAN'};
STORM_LIST = {'KYARR','AMPHAN','FANI','TAUKTAE'};

for iSTORM = 1:numel(STORM_LIST)
    STORM_ID = STORM_LIST{iSTORM};

    % clear per-storm workspace so nothing leaks between storms
    clearvars -except STORM_LIST iSTORM STORM_ID

    cfg = get_post_config(STORM_ID);

    fprintf('\n============================================================\n');
    fprintf(' UNIFIED CYCLONE POST-PROCESSING v7 — STORM: %s\n', cfg.name);
    fprintf('============================================================\n\n');

    if ~isfile(cfg.mat_file)
        error('Main output not found: %s\n  Run Main_combined_v7_FINAL.m first.', cfg.mat_file);
    end
    S  = load(cfg.mat_file);
    z  = S.z(:);
    nz = numel(z);
    nt = numel(S.time_tc);
    rho0 = 1025;

    %% DERIVED STORAGE
    Frac_NIKE_below_MLD    = nan(nt,1);
    NIW_Penetration_Frac   = nan(nt,1);   % deep coherent / total broadband
    Centroid_MLD_diff_m    = nan(nt,1);
    Penetration_Efficiency = nan(nt,1);
    NIW_Penetration_right  = nan(nt,1);
    NIW_Penetration_left   = nan(nt,1);
    Asymmetry_Class        = strings(nt,1);
    isStrongLC             = false(nt,1);
    isWeakLC               = false(nt,1);
    WindPower_SpinUpTime_h = nan(nt,1);
    WindPower_neg_flag     = false(nt,1);
    Horiz_wavenumber_perkm = nan(nt,1);
    BLT_flag_barrier       = false(nt,1);

    %% PER-TIMESTEP LOOP
    for it = 1:nt
        MLD_it = S.MLD_m(it);
        if isnan(MLD_it), continue; end

        Ecoh = S.NIKE_coherent_Jm2(it);
        Ebrd = S.NIKE_broadband_Jm2(it);
        if isnan(Ecoh) || Ecoh <= 0, continue; end

        % ---- penetration ----
        if isfinite(S.NIKE_below_MLD_Jm2(it))
            Frac_NIKE_below_MLD(it) = S.NIKE_below_MLD_Jm2(it) / (Ecoh + eps);
        end
        if isfinite(Ebrd) && Ebrd > 0 && isfinite(S.NIKE_below_MLD_Jm2(it))
            NIW_Penetration_Frac(it) = S.NIKE_below_MLD_Jm2(it) / (Ebrd + eps);
        end

        Centroid_MLD_diff_m(it)    = S.NIKE_centroid_m(it) - MLD_it;
        Penetration_Efficiency(it) = max(-2, min(2, Centroid_MLD_diff_m(it)/(MLD_it+eps)));

        % ---- right/left deep penetration (from saved deep integrals) ----
        if isfinite(S.NIKE_deep_right_Jm2(it)) && isfinite(S.NIKE_coherent_right_Jm2(it)) ...
                && S.NIKE_coherent_right_Jm2(it) > 0
            NIW_Penetration_right(it) = S.NIKE_deep_right_Jm2(it) / S.NIKE_coherent_right_Jm2(it);
        end
        if isfinite(S.NIKE_deep_left_Jm2(it)) && isfinite(S.NIKE_coherent_left_Jm2(it)) ...
                && S.NIKE_coherent_left_Jm2(it) > 0
            NIW_Penetration_left(it) = S.NIKE_deep_left_Jm2(it) / S.NIKE_coherent_left_Jm2(it);
        end

        % ---- Langmuir ----
        if isfinite(S.La_t(it))
            isStrongLC(it) = S.La_t(it) < 0.35;
            isWeakLC(it)   = S.La_t(it) > 0.70;
        end

        % ---- spin-up time = coherent NIKE / positive NI wind power ----
        wp = S.WindPower_NI_Wm2(it);
        if isfinite(wp)
            if wp < 0, WindPower_neg_flag(it) = true; end
            if wp > 0
                WindPower_SpinUpTime_h(it) = Ecoh / wp / 3600;
            end
        end

        % ---- barrier layer flag ----
        if isfinite(S.BLT_m(it)) && S.BLT_m(it) > 5
            BLT_flag_barrier(it) = true;
        end

        % ---- asymmetry classification ----
        ar = S.Asymmetry_Ratio(it);
        if isfinite(ar)
            if     ar >= 1.50,                 Asymmetry_Class(it) = "right-bias";
            elseif ar >= 1.18 && ar < 1.50,    Asymmetry_Class(it) = "transitional";
            elseif ar >= 0.85 && ar < 1.18,    Asymmetry_Class(it) = "symmetric";
            elseif ar >  0.67 && ar < 0.85,    Asymmetry_Class(it) = "transitional";
            else,                              Asymmetry_Class(it) = "left-bias";
            end
        else
            Asymmetry_Class(it) = "undefined";
        end
    end

    %% horizontal wavenumber from along-track coherent NIKE
    % k_h ~ 2*pi / (2 * distance between successive NIKE maxima along track)
    % Estimated once from the along-track structure, assigned per segment.
    lat = S.lat_tc; lon = S.lon_tc;
    seg_km = zeros(nt,1);
    for it = 2:nt
        seg_km(it) = distdim(distance(lat(it-1),lon(it-1),lat(it),lon(it)),'deg','km');
    end
    cumdist = cumsum(seg_km);
    nk = S.NIKE_coherent_Jm2;
    good = isfinite(nk) & isfinite(cumdist);
    % findpeaks requires a strictly increasing abscissa, which cumulative
    % along-track distance is not when a TC recurves, stalls or loops, so the
    % maxima are located by direct neighbour comparison on the retained points.
    if sum(good) > 8
        xg = cumdist(good);  yg = nk(good);
        [xg, iu] = unique(xg, 'stable');  yg = yg(iu);   % drop repeated positions
        [xg, is] = sort(xg);              yg = yg(is);   % enforce increasing x
        locs = [];
        for jj = 2:numel(yg)-1
            if yg(jj) > yg(jj-1) && yg(jj) >= yg(jj+1)
                locs(end+1,1) = xg(jj); %#ok<AGROW>
            end
        end
        if numel(locs) >= 2
            d = diff(locs);
            d = d(d > 0);
            if ~isempty(d)
                lam_h = 2*mean(d);                 % km, crude
                Horiz_wavenumber_perkm(:) = 2*pi/max(lam_h,eps);
            end
        end
    end

    %% TRACK-LEVEL STATISTICS
    st = struct();
    grab = @(x) x(isfinite(x));
    % Final finite value of a cumulative series. Cumulative wind work must be
    % read at the end of the record, not at its maximum: a curve that ends below
    % its own peak (net negative transfer) is otherwise reported as the peak.
    lastfin = @(x) x(find(isfinite(x),1,'last'));
    st.peak_NIKE_coh_kJ   = max(S.NIKE_coherent_Jm2)/1000;
    st.peak_NIKE_brd_kJ   = max(S.NIKE_broadband_Jm2)/1000;
    st.mean_NIKE_coh_kJ   = mean(grab(S.NIKE_coherent_Jm2))/1000;
    st.baseline_kJ        = S.NIKE_base_col/1000;
    st.mean_coh_ratio     = mean(grab(S.Coherence_Ratio));
    st.mean_MLD           = mean(grab(S.MLD_m));
    % [P4] Mixed-layer deepening (R1 comment 5). Pre-storm value is the median
    % of the first five valid timesteps; the maximum is taken over the record.
    % Reported in Table 2, so it must come from the committed code rather than
    % from an ad-hoc script -- the Open Research statement claims as much.
    mld_v                 = S.MLD_m(isfinite(S.MLD_m));
    st.MLD_prestorm       = median(mld_v(1:min(5,numel(mld_v))));
    st.MLD_max            = max(mld_v);
    [~, i_pk]             = max(S.NIKE_coherent_Jm2);
    st.MLD_at_NIKEpeak    = S.MLD_m(i_pk);
    st.MLD_deepening_m    = st.MLD_max - st.MLD_prestorm;
    st.MLD_deepening_pct  = 100*st.MLD_deepening_m/max(st.MLD_prestorm,eps);
    % [P5] Near-inertial share of the windowed kinetic energy of Eq. (6),
    % measured window by window on the mask-averaged surface velocity. This is
    % a LOWER BOUND: the disk mean cancels near-inertial phase, whereas Eq. (6)
    % is evaluated pointwise before averaging.
    st.NI_share_median    = NaN;
    if isfield(S,'u_surf_save') && isfield(S,'f_tc')
        i0 = ceil(size(S.u_surf_save,2)/2);
        us = S.u_surf_save(:,i0); vs = S.v_surf_save(:,i0);
        gd = isfinite(us) & isfinite(vs);
        % Need enough rows for one centred 96-h window (2*hw+1) plus a margin,
        % and most of them finite. A flat count of >40 wrongly excluded Amphan
        % (45 steps) and Tauktae (40 steps), which are the two short records.
        if numel(us) >= 2*16+3 && sum(gd) >= 0.6*numel(us)
            us = fillmissing(us,'linear','EndValues','nearest');
            vs = fillmissing(vs,'linear','EndValues','nearest');
            f_cpd = abs(mean(S.f_tc,'omitnan'))*86400/(2*pi);
            Fs = 8; fn = Fs/2;
            [bb,aa] = butter(2,[0.5*f_cpd/fn min(1.5*f_cpd/fn,0.95)],'bandpass');
            ub = filtfilt(bb,aa,detrend(us)); vb = filtfilt(bb,aa,detrend(vs));
            hw = 16; r = nan(numel(us),1);          % +/-48 h at 3 h = 33 points
            for iw = 1+hw : numel(us)-hw
                idx = iw-hw : iw+hw;
                if sum(gd(idx)) < 0.75*numel(idx), continue; end
                uw = detrend(us(idx)); vw = detrend(vs(idx));   % 96-h detrend, as Eq. (6)
                r(iw) = mean(ub(idx).^2+vb(idx).^2)/mean(uw.^2+vw.^2);
            end
            st.NI_share_median = median(r,'omitnan');
            st.NI_share_n      = sum(isfinite(r));
        end
    end
    st.mean_ILD           = mean(grab(S.ILD_m));
    st.mean_BLT           = mean(grab(S.BLT_m));
    % [P2] Barrier-layer fraction on the valid-BLT count, not track length.
    % BLT_flag_barrier can only be true where BLT is finite, so dividing by nt
    % mixes a valid-only numerator with a full-record denominator.
    st.n_valid_BLT        = sum(isfinite(S.BLT_m));
    st.barrier_frac       = 100*sum(BLT_flag_barrier)/max(st.n_valid_BLT,1);
    st.below_ML_frac      = 100*mean(grab(Frac_NIKE_below_MLD));
    st.mean_penetration   = mean(grab(NIW_Penetration_Frac));
    st.median_penetration = median(grab(NIW_Penetration_Frac));
    st.pen_right          = mean(grab(NIW_Penetration_right));
    st.pen_left           = mean(grab(NIW_Penetration_left));
    st.mean_windstress    = mean(grab(S.WindStress_Nm2));
    st.max_windstress     = max(grab(S.WindStress_Nm2));
    st.WW_signed_kJ       = lastfin(S.WindWork_signed_Jm2)/1000;   % net, end of record
    st.WW_pos_kJ          = lastfin(S.WindWork_pos_Jm2)/1000;      % gross, end of record
    st.mean_tau_uI        = mean(grab(S.WindPower_NI_Wm2));
    st.neg_WP_steps       = sum(WindPower_neg_flag);
    st.median_spinup_h    = median(grab(WindPower_SpinUpTime_h));
    st.median_phase_deg   = median(grab(S.phase_tau_uI_deg));
    st.frac_phase45       = 100*sum(abs(grab(S.phase_tau_uI_deg))<45)/numel(grab(S.phase_tau_uI_deg));
    st.mean_zeta_over_f   = mean(grab(abs(S.zeta_save)./abs(S.f_tc)));
    st.mean_strain        = mean(grab(S.strain_rate));
    st.mean_feff_over_f   = mean(grab(S.f_eff./S.f_tc));
    st.n_track            = nt;
    st.n_valid_NIKE       = numel(grab(S.NIKE_coherent_Jm2));
    st.n_valid_WP         = numel(grab(S.WindPower_NI_Wm2));
    st.n_valid_phase      = numel(grab(S.phase_tau_uI_deg));
    st.n_valid_MLD        = numel(grab(S.MLD_m));
    st.n_valid_Ri         = numel(grab(S.Ri_min_100m));
    st.n_valid_AR         = numel(grab(S.Asymmetry_Ratio));
    st.mean_Rmin100       = median(grab(S.Ri_min_100m));
    st.frac_Ri_subcrit    = 100*sum(grab(S.Ri_min_100m)<0.25)/numel(grab(S.Ri_min_100m));
    st.mean_Lat           = mean(grab(S.La_t));
    st.n_valid_Lat        = numel(grab(S.La_t));
    st.strongLC_frac      = 100*sum(isStrongLC)/max(st.n_valid_Lat,1);
    st.mean_Rratio        = mean(grab(S.Asymmetry_Ratio));
    st.median_Rratio      = median(grab(S.Asymmetry_Ratio));
    % [P3] Asymmetry class fractions on the count of timesteps where AR is
    % DEFINED, not on track length. AR is NaN wherever the bilateral 10% guard
    % fails; an undefined timestep is not evidence of symmetry and must not
    % dilute the right-bias fraction. 'transitional' is a defined class and is
    % therefore included, so the four fractions sum to 100%.
    % NOTE: count from the numeric array, NOT from Asymmetry_Class. The class
    % array is initialised as strings(nt,1), so any timestep skipped by a
    % 'continue' in the main loop stays as an empty string "" rather than
    % "undefined", and a  ~= "undefined"  test would wrongly count it.
    st.n_defined_AR       = sum(isfinite(S.Asymmetry_Ratio));
    st.rightbias_frac     = 100*sum(Asymmetry_Class=="right-bias")/max(st.n_defined_AR,1);
    st.symmetric_frac     = 100*sum(Asymmetry_Class=="symmetric")/max(st.n_defined_AR,1);
    st.leftbias_frac      = 100*sum(Asymmetry_Class=="left-bias")/max(st.n_defined_AR,1);
    st.transitional_frac  = 100*sum(Asymmetry_Class=="transitional")/max(st.n_defined_AR,1);
    st.phase_sd_deg       = mean(grab(S.phase_sd_deg));
    st.mean_SST           = mean(grab(S.SST_C));
    st.mean_OHC_kJcm2     = mean(grab(S.OHC_Jm2))/1e7;

    %% BUILD POST TABLE
    Tpost = table( ...
        S.time_tc, S.lat_tc, S.lon_tc, S.heading_deg_save, ...
        S.MLD_m, S.ILD_m, S.BLT_m, BLT_flag_barrier, ...
        S.NIKE_coherent_Jm2, S.NIKE_broadband_Jm2, ...
        S.NIKE_above_MLD_Jm2, S.NIKE_below_MLD_Jm2, ...
        Frac_NIKE_below_MLD, NIW_Penetration_Frac, ...
        S.NIKE_centroid_m, Centroid_MLD_diff_m, Penetration_Efficiency, ...
        S.NIKE_coherent_right_Jm2, S.NIKE_coherent_left_Jm2, ...
        S.NIKE_broadband_right_Jm2, S.NIKE_broadband_left_Jm2, ...
        S.NIKE_deep_right_Jm2, S.NIKE_deep_left_Jm2, ...
        NIW_Penetration_right, NIW_Penetration_left, ...
        S.Asymmetry_Ratio, Asymmetry_Class, ...
        S.Coherence_Ratio, ...
        S.WindStress_Nm2, S.WindStress_right_Nm2, S.WindStress_left_Nm2, ...
        S.WindPower_Wm2, S.WindPower_NI_Wm2, ...
        S.WindPower_right_Wm2, S.WindPower_left_Wm2, ...
        S.WindWork_signed_Jm2, S.WindWork_pos_Jm2, ...
        WindPower_SpinUpTime_h, WindPower_neg_flag, ...
        S.phase_tau_uI_deg, S.uI_amp_surf, S.phase_sd_deg, ...
        S.f_tc, S.f_eff, S.zeta_save, S.strain_rate, Horiz_wavenumber_perkm, ...
        S.Ri_Bulk_MLD, S.Ri_min_100m, S.Ri_med_100m, ...
        S.La_t, isStrongLC, isWeakLC, ...
        S.SST_C, S.OHC_Jm2, S.Stokes_ms, S.ML_Buoyancy, ...
        S.Inertial_hr, S.n_cells_used);

    Tpost.Properties.VariableNames = { ...
        'time','lat','lon','heading_deg', ...
        'MLD_m','ILD_m','BLT_m','has_barrier_layer', ...
        'NIKE_coherent_Jm2','NIKE_broadband_Jm2', ...
        'NIKE_above_MLD_Jm2','NIKE_below_MLD_Jm2', ...
        'Frac_below_MLD_coh','NIW_Penetration_Frac', ...
        'NIKE_centroid_m','Centroid_MLD_diff_m','Penetration_Efficiency', ...
        'NIKE_coh_right_Jm2','NIKE_coh_left_Jm2', ...
        'NIKE_brd_right_Jm2','NIKE_brd_left_Jm2', ...
        'NIKE_deep_right_Jm2','NIKE_deep_left_Jm2', ...
        'Penetration_right','Penetration_left', ...
        'Asymmetry_Ratio','Asymmetry_Class', ...
        'Coherence_Ratio', ...
        'WindStress_Nm2','WindStress_right_Nm2','WindStress_left_Nm2', ...
        'WindPower_total_Wm2','WindPower_NI_Wm2', ...
        'WindPower_right_Wm2','WindPower_left_Wm2', ...
        'WindWork_signed_Jm2','WindWork_pos_Jm2', ...
        'SpinUpTime_h','WindPower_negative', ...
        'phase_tau_uI_deg','uI_amp_surf_ms','phase_sd_deg', ...
        'f_tc','f_eff','zeta','strain_rate','k_horiz_perkm', ...
        'Ri_Bulk_MLD','Ri_min_100m','Ri_med_100m', ...
        'La_t','isStrongLC','isWeakLC', ...
        'SST_C','OHC_Jm2','Stokes_ms','ML_Buoyancy', ...
        'Inertial_hr','n_cells_used'};

    %% SAVE
    writetable(Tpost, cfg.out_xlsx);
    save(cfg.out_mat, 'Tpost', 'st', 'z', '-v7.3');

    fprintf('\nsaved: %s\n', cfg.out_xlsx);
    fprintf('saved: %s\n', cfg.out_mat);
    fprintf('\n=== %s POST-PROCESSING v7 : COMPLETE ===\n', cfg.name);


end   % STORM_LIST loop

fprintf('\n============================================================\n');
fprintf(' ALL STORMS COMPLETE\n');
fprintf('============================================================\n');

%% CONFIG
function cfg = get_post_config(storm_id)
    switch upper(storm_id)
        case 'KYARR'
            cfg.name='KYARR';   cfg.mat_file='KYARR_V17.mat';
        case 'AMPHAN'
            cfg.name='AMPHAN';  cfg.mat_file='AMPHAN_V17.mat';
        case 'FANI'
            cfg.name='FANI';    cfg.mat_file='FANI_V17.mat';
        case 'TAUKTAE'
            cfg.name='TAUKTAE'; cfg.mat_file='TAUKTAE_V17.mat';
        otherwise
            error('Unknown STORM_ID: %s', storm_id);
    end
    cfg.out_xlsx    = sprintf('%s_POST_V17.xlsx', cfg.name);
    cfg.out_mat     = sprintf('%s_POST_V17.mat',  cfg.name);
end