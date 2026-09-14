% Main analysis: reads HYCOM and ERA5, fits the near-inertial oscillation at
% each grid cell. Run this first.

clc; clear; close all; warning('off','all');

DATA_DIR = 'path/to/data';        % folder with HYCOM, IMD, Stress, ARGO etc.
GSW_PATH = 'path/to/GSW-Matlab';  % Gibbs SeaWater toolbox


%% USER SETTINGS
%% RUN ALL STORMS
% Processes every storm in one run. For a single storm set
% STORM_LIST = {'AMPHAN'};
STORM_LIST    = {'KYARR','AMPHAN','FANI','TAUKTAE'};
WIN_HR        = 48;
BASELINE_DAYS = 5;
DO_POINTWISE  = true;

for iSTORM = 1:numel(STORM_LIST)
    STORM_ID = STORM_LIST{iSTORM};
    clearvars -except STORM_LIST iSTORM STORM_ID WIN_HR BASELINE_DAYS DO_POINTWISE


    cfg = get_storm_config(STORM_ID);

    fprintf('\n============================================================\n');
    fprintf(' UNIFIED CYCLONE NIKE CORE v7 FINAL — STORM: %s\n', cfg.name);
    fprintf(' window = +/-%d h | pointwise = %d | baseline = %d d\n', ...
        WIN_HR, DO_POINTWISE, BASELINE_DAYS);
    fprintf('============================================================\n\n');

    %% PATHS & CONSTANTS
    GSW_PATH = GSW_PATH;
    if exist(GSW_PATH,'dir'), addpath(genpath(GSW_PATH)); end

    g          = 9.81;
    Omega      = 7.2921e-5;
    rho0       = 1025;
    DELTA_RHO  = 0.125;
    DELTA_T    = 0.2;
    Hmax       = 700;
    R_NI       = 150;
    twin_hours = -WIN_HR:3:WIN_HR;
    tsec       = twin_hours(:)*3600;
    nw         = numel(twin_hours);

    %% TRACK
    Ttrk    = readtable(cfg.track_path);
    time_tc = datetime(Ttrk.year,Ttrk.month,Ttrk.day,Ttrk.hour,0,0);
    lat_tc  = Ttrk.lat;  lon_tc = Ttrk.lon;
    nt      = numel(time_tc);

    heading_rad = compute_heading(lat_tc, lon_tc, cfg.heading_smooth_window);

    %% HYCOM GRID
    uv_files = dir(fullfile(cfg.uv_path,'*.nc'));
    ts_files = dir(fullfile(cfg.ts_path,'*.nc'));
    if isempty(uv_files), error('No UV3Z files in %s', cfg.uv_path); end
    if isempty(ts_files), error('No TS3Z files in %s', cfg.ts_path); end

    fn0   = fullfile(uv_files(1).folder, uv_files(1).name);
    lon_g = ncread(fn0,'lon');  lat_g = ncread(fn0,'lat');
    depth = ncread(fn0,'depth');
    idz   = find(depth <= Hmax);
    z     = depth(idz);
    nz    = numel(z);

    fprintf('HYCOM grid : %d lon x %d lat, %d levels <= %d m\n', ...
        numel(lon_g), numel(lat_g), nz, Hmax);
    fprintf('levels     : '); fprintf('%g ', z); fprintf('\n\n');

    [uv_times, uv_fidx, uv_tidx] = index_time_axis(uv_files);
    [ts_times, ts_fidx, ts_tidx] = index_time_axis(ts_files);
    fprintf('UV3Z timestamps: %d  (%s to %s)\n', numel(uv_times), ...
        datestr(uv_times(1)), datestr(uv_times(end)));
    fprintf('TS3Z timestamps: %d  (%s to %s)\n\n', numel(ts_times), ...
        datestr(ts_times(1)), datestr(ts_times(end)));

    %% EXTERNAL DATA
    t_ref_dt  = datetime(cfg.year_ref,1,1);
    t_e       = t_ref_dt + hours(ncread(cfg.stress_file,'time'));
    tau_x_all = ncread(cfg.stress_file,'iews');
    tau_y_all = ncread(cfg.stress_file,'inss');
    try
        lon_e = ncread(cfg.stress_file,'lon');  lat_e = ncread(cfg.stress_file,'lat');
    catch
        lon_e = ncread(cfg.stress_file,'longitude'); lat_e = ncread(cfg.stress_file,'latitude');
    end
    [LonE,LatE] = meshgrid(lon_e, lat_e);

    [sst_all, LonS, LatS, t_sst] = load_aux(cfg.sst_file, 'sst', t_ref_dt, t_e);
    [ohc_all, LonO, LatO, t_ohc] = load_aux(cfg.ohc_file, {'ohc','ohc700'}, t_ref_dt, t_e);
    [ust_all, LonSt, LatSt, t_st] = load_aux(cfg.stokes_file, {'u_stokes','ust'}, t_ref_dt, t_e);
    [vst_all, ~, ~, ~]            = load_aux(cfg.stokes_file, {'v_stokes','vst'}, t_ref_dt, t_e);

    %% STORAGE (all outputs)
    NIKE_coherent_Jm2       = nan(nt,1);   NIKE_broadband_Jm2      = nan(nt,1);
    NIKE_above_MLD_Jm2      = nan(nt,1);   NIKE_below_MLD_Jm2      = nan(nt,1);
    NIKE_centroid_m         = nan(nt,1);   NIKE_surface_frac       = nan(nt,1);
    NIKE_coherent_z_save    = nan(nt,nz);
    NIKE_broadband_z_save   = nan(nt,nz);
    NIKE_coherent_z_right   = nan(nt,nz);
    NIKE_coherent_z_left    = nan(nt,nz);
    NIKE_broadband_z_right  = nan(nt,nz);
    NIKE_broadband_z_left   = nan(nt,nz);
    NIKE_coherent_right_Jm2 = nan(nt,1);   NIKE_coherent_left_Jm2  = nan(nt,1);
    NIKE_broadband_right_Jm2= nan(nt,1);   NIKE_broadband_left_Jm2 = nan(nt,1);
    NIKE_deep_right_Jm2     = nan(nt,1);   NIKE_deep_left_Jm2      = nan(nt,1);
    Asymmetry_Ratio         = nan(nt,1);
    Coherence_Ratio         = nan(nt,1);
    MLD_m                   = nan(nt,1);   ILD_m                   = nan(nt,1);
    BLT_m                   = nan(nt,1);
    rho_save                = nan(nt,nz);  N2_save                 = nan(nt,nz-1);
    Ri_save                 = nan(nt,nz-1);S2_save                 = nan(nt,nz-1);
    SST_C                   = nan(nt,1);   OHC_Jm2                 = nan(nt,1);
    Stokes_ms               = nan(nt,1);   La_t                    = nan(nt,1);
    WindStress_Nm2          = nan(nt,1);
    WindStress_right_Nm2    = nan(nt,1);   WindStress_left_Nm2     = nan(nt,1);
    WindPower_Wm2           = nan(nt,1);
    WindPower_right_Wm2     = nan(nt,1);   WindPower_left_Wm2      = nan(nt,1);
    WindPower_NI_Wm2        = nan(nt,1);   % tau . u_I (PF06 form)
    f_tc                    = nan(nt,1);   f_eff                   = nan(nt,1);
    zeta_save               = nan(nt,1);   strain_rate             = nan(nt,1);
    Inertial_hr             = nan(nt,1);   ML_Buoyancy             = nan(nt,1);
    Ri_Bulk_MLD             = nan(nt,1);
    Ri_min_100m             = nan(nt,1);   Ri_med_100m             = nan(nt,1);
    phase_tau_uI_deg        = nan(nt,1);
    uI_amp_surf             = nan(nt,1);
    heading_deg_save        = nan(nt,1);
    n_cells_used            = nan(nt,1);
    phase_sd_deg            = nan(nt,1);
    u_surf_save             = nan(nt,nw);  v_surf_save             = nan(nt,nw);
    PE_to_KE                = nan(nt,1);    NIW_PE_Jm2              = nan(nt,1);

    %% PRE-STORM BASELINE
    fprintf('--- computing pre-storm baseline ---\n');
    t_base = time_tc(1) - days(BASELINE_DAYS);
    [NIKE_base_z, NIKE_base_col] = nike_at_point( ...
        t_base, lat_tc(1), lon_tc(1), ...
        lon_g, lat_g, z, idz, R_NI, twin_hours, tsec, ...
        uv_times, uv_fidx, uv_tidx, uv_files, rho0, DO_POINTWISE);
    fprintf('baseline NIKE (column) : %.2f kJ/m2  at %s\n\n', ...
        NIKE_base_col/1000, datestr(t_base));

    %% MAIN LOOP
    for it = 1:nt

        f_raw = 2*Omega*sind(lat_tc(it));
        f_val = sign(f_raw)*max(abs(f_raw),1e-5);
        if f_raw == 0, f_val = 1e-5; end
        f_tc(it) = f_val;
        Inertial_hr(it) = 2*pi/abs(f_val)/3600;

        h_rad = heading_rad(it);
        heading_deg_save(it) = h_rad*180/pi;

        % ---- subgrid box ----
        ix = find(lon_g >= lon_tc(it)-2 & lon_g <= lon_tc(it)+2);
        iy = find(lat_g >= lat_tc(it)-2 & lat_g <= lat_tc(it)+2);
        if isempty(ix)||isempty(iy), continue; end
        [LonSub,LatSub] = meshgrid(lon_g(ix), lat_g(iy));
        dist = distdim(distance(lat_tc(it),lon_tc(it),LatSub,LonSub),'deg','km');

        mask   = dist <= R_NI;
        dy_pt  = (LatSub - lat_tc(it))*111;
        dx_pt  = (LonSub - lon_tc(it))*111*cosd(lat_tc(it));
        cross_z = dx_pt*sin(h_rad) - dy_pt*cos(h_rad);
        mask_R = mask & (cross_z >= 0);
        mask_L = mask & (cross_z <  0);

        ny = numel(iy);  nx = numel(ix);  ncell = ny*nx;

        % ---- read velocity cube over the window ----
        U = nan(nw, ncell, nz);
        V = nan(nw, ncell, nz);
        got = false(nw,1);
        for jt = 1:nw
            tgt = time_tc(it) + hours(twin_hours(jt));
            [fi, ti] = find_time(uv_times, uv_fidx, uv_tidx, tgt);
            if isempty(fi), continue; end
            fn = fullfile(uv_files(fi).folder, uv_files(fi).name);
            try
                ub = ncread(fn,'water_u',[ix(1) iy(1) idz(1) ti],[nx ny nz 1]);
                vb = ncread(fn,'water_v',[ix(1) iy(1) idz(1) ti],[nx ny nz 1]);
            catch
                continue;
            end
            for k = 1:nz
                su = ub(:,:,k).';   sv = vb(:,:,k).';
                U(jt,:,k) = su(:).';
                V(jt,:,k) = sv(:).';
            end
            got(jt) = true;
        end
        if sum(got) < 11, continue; end

        % ---- thermodynamics (TS3Z) ----
        [fi, ti] = find_time(ts_times, ts_fidx, ts_tidx, time_tc(it));
        if ~isempty(fi)
            fnT = fullfile(ts_files(fi).folder, ts_files(fi).name);
            tb = ncread(fnT,'water_temp',[ix(1) iy(1) idz(1) ti],[nx ny nz 1]);
            sb = ncread(fnT,'salinity'  ,[ix(1) iy(1) idz(1) ti],[nx ny nz 1]);
            T_p = squeeze(mean(tb,[1 2],'omitnan'));
            S_p = squeeze(mean(sb,[1 2],'omitnan'));
            if size(T_p,1)==1, T_p = T_p.'; S_p = S_p.'; end

            P   = gsw_p_from_z(-z, lat_tc(it));
            SA  = gsw_SA_from_SP(S_p,P,lon_tc(it),lat_tc(it));
            CT  = gsw_CT_from_pt(SA,T_p);
            rho = gsw_rho(SA,CT,P);
            rho_save(it,:) = rho;

            % [v7-G] shallow / coastal guard
            if sum(isfinite(rho)) < 8
                MLD_m(it) = NaN; ILD_m(it) = NaN; BLT_m(it) = NaN;
                continue;
            end
            [N2,~] = gsw_Nsquared(SA,CT,P,lat_tc(it));
            N2_save(it,:) = max(N2.',1e-9);

            MLD_m(it) = mld_density(z, rho, DELTA_RHO);
            ILD_m(it) = mld_temperature(z, T_p, DELTA_T);
            BLT_m(it) = ILD_m(it) - MLD_m(it);

            [~,imb] = min(abs(z - MLD_m(it)));
            ML_Buoyancy(it) = g*(rho(imb)-rho(1))/rho0;
        end
        if ~isfinite(MLD_m(it)), continue; end
        % [v7-B/C] POINTWISE NIKE
        gt = find(got);
        X  = [cos(f_val*tsec(gt)) sin(f_val*tsec(gt))];

        Ecoh = nan(ncell, nz);
        Ebrd = nan(ncell, nz);
        ph_surface = nan(ncell,1);

        for k = 1:nz
            Uk = U(gt,:,k);   Vk = V(gt,:,k);
            good_t = ~all(~isfinite(Uk),2) & ~all(~isfinite(Vk),2);
            Uk = Uk(good_t,:);  Vk = Vk(good_t,:);  Xk = X(good_t,:);
            if size(Uk,1) < 11, continue; end
            Uk = detrend_cols(Uk);
            Vk = detrend_cols(Vk);
            ok = all(isfinite(Uk),1) & all(isfinite(Vk),1);
            if ~any(ok), continue; end

            if DO_POINTWISE
                CU = Xk \ Uk(:,ok);
                CV = Xk \ Vk(:,ok);
                Au = hypot(CU(1,:),CU(2,:));
                Av = hypot(CV(1,:),CV(2,:));
                rk = rho_save(it,k); if ~isfinite(rk), rk = rho0; end
                Ecoh(ok,k) = 0.25*rk*(Au.^2 + Av.^2).';
                Ebrd(ok,k) = 0.5 *rk*mean(Uk(:,ok).^2 + Vk(:,ok).^2, 1).';
                if k == 1
                    ph_surface(ok) = atan2(CU(2,:),CU(1,:)).';
                end
            else
                um = mean(Uk(:,ok),2);  vm = mean(Vk(:,ok),2);
                cu = Xk\um;  cv = Xk\vm;
                Au = hypot(cu(1),cu(2));  Av = hypot(cv(1),cv(2));
                rk = rho_save(it,k); if ~isfinite(rk), rk = rho0; end
                Ecoh(ok,k) = 0.25*rk*(Au^2 + Av^2);
                Ebrd(ok,k) = 0.5 *rk*mean(um.^2 + vm.^2);
            end
        end

        % spatial coherence diagnostic
        pv = ph_surface(isfinite(ph_surface) & mask(:));
        if ~isempty(pv)
            Rbar = abs(mean(exp(1i*pv)));
            phase_sd_deg(it) = sqrt(max(-2*log(Rbar),0))*180/pi;
        end

        % ---- spatial means, then depth integral ----
        mk = mask(:);  mR = mask_R(:);  mL = mask_L(:);
        n_cells_used(it) = sum(mk);

        % disk-mean surface velocity time series (rotary spectra, Fig 6/8)
        u_surf_save(it,:) = squeeze(mean(U(:,mk,1), 2, 'omitnan')).';
        v_surf_save(it,:) = squeeze(mean(V(:,mk,1), 2, 'omitnan')).';

        prof_coh   = squeeze(mean(Ecoh(mk,:), 1, 'omitnan')).';
        prof_brd   = squeeze(mean(Ebrd(mk,:), 1, 'omitnan')).';
        prof_cohR  = squeeze(mean(Ecoh(mR,:), 1, 'omitnan')).';
        prof_cohL  = squeeze(mean(Ecoh(mL,:), 1, 'omitnan')).';
        prof_brdR  = squeeze(mean(Ebrd(mR,:), 1, 'omitnan')).';
        prof_brdL  = squeeze(mean(Ebrd(mL,:), 1, 'omitnan')).';

        % [v7-D] baseline subtraction — SURFACE/ML ONLY
        if ~isempty(NIKE_base_z) && numel(NIKE_base_z)==nz
            bl = NIKE_base_z(:);
            bl(z > MLD_m(it)) = 0;              % no subtraction below the ML
            prof_coh  = max(prof_coh  - bl, 0);
            prof_cohR = max(prof_cohR - bl, 0);
            prof_cohL = max(prof_cohL - bl, 0);
        end

        NIKE_coherent_z_save(it,:)   = prof_coh.';
        NIKE_broadband_z_save(it,:)  = prof_brd.';
        NIKE_coherent_z_right(it,:)  = prof_cohR.';
        NIKE_coherent_z_left(it,:)   = prof_cohL.';
        NIKE_broadband_z_right(it,:) = prof_brdR.';
        NIKE_broadband_z_left(it,:)  = prof_brdL.';

        z_col = z(:);
        NIKE_coherent_Jm2(it)        = safe_trapz(z_col, prof_coh);
        NIKE_broadband_Jm2(it)       = safe_trapz(z_col, prof_brd);
        NIKE_coherent_right_Jm2(it)  = safe_trapz(z_col, prof_cohR);
        NIKE_coherent_left_Jm2(it)   = safe_trapz(z_col, prof_cohL);
        NIKE_broadband_right_Jm2(it) = safe_trapz(z_col, prof_brdR);
        NIKE_broadband_left_Jm2(it)  = safe_trapz(z_col, prof_brdL);

        Coherence_Ratio(it) = NIKE_coherent_Jm2(it)/(NIKE_broadband_Jm2(it)+eps);
        if Coherence_Ratio(it) > 1.0
            warning('R_C = %.3f > 1 at it=%d', Coherence_Ratio(it), it);
        end

        % ML / deep split (total + right/left)
        vv = isfinite(prof_coh);
        iml = (z_col <= MLD_m(it)) & vv;
        if sum(iml) >= 2
            NIKE_above_MLD_Jm2(it) = trapz(z_col(iml), prof_coh(iml));
            NIKE_below_MLD_Jm2(it) = NIKE_coherent_Jm2(it) - NIKE_above_MLD_Jm2(it);
        end
        ivR = isfinite(prof_cohR) & (z_col > MLD_m(it));
        if sum(ivR) >= 2, NIKE_deep_right_Jm2(it) = trapz(z_col(ivR), prof_cohR(ivR)); end
        ivL = isfinite(prof_cohL) & (z_col > MLD_m(it));
        if sum(ivL) >= 2, NIKE_deep_left_Jm2(it) = trapz(z_col(ivL), prof_cohL(ivL)); end

        if sum(vv) >= 2
            NIKE_centroid_m(it) = trapz(z_col(vv), z_col(vv).*prof_coh(vv)) / ...
                                  (NIKE_coherent_Jm2(it)+eps);
        end
        i50 = (z_col<=50) & isfinite(prof_brd);
        if sum(i50) >= 2
            NIKE_surface_frac(it) = trapz(z_col(i50), prof_brd(i50)) / ...
                                    (NIKE_broadband_Jm2(it)+eps);
        end

        % asymmetry (coherent, bilateral 10% guard)
        if isfinite(NIKE_coherent_right_Jm2(it)) && isfinite(NIKE_coherent_left_Jm2(it)) && ...
           NIKE_coherent_left_Jm2(it)  > 0.10*NIKE_coherent_Jm2(it) && ...
           NIKE_coherent_right_Jm2(it) > 0.10*NIKE_coherent_Jm2(it)
            Asymmetry_Ratio(it) = min(NIKE_coherent_right_Jm2(it)/NIKE_coherent_left_Jm2(it), 10);
        end
        %  WIND STRESS / WIND WORK
        i0 = find(twin_hours==0);
        [~,ie] = min(abs(t_e - time_tc(it)));
        tx = tau_x_all(:,:,ie).';  ty = tau_y_all(:,:,ie).';
        dE = distdim(distance(lat_tc(it),lon_tc(it),LatE,LonE),'deg','km');
        dyE = (LatE-lat_tc(it))*111;  dxE = (LonE-lon_tc(it))*111*cosd(lat_tc(it));
        crE = dxE*sin(h_rad) - dyE*cos(h_rad);
        mE  = dE <= R_NI;  mER = mE & crE>=0;  mEL = mE & crE<0;

        tmag = hypot(tx,ty);
        WindStress_Nm2(it)       = mean(tmag(mE), 'omitnan');
        WindStress_right_Nm2(it) = mean(tmag(mER),'omitnan');
        WindStress_left_Nm2(it)  = mean(tmag(mEL),'omitnan');

        txh = interp2(LonE,LatE,tx,LonSub,LatSub,'linear');
        tyh = interp2(LonE,LatE,ty,LonSub,LatSub,'linear');
        u0 = reshape(U(i0,:,1), ny, nx);
        v0 = reshape(V(i0,:,1), ny, nx);
        wp = txh.*u0 + tyh.*v0;
        WindPower_Wm2(it)       = mean(wp(mask), 'omitnan');
        WindPower_right_Wm2(it) = mean(wp(mask_R),'omitnan');
        WindPower_left_Wm2(it)  = mean(wp(mask_L),'omitnan');

        % PF06 form: tau . u_I
        Uk = detrend_cols(U(gt,:,1));  Vk = detrend_cols(V(gt,:,1));
        ok = all(isfinite(Uk),1) & all(isfinite(Vk),1);
        if any(ok)
            CU = X \ Uk(:,ok);  CV = X \ Vk(:,ok);
            uI0 = nan(ncell,1); vI0 = nan(ncell,1);
            c0  = [cos(f_val*tsec(i0)) sin(f_val*tsec(i0))];
            uI0(ok) = (c0*CU).';   vI0(ok) = (c0*CV).';
            uI0m = reshape(uI0, ny, nx);  vI0m = reshape(vI0, ny, nx);
            wpNI = txh.*uI0m + tyh.*vI0m;
            WindPower_NI_Wm2(it) = mean(wpNI(mask),'omitnan');

            uI_amp_surf(it) = mean(hypot(hypot(CU(1,:),CU(2,:)), ...
                                         hypot(CV(1,:),CV(2,:))),'omitnan');
            a1 = atan2(mean(tyh(mask),'omitnan'), mean(txh(mask),'omitnan'));
            a2 = atan2(mean(vI0m(mask),'omitnan'), mean(uI0m(mask),'omitnan'));
            phase_tau_uI_deg(it) = (mod(a1-a2+pi,2*pi)-pi)*180/pi;
        end

        % vorticity, strain, effective Coriolis (R1#2 / R2#7)
        dym = 111e3*mean(diff(lat_g(iy)));
        dxm = 111e3*cosd(lat_tc(it))*mean(diff(lon_g(ix)));
        [du_dx, du_dy] = gradient(u0, dxm, dym);
        [dv_dx, dv_dy] = gradient(v0, dxm, dym);
        zeta = mean(dv_dx - du_dy, 'all', 'omitnan');
        strain_n = mean(du_dx - dv_dy, 'all', 'omitnan');
        strain_s = mean(dv_dx + du_dy, 'all', 'omitnan');
        strain_rate(it) = hypot(strain_n, strain_s);
        zeta_save(it)   = zeta;
        f_eff(it)       = f_val + 0.5*zeta;

        % Richardson number
        u0p = squeeze(mean(U(i0,mk,:),2,'omitnan'));
        v0p = squeeze(mean(V(i0,mk,:),2,'omitnan'));
        du = diff(u0p)./diff(z);  dv = diff(v0p)./diff(z);
        S2raw = du.^2 + dv.^2;
        S2    = max(S2raw, 1e-7);
        S2_save(it,:) = S2.';
        N2row = N2_save(it,:).';
        Rirow = N2row./S2;
        bad = (N2row<1e-7)|(S2raw<1e-8)|~isfinite(Rirow)|(Rirow<1e-3)|(Rirow>1e4);
        Rirow(bad) = NaN;
        Ri_save(it,:) = Rirow.';
        [~,iml2] = min(abs(z-MLD_m(it)));
        iml2 = min(nz-1, max(1,iml2));
        Ri_Bulk_MLD(it) = Rirow(iml2);

        % upper-100 m robust Ri percentiles (headline; avoids single bad cell)
        zc = z(1:nz-1);
        up = zc <= 100;
        Ri_u = Rirow(up); Ri_u = Ri_u(isfinite(Ri_u));
        if numel(Ri_u) >= 5
            Ri_min_100m(it) = prctile(Ri_u,5);
            Ri_med_100m(it) = median(Ri_u);
        elseif ~isempty(Ri_u)
            Ri_min_100m(it) = min(Ri_u);
            Ri_med_100m(it) = median(Ri_u);
        end

        % ---- near-inertial PE and PE/KE ratio ----
        % PE from isopycnal heave xi = (d u_I / dz)/N ; PE = 0.5*rho*N^2*xi^2.
        % KE here = coherent NIKE profile. Ratio -> ~0.5 forced slab, ~1 free NIW.
        prof_ke = prof_coh(:);                      % J/m^3 coherent KE profile
        if any(isfinite(prof_ke))
            uI_col = sqrt(max(prof_ke,0)/(0.5*rho0));   % proxy inertial speed profile
            duIdz  = [diff(uI_col)./diff(z); NaN];
            Nmid   = sqrt(max(N2_save(it,:).',1e-9));
            Nmid(end+1) = Nmid(end);
            xi     = duIdz ./ (Nmid + eps);             % isopycnal displacement
            pe_z   = 0.5*rho0*(Nmid.^2).*(xi.^2);       % J/m^3
            NIW_PE_Jm2(it) = safe_trapz(z, pe_z);
            PE_to_KE(it)   = NIW_PE_Jm2(it) / (NIKE_coherent_Jm2(it) + eps);
            if PE_to_KE(it) > 5, PE_to_KE(it) = NaN; end   % guard nonphysical
        end

        % SST / OHC / Stokes / Langmuir
        SST_C(it)   = sample_aux(sst_all, LonS, LatS, t_sst, time_tc(it), lon_tc(it), lat_tc(it), -273.15);
        OHC_Jm2(it) = sample_aux(ohc_all, LonO, LatO, t_ohc, time_tc(it), lon_tc(it), lat_tc(it), 0)*1e7;
        us = sample_aux(ust_all, LonSt, LatSt, t_st, time_tc(it), lon_tc(it), lat_tc(it), 0);
        vs = sample_aux(vst_all, LonSt, LatSt, t_st, time_tc(it), lon_tc(it), lat_tc(it), 0);
        Stokes_ms(it) = hypot(us,vs);
        ustar = sqrt(WindStress_Nm2(it)/rho0);
        if isfinite(Stokes_ms(it)) && isfinite(ustar)
            La_t(it) = sqrt(ustar/(Stokes_ms(it)+eps));
        end

        if mod(it,5)==0 || it==nt
            fprintf('  it %3d/%d  NIKE=%6.2f kJ/m2  MLD=%5.1f  BLT=%5.1f  phSD=%5.1f deg\n', ...
                it, nt, NIKE_coherent_Jm2(it)/1000, MLD_m(it), BLT_m(it), phase_sd_deg(it));
        end
    end

    %% WIND WORK BUDGET
    t_sec = seconds(time_tc - time_tc(1));
    WindPower_signed = WindPower_NI_Wm2;
    WindPower_signed(~isfinite(WindPower_signed)) = 0;
    WindWork_signed_Jm2 = cumtrapz(t_sec, WindPower_signed);

    WindPower_pos = WindPower_signed;  WindPower_pos(WindPower_pos<0) = 0;
    WindWork_pos_Jm2 = cumtrapz(t_sec, WindPower_pos);

    dNIKE_dt_Wm2 = gradient(NIKE_broadband_Jm2, t_sec);

    %% REPORT
    fprintf('\n=================== %s v7 RESULTS ===================\n', cfg.name);
    fprintf('Peak NIKE coherent     : %.2f kJ/m2\n', max(NIKE_coherent_Jm2)/1000);
    fprintf('Peak NIKE broadband    : %.2f kJ/m2\n', max(NIKE_broadband_Jm2)/1000);
    fprintf('Baseline subtracted    : %.2f kJ/m2\n', NIKE_base_col/1000);
    fprintf('Mean coherence ratio   : %.3f\n', mean(Coherence_Ratio,'omitnan'));
    fprintf('Mean MLD (density)     : %.1f m\n', mean(MLD_m,'omitnan'));
    fprintf('Mean ILD (temperature) : %.1f m\n', mean(ILD_m,'omitnan'));
    fprintf('Mean BLT               : %.1f m\n', mean(BLT_m,'omitnan'));
    fprintf('Below-ML fraction      : %.1f %%\n', ...
        100*mean(NIKE_below_MLD_Jm2./(NIKE_coherent_Jm2+eps),'omitnan'));
    fprintf('Mean spatial phase sd  : %.1f deg\n', mean(phase_sd_deg,'omitnan'));
    fprintf('Mean |u_I| surface     : %.4f m/s\n', mean(uI_amp_surf,'omitnan'));
    fprintf('Mean |zeta|/f          : %.3f\n', mean(abs(zeta_save)./abs(f_tc),'omitnan'));
    fprintf('Mean strain rate       : %.3e 1/s\n', mean(strain_rate,'omitnan'));
    fprintf('\n--- WIND WORK ---\n');
    fprintf('Mean WindStress        : %.4f N/m2\n', mean(WindStress_Nm2,'omitnan'));
    fprintf('Mean tau.u_total       : %.4f W/m2\n', mean(WindPower_Wm2,'omitnan'));
    fprintf('Mean tau.u_I  (PF06)   : %.4f W/m2\n', mean(WindPower_NI_Wm2,'omitnan'));
    iw = find(isfinite(WindWork_signed_Jm2),1,'last');
    fprintf('Cum wind work (signed) : %+.2f kJ/m2   [net, end of record]\n', WindWork_signed_Jm2(iw)/1000);
    fprintf('Cum wind work (pos)    : %+.2f kJ/m2   [gross]\n', WindWork_pos_Jm2(iw)/1000);
    fprintf('  retention (signed/pos): %.0f %%\n', 100*WindWork_signed_Jm2(iw)/max(WindWork_pos_Jm2(iw),eps));
    fprintf('\n--- PHASE (R1 #3) ---\n');
    pvp = phase_tau_uI_deg(isfinite(phase_tau_uI_deg));
    fprintf('median phase tau vs uI : %.1f deg\n', median(pvp));
    fprintf('frac within +/-45 deg  : %.1f %%\n', 100*sum(abs(pvp)<45)/numel(pvp));
    fprintf('\n--- ASYMMETRY ---\n');
    fprintf('Mean R/L ratio         : %.2f\n', mean(Asymmetry_Ratio,'omitnan'));
    fprintf('Median R/L ratio       : %.2f\n', median(Asymmetry_Ratio,'omitnan'));
    fprintf('=====================================================\n');

    %% SAVE (exhaustive)
    save(cfg.out_mat, ...
        'z','time_tc','lat_tc','lon_tc','heading_deg_save', ...
        'NIKE_coherent_Jm2','NIKE_broadband_Jm2', ...
        'NIKE_coherent_z_save','NIKE_broadband_z_save', ...
        'NIKE_coherent_z_right','NIKE_coherent_z_left', ...
        'NIKE_broadband_z_right','NIKE_broadband_z_left', ...
        'NIKE_coherent_right_Jm2','NIKE_coherent_left_Jm2', ...
        'NIKE_broadband_right_Jm2','NIKE_broadband_left_Jm2', ...
        'NIKE_deep_right_Jm2','NIKE_deep_left_Jm2', ...
        'NIKE_above_MLD_Jm2','NIKE_below_MLD_Jm2', ...
        'NIKE_centroid_m','NIKE_surface_frac', ...
        'MLD_m','ILD_m','BLT_m','ML_Buoyancy', ...
        'rho_save','N2_save','Ri_save','S2_save', ...
        'Ri_Bulk_MLD','Ri_min_100m','Ri_med_100m', ...
        'SST_C','OHC_Jm2','Stokes_ms','La_t', ...
        'WindStress_Nm2','WindStress_right_Nm2','WindStress_left_Nm2', ...
        'WindPower_Wm2','WindPower_right_Wm2','WindPower_left_Wm2', ...
        'WindPower_NI_Wm2','WindPower_signed','WindPower_pos', ...
        'WindWork_signed_Jm2','WindWork_pos_Jm2','dNIKE_dt_Wm2', ...
        'Coherence_Ratio','Asymmetry_Ratio', ...
        'f_tc','f_eff','zeta_save','strain_rate','Inertial_hr', ...
        'phase_tau_uI_deg','uI_amp_surf','phase_sd_deg','n_cells_used', ...
        'u_surf_save','v_surf_save','PE_to_KE','NIW_PE_Jm2', ...
        'NIKE_base_z','NIKE_base_col', ...
        'WIN_HR','DO_POINTWISE','BASELINE_DAYS','R_NI','DELTA_RHO','DELTA_T', ...
        '-v7.3');
    fprintf('saved: %s\n', cfg.out_mat);


end   % STORM_LIST loop

fprintf('\n============================================================\n');
fprintf(' ALL STORMS COMPLETE\n');
fprintf('============================================================\n');

%% LOCAL FUNCTIONS
function hr = compute_heading(lat, lon, win)
    n = numel(lat);  hr = nan(n,1);
    for i = 1:n
        if i==1,      dla = lat(2)-lat(1);      dlo = lon(2)-lon(1);
        elseif i==n,  dla = lat(n)-lat(n-1);    dlo = lon(n)-lon(n-1);
        else,         dla = lat(i+1)-lat(i-1);  dlo = lon(i+1)-lon(i-1);
        end
        hr(i) = atan2(dla*111, dlo*111*cosd(lat(i)));
    end
    hr = movmean(unwrap(hr), win, 'omitnan');
end

function [times, fidx, tidx] = index_time_axis(files)
    times = datetime.empty(0,1);  fidx = [];  tidx = [];
    for i = 1:numel(files)
        fn = fullfile(files(i).folder, files(i).name);
        try
            tv = ncread(fn,'time');
            ti = ncinfo(fn,'time');
            un = '';
            for a = 1:numel(ti.Attributes)
                if strcmpi(ti.Attributes(a).Name,'units')
                    un = ti.Attributes(a).Value;
                end
            end
            tok = regexp(un,'hours since\s+([\d\-]+)[ T]([\d:]+)','tokens','once');
            if isempty(tok)
                tok = regexp(un,'hours since\s+([\d\-]+)','tokens','once');
                t0 = datetime(tok{1},'InputFormat','yyyy-MM-dd');
            else
                t0 = datetime([tok{1} ' ' tok{2}],'InputFormat','yyyy-MM-dd HH:mm:ss');
            end
            tt = t0 + hours(double(tv(:)));
            times = [times; tt];             %#ok<AGROW>
            fidx  = [fidx;  repmat(i,numel(tt),1)];  %#ok<AGROW>
            tidx  = [tidx;  (1:numel(tt)).'];        %#ok<AGROW>
        catch
        end
    end
    [times, o] = sort(times);
    fidx = fidx(o);  tidx = tidx(o);
end

function [fi, ti] = find_time(times, fidx, tidx, target)
    fi = []; ti = [];
    if isempty(times), return; end
    [d, k] = min(abs(times - target));
    if d <= minutes(90)
        fi = fidx(k);  ti = tidx(k);
    end
end

function Y = detrend_cols(X)
    Y = X;  n = size(X,1);  t = (1:n).';
    for c = 1:size(X,2)
        gg = isfinite(X(:,c));
        if sum(gg) >= 3
            p = polyfit(t(gg), X(gg,c), 1);
            Y(gg,c) = X(gg,c) - polyval(p, t(gg));
        end
    end
end

function d = mld_density(z, rho, dcrit)
    [~,i0] = min(abs(z-10));
    while i0 <= numel(rho) && ~isfinite(rho(i0)), i0 = i0+1; end
    if i0 > numel(rho), d = NaN; return; end
    k = find(z > z(i0) & rho > rho(i0)+dcrit, 1, 'first');
    if isempty(k)
        k = find(isfinite(rho),1,'last');  d = z(k);
    else
        if k > 1
            f = (rho(i0)+dcrit - rho(k-1)) / (rho(k)-rho(k-1)+eps);
            d = z(k-1) + f*(z(k)-z(k-1));
        else
            d = z(k);
        end
    end
end

function d = mld_temperature(z, T, dT)
    [~,i0] = min(abs(z-10));
    while i0 <= numel(T) && ~isfinite(T(i0)), i0 = i0+1; end
    if i0 > numel(T), d = NaN; return; end
    k = find(z > z(i0) & T < T(i0)-dT, 1, 'first');
    if isempty(k)
        k = find(isfinite(T),1,'last');  d = z(k);
    else
        if k > 1
            f = (T(i0)-dT - T(k-1)) / (T(k)-T(k-1)+eps);
            d = z(k-1) + f*(z(k)-z(k-1));
        else
            d = z(k);
        end
    end
end

function I = safe_trapz(z, p)
    v = isfinite(p);
    if sum(v) >= 2, I = trapz(z(v), p(v)); else, I = NaN; end
end

function [A, Lo, La, tt] = load_aux(file, vnames, t_ref, t_fallback)
    if ischar(vnames), vnames = {vnames}; end
    A = [];
    for i = 1:numel(vnames)
        try, A = ncread(file, vnames{i}); break; catch, end
    end
    try
        lo = ncread(file,'lon');  la = ncread(file,'lat');
    catch
        lo = ncread(file,'longitude'); la = ncread(file,'latitude');
    end
    [Lo,La] = meshgrid(lo,la);
    try
        tt = t_ref + hours(ncread(file,'time'));
    catch
        tt = t_fallback;
    end
end

function val = sample_aux(A, Lo, La, tt, target, lon_q, lat_q, offset)
    val = NaN;
    if isempty(A), return; end
    if ndims(A) == 3
        [~,k] = min(abs(tt - target));
        Ssl = A(:,:,k).';
    else
        Ssl = A.';
    end
    if min(Lo(:)) >= 0, lon_q = mod(lon_q,360); end
    val = interp2(Lo, La, Ssl, lon_q, lat_q) + offset;
    if isnan(val)
        val = interp2(Lo, La, Ssl, lon_q, lat_q, 'nearest') + offset;
    end
end

function [prof, col] = nike_at_point(t_center, lat_c, lon_c, lon_g, lat_g, z, idz, ...
                                     R_NI, twin_hours, tsec, times, fidx, tidx, ...
                                     files, rho0, pointwise)
    nz = numel(z);  nw = numel(twin_hours);
    prof = nan(nz,1);  col = NaN;
    Omega = 7.2921e-5;
    f_val = 2*Omega*sind(lat_c);
    if abs(f_val) < 1e-5, f_val = sign(f_val)*1e-5; end

    ix = find(lon_g >= lon_c-2 & lon_g <= lon_c+2);
    iy = find(lat_g >= lat_c-2 & lat_g <= lat_c+2);
    if isempty(ix)||isempty(iy), return; end
    [LonSub,LatSub] = meshgrid(lon_g(ix), lat_g(iy));
    dist = distdim(distance(lat_c,lon_c,LatSub,LonSub),'deg','km');
    mk = dist(:) <= R_NI;
    nx = numel(ix); ny = numel(iy); ncell = nx*ny;

    U = nan(nw,ncell,nz); V = nan(nw,ncell,nz); got = false(nw,1);
    for jt = 1:nw
        tgt = t_center + hours(twin_hours(jt));
        [fi,ti] = find_time(times, fidx, tidx, tgt);
        if isempty(fi), continue; end
        fn = fullfile(files(fi).folder, files(fi).name);
        try
            ub = ncread(fn,'water_u',[ix(1) iy(1) idz(1) ti],[nx ny nz 1]);
            vb = ncread(fn,'water_v',[ix(1) iy(1) idz(1) ti],[nx ny nz 1]);
        catch, continue; end
        for k = 1:nz
            su = ub(:,:,k).'; sv = vb(:,:,k).';
            U(jt,:,k) = su(:).'; V(jt,:,k) = sv(:).';
        end
        got(jt) = true;
    end
    if sum(got) < 11, return; end

    gt = find(got);
    X = [cos(f_val*tsec(gt)) sin(f_val*tsec(gt))];
    for k = 1:nz
        Uk = detrend_cols(U(gt,:,k));  Vk = detrend_cols(V(gt,:,k));
        ok = all(isfinite(Uk),1) & all(isfinite(Vk),1) & mk.';
        if ~any(ok), continue; end
        if pointwise
            CU = X\Uk(:,ok); CV = X\Vk(:,ok);
            Au = hypot(CU(1,:),CU(2,:)); Av = hypot(CV(1,:),CV(2,:));
            prof(k) = mean(0.25*rho0*(Au.^2+Av.^2));
        else
            um = mean(Uk(:,ok),2); vm = mean(Vk(:,ok),2);
            cu = X\um; cv = X\vm;
            prof(k) = 0.25*rho0*(hypot(cu(1),cu(2))^2 + hypot(cv(1),cv(2))^2);
        end
    end
    col = safe_trapz(z(:), prof);
end

function cfg = get_storm_config(storm_id)
    base_hycom  = [fullfile(DATA_DIR,'HYCOM') filesep];
    base_track  = [fullfile(DATA_DIR,'IMD') filesep];
    base_stress = [fullfile(DATA_DIR,'Stress') filesep];
    base_curr   = [fullfile(DATA_DIR,'Currents') filesep];
    base_ohc    = [fullfile(DATA_DIR,'INCOIS_OHC') filesep];

    switch upper(storm_id)
        case 'KYARR'
            cfg.name='KYARR'; cfg.year_ref=2019; cfg.heading_smooth_window=3;
            cfg.uv_path=[base_hycom 'Kyarr_UV3Z\']; cfg.ts_path=[base_hycom 'Kyarr_TS3Z\'];
            cfg.track_path=[base_track 'Kyarr.csv'];
            cfg.stress_file=[base_stress 'Kyarr\era5_hourly_stress_KYARR_FINAL.nc'];
            cfg.stokes_file=[base_curr 'Kyarr\ocean_stokes_kyarr.nc'];
            cfg.sst_file   =[base_curr 'Kyarr\ocean_sst_kyarr.nc'];
            cfg.ohc_file   =[base_ohc  'Kyarr\2019_20_ohc700.nc'];
            cfg.out_mat='KYARR_V17.mat';
        case 'AMPHAN'
            cfg.name='AMPHAN'; cfg.year_ref=2020; cfg.heading_smooth_window=3;
            cfg.uv_path=[base_hycom 'Amphan_UV3Z\']; cfg.ts_path=[base_hycom 'Amphan_TS3Z\'];
            cfg.track_path=[base_track 'Amphan.csv'];
            cfg.stress_file=[base_stress 'Amphan\era5_hourly_stress_FINAL.nc'];
            cfg.stokes_file=[base_curr 'Amphan\ocean_stokes.nc'];
            cfg.sst_file   =[base_curr 'Amphan\ocean_sst.nc'];
            cfg.ohc_file   =[base_ohc  'Amphan\ohc700_May2020.nc'];
            cfg.out_mat='AMPHAN_V17.mat';
        case 'FANI'
            cfg.name='FANI'; cfg.year_ref=2019; cfg.heading_smooth_window=5;
            cfg.uv_path=[base_hycom 'Fani_UV3Z\']; cfg.ts_path=[base_hycom 'Fani_TS3Z\'];
            cfg.track_path=[base_track 'Fani.csv'];
            cfg.stress_file=[base_stress 'Fani\era5_hourly_stress_FANI_FINAL.nc'];
            cfg.stokes_file=[base_curr 'Fani\ocean_stokes_fani.nc'];
            cfg.sst_file   =[base_curr 'Fani\ocean_sst_fani.nc'];
            cfg.ohc_file   =[base_ohc  'Fani\2019_20_ohc700.nc'];
            cfg.out_mat='FANI_V17.mat';
        case 'TAUKTAE'
            cfg.name='TAUKTAE'; cfg.year_ref=2021; cfg.heading_smooth_window=5;
            cfg.uv_path=[base_hycom 'Tauktae_UV3Z\']; cfg.ts_path=[base_hycom 'Tauktae_TS3Z\'];
            cfg.track_path=[base_track 'Tauktae.csv'];
            cfg.stress_file=[base_stress 'Tauktae\era5_hourly_stress_TAUKTAE_FINAL.nc'];
            cfg.stokes_file=[base_curr 'tauktae\ocean_stokes_tauktae.nc'];
            cfg.sst_file   =[base_curr 'tauktae\ocean_sst_tauktae.nc'];
            cfg.ohc_file   =[base_ohc  'tauktae\2021_22_ohc700.nc'];
            cfg.out_mat='TAUKTAE_V17.mat';
        otherwise
            error('Unknown STORM_ID: %s', storm_id);
    end
end
