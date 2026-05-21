clc; clear; close all; warning('off','all');

% Unified cyclone NIKE core. Set STORM_ID below to run one of the four storms.
% Output: <STORM>_V16_SUPERCHARGED.mat + .xlsx, and a *_SUMMARY_v6.txt report.

STORM_ID = 'TAUKTAE';   % KYARR | AMPHAN | FANI | TAUKTAE

cfg = get_storm_config(STORM_ID);

fprintf('\n=== %s | year_ref=%d | heading_smooth=%d | R_NI=150 km ===\n\n', ...
    cfg.name, cfg.year_ref, cfg.heading_smooth_window);

% Paths and constants
GSW_PATH = fullfile(pwd, 'GSW-Matlab');     % adjust to your local GSW install
if exist(GSW_PATH,'dir'), addpath(genpath(GSW_PATH)); end

track_path   = cfg.track_path;
hycom_path   = cfg.hycom_path;
stress_file  = cfg.stress_file;
wind_file    = cfg.wind_file;
stokes_file  = cfg.stokes_file;
sst_file     = cfg.sst_file;
ohc_file     = cfg.ohc_file;
year_ref     = cfg.year_ref;
out_mat      = cfg.out_mat;
out_xlsx     = cfg.out_xlsx;
hycom_prefix = cfg.hycom_prefix;
heading_window = cfg.heading_smooth_window;

g         = 9.81;
Omega     = 7.2921e-5;
rho0      = 1025;
DELTA_RHO = 0.03;
Hmax      = 700;
R_NI      = 150;                  % km, fixed mask radius (Brizuela 2023)
twin_hours = -24:3:24;
tsec      = twin_hours(:)*3600;
dt_win    = 3*3600;

% Track and heading
Ttrk    = readtable(track_path);
time_tc = datetime(Ttrk.year,Ttrk.month,Ttrk.day,Ttrk.hour,0,0);
lat_tc  = Ttrk.lat;
lon_tc  = Ttrk.lon;
nt      = length(time_tc);

heading_rad = nan(nt,1);
for it = 1:nt
    if it == 1
        dlat_h = lat_tc(2) - lat_tc(1);
        dlon_h = lon_tc(2) - lon_tc(1);
    elseif it == nt
        dlat_h = lat_tc(nt) - lat_tc(nt-1);
        dlon_h = lon_tc(nt) - lon_tc(nt-1);
    else
        dlat_h = lat_tc(it+1) - lat_tc(it-1);
        dlon_h = lon_tc(it+1) - lon_tc(it-1);
    end
    dy_h = dlat_h * 111;
    dx_h = dlon_h * 111 * cosd(lat_tc(it));
    heading_rad(it) = atan2(dy_h, dx_h);
end
heading_rad = unwrap(heading_rad);
heading_rad = movmean(heading_rad, heading_window, 'omitnan');   % 3 straight, 5 recurving

% HYCOM grid reference
f0 = dir(fullfile(hycom_path,'*.nc'));
if isempty(f0), error('No HYCOM files found in %s', hycom_path); end
fn_ref = fullfile(f0(1).folder,f0(1).name);
try
    lon_g = ncread(fn_ref,'lon'); lat_g = ncread(fn_ref,'lat');
catch
    lon_g = ncread(fn_ref,'longitude'); lat_g = ncread(fn_ref,'latitude');
end
depth = ncread(fn_ref,'depth');
idz   = find(depth<=Hmax);
z     = depth(idz);
nz    = length(z);
[~,idx10] = min(abs(z-10));

% Load external data: stress, SST, OHC, Stokes
t_ref_dt  = datetime(year_ref,1,1);
t_e       = t_ref_dt + hours(ncread(stress_file,'time'));
tau_x_all = ncread(stress_file,'iews');
tau_y_all = ncread(stress_file,'inss');
try
    lon_e = ncread(stress_file,'lon'); lat_e = ncread(stress_file,'lat');
catch
    lon_e = ncread(stress_file,'longitude'); lat_e = ncread(stress_file,'latitude');
end
[LonE,LatE] = meshgrid(lon_e, lat_e);

try
    lon_s = ncread(sst_file,'lon'); lat_s = ncread(sst_file,'lat');
catch
    lon_s = ncread(sst_file,'longitude'); lat_s = ncread(sst_file,'latitude');
end
sst_all      = ncread(sst_file,'sst');
[LonS, LatS] = meshgrid(lon_s, lat_s);
try
    t_sst = t_ref_dt + hours(ncread(sst_file,'time'));
catch
    t_sst = t_e;
end

try
    lon_o = ncread(ohc_file,'lon'); lat_o = ncread(ohc_file,'lat');
catch
    lon_o = ncread(ohc_file,'longitude'); lat_o = ncread(ohc_file,'latitude');
end
try
    ohc_all = ncread(ohc_file,'ohc');
catch
    ohc_all = ncread(ohc_file,'ohc700');
end
[LonO, LatO] = meshgrid(lon_o, lat_o);
try
    t_ohc = t_ref_dt + hours(ncread(ohc_file,'time'));
catch
    t_ohc = t_e;
end

try
    lon_st = ncread(stokes_file,'lon'); lat_st = ncread(stokes_file,'lat');
catch
    lon_st = ncread(stokes_file,'longitude'); lat_st = ncread(stokes_file,'latitude');
end
try
    ust_all = ncread(stokes_file,'u_stokes');
    vst_all = ncread(stokes_file,'v_stokes');
catch
    ust_all = ncread(stokes_file,'ust');
    vst_all = ncread(stokes_file,'vst');
end
[LonSt, LatSt] = meshgrid(lon_st, lat_st);
try
    t_st = t_ref_dt + hours(ncread(stokes_file,'time'));
catch
    t_st = t_e;
end

% Storage
NIKE_coherent_Jm2          = nan(nt,1);
NIKE_broadband_Jm2         = nan(nt,1);
NIKE_surface_frac          = nan(nt,1);
NIKE_centroid_m            = nan(nt,1);
NIKE_above_MLD_Jm2         = nan(nt,1);
NIKE_below_MLD_Jm2         = nan(nt,1);
WindStress_Nm2             = nan(nt,1);
SST_C                      = nan(nt,1);
OHC_Jm2                    = nan(nt,1);
Stokes_ms                  = nan(nt,1);
La_t                       = nan(nt,1);
ML_Buoyancy                = nan(nt,1);
MLD_m                      = nan(nt,1);
Inertial_hr                = nan(nt,1);
f_tc                       = nan(nt,1);
f_eff                      = nan(nt,1);
Ri_Bulk_MLD                = nan(nt,1);
Coherence_Ratio            = nan(nt,1);
NIKE_WKB_Jm2               = nan(nt,1);
NIW_PE_Jm2                 = nan(nt,1);
PE_to_KE_Ratio             = nan(nt,1);
Inferred_Cgz_m_day         = nan(nt,1);
Centroid_Trend_m_day       = nan(nt,1);
Entrainment_Velocity_m_day = nan(nt,1);
NIKE_coherent_z_save       = nan(nt,nz);
NIKE_broadband_z_save      = nan(nt,nz);
rho_save                   = nan(nt,nz);
N2_save                    = nan(nt,nz-1);
Ri_save                    = nan(nt,nz-1);
S2_save                    = nan(nt,nz-1);
u_surf_save                = nan(nt,length(twin_hours));
v_surf_save                = nan(nt,length(twin_hours));

WindPower_Wm2 = nan(nt,1);
tau_x_save    = nan(nt,1);
tau_y_save    = nan(nt,1);
u_sfc_save    = nan(nt,1);
v_sfc_save    = nan(nt,1);

NIKE_coherent_right_Jm2   = nan(nt,1);
NIKE_coherent_left_Jm2    = nan(nt,1);
NIKE_broadband_right_Jm2  = nan(nt,1);
NIKE_broadband_left_Jm2   = nan(nt,1);
NIKE_coherent_z_right     = nan(nt,nz);
NIKE_coherent_z_left      = nan(nt,nz);
NIKE_broadband_z_right    = nan(nt,nz);
NIKE_broadband_z_left     = nan(nt,nz);
Asymmetry_Ratio           = nan(nt,1);
WindStress_right_Nm2      = nan(nt,1);
WindStress_left_Nm2       = nan(nt,1);
WindPower_right_Wm2       = nan(nt,1);
WindPower_left_Wm2        = nan(nt,1);
heading_deg_save          = nan(nt,1);
R_NI_used_km              = nan(nt,1);

% Main loop
for it = 1:nt

    ix = find(lon_g>=lon_tc(it)-2 & lon_g<=lon_tc(it)+2);
    iy = find(lat_g>=lat_tc(it)-2 & lat_g<=lat_tc(it)+2);
    if isempty(ix)||isempty(iy), continue; end

    [LonSub,LatSub] = meshgrid(lon_g(ix),lat_g(iy));
    dist = distdim(distance(lat_tc(it),lon_tc(it),LatSub,LonSub),'deg','km');
    mask = double(dist<=R_NI); mask(mask==0)=NaN;
    R_NI_used_km(it) = R_NI;

    % Right/left half-plane masks: cross_z = dx*sin(h) - dy*cos(h)
    % Storm moving north (h=+pi/2): a point to the east has cross_z>0 = RIGHT
    h_rad = heading_rad(it);
    heading_deg_save(it) = h_rad * 180/pi;

    dy_pt = (LatSub - lat_tc(it)) * 111;
    dx_pt = (LonSub - lon_tc(it)) * 111 * cosd(lat_tc(it));
    cross_z = dx_pt * sin(h_rad) - dy_pt * cos(h_rad);

    mask_right = mask;  mask_right(cross_z <  0) = NaN;
    mask_left  = mask;  mask_left( cross_z >= 0) = NaN;

    distE  = distdim(distance(lat_tc(it),lon_tc(it),LatE,LonE),'deg','km');
    dy_E   = (LatE - lat_tc(it)) * 111;
    dx_E   = (LonE - lon_tc(it)) * 111 * cosd(lat_tc(it));
    cross_E = dx_E * sin(h_rad) - dy_E * cos(h_rad);
    mask_e        = double(distE <= R_NI); mask_e(mask_e==0) = NaN;
    mask_e_right  = mask_e; mask_e_right(cross_E <  0) = NaN;
    mask_e_left   = mask_e; mask_e_left( cross_E >= 0) = NaN;

    f_raw = 2*Omega*sind(lat_tc(it));
    if f_raw == 0
        f_val = 1e-5;
    else
        f_val = sign(f_raw) * max(abs(f_raw), 1e-5);
    end
    f_tc(it)        = f_val;
    Inertial_hr(it) = 2*pi/abs(f_val)/3600;

    u_ts        = nan(length(twin_hours),nz);
    v_ts        = nan(length(twin_hours),nz);
    u_ts_right  = nan(length(twin_hours),nz);
    v_ts_right  = nan(length(twin_hours),nz);
    u_ts_left   = nan(length(twin_hours),nz);
    v_ts_left   = nan(length(twin_hours),nz);
    tau_ts      = nan(length(twin_hours),2);

    avail_jt = [];
    for jt = 1:length(twin_hours)
        tgt = time_tc(it)+hours(twin_hours(jt));
        fn  = fullfile(hycom_path, sprintf('%s_DATA_%s.nc', hycom_prefix, datestr(tgt,'yyyymmdd_HH')));
        if isfile(fn), avail_jt(end+1) = jt; end %#ok<AGROW>
    end
    if isempty(avail_jt), continue; end
    [~, best_idx] = min(abs(twin_hours(avail_jt)));
    best_jt_idx   = avail_jt(best_idx);

    for jt = 1:length(twin_hours)
        tgt = time_tc(it)+hours(twin_hours(jt));
        fn  = fullfile(hycom_path, sprintf('%s_DATA_%s.nc', hycom_prefix, datestr(tgt,'yyyymmdd_HH')));
        if ~isfile(fn), continue; end

        try
            u_box = ncread(fn,'water_u',[ix(1) iy(1) 1],[length(ix) length(iy) nz]);
            v_box = ncread(fn,'water_v',[ix(1) iy(1) 1],[length(ix) length(iy) nz]);
        catch
            u_box = ncread(fn,'water_u',[ix(1) iy(1) 1 1],[length(ix) length(iy) nz 1]);
            v_box = ncread(fn,'water_v',[ix(1) iy(1) 1 1],[length(ix) length(iy) nz 1]);
        end

        for k=1:nz
            sliceU = u_box(:,:,k); sliceV = v_box(:,:,k);
            if ~isequal(size(sliceU),size(mask)), sliceU=sliceU'; sliceV=sliceV'; end

            u_ts(jt,k) = mean(sliceU.*mask,'all','omitnan');
            v_ts(jt,k) = mean(sliceV.*mask,'all','omitnan');

            u_ts_right(jt,k) = mean(sliceU.*mask_right,'all','omitnan');
            v_ts_right(jt,k) = mean(sliceV.*mask_right,'all','omitnan');
            u_ts_left(jt,k)  = mean(sliceU.*mask_left, 'all','omitnan');
            v_ts_left(jt,k)  = mean(sliceV.*mask_left, 'all','omitnan');
        end

        [~,ie_idx] = min(abs(t_e-tgt));
        tx_s = tau_x_all(:,:,ie_idx)'; ty_s = tau_y_all(:,:,ie_idx)';
        tau_ts(jt,1) = mean(tx_s .* mask_e, 'all', 'omitnan');
        tau_ts(jt,2) = mean(ty_s .* mask_e, 'all', 'omitnan');

        if jt == best_jt_idx
            % SST
            [~,ie_sst_idx] = min(abs(t_sst-tgt));
            if ndims(sst_all) == 3
                SST_K = interp2(LonS, LatS, sst_all(:,:,ie_sst_idx)', lon_tc(it), lat_tc(it));
                SST_C(it) = SST_K - 273.15;
            else
                SST_C(it) = interp2(LonS, LatS, sst_all', lon_tc(it), lat_tc(it));
            end

            % OHC
            [~,ie_ohc_idx] = min(abs(t_ohc-tgt));
            if ndims(ohc_all) == 3
                OHC_MJm2 = interp2(LonO, LatO, ohc_all(:,:,ie_ohc_idx)', lon_tc(it), lat_tc(it));
                OHC_Jm2(it) = OHC_MJm2 * 1e6;
            else
                OHC_Jm2(it) = interp2(LonO, LatO, ohc_all', lon_tc(it), lat_tc(it));
            end

            % Stokes drift
            [~,ie_st_idx] = min(abs(t_st-tgt));
            lon_query = lon_tc(it);
            if min(lon_st) >= 0
                lon_query = mod(lon_query, 360);
            else
                lon_query = mod(lon_query + 180, 360) - 180;
            end
            if ndims(ust_all) == 3
                us_val = interp2(LonSt, LatSt, ust_all(:,:,ie_st_idx)', lon_query, lat_tc(it));
                vs_val = interp2(LonSt, LatSt, vst_all(:,:,ie_st_idx)', lon_query, lat_tc(it));
            else
                us_val = interp2(LonSt, LatSt, ust_all', lon_query, lat_tc(it));
                vs_val = interp2(LonSt, LatSt, vst_all', lon_query, lat_tc(it));
            end
            Stokes_ms(it) = hypot(us_val, vs_val);
            if isnan(Stokes_ms(it))
                if ndims(ust_all) == 3
                    us_val = interp2(LonSt, LatSt, ust_all(:,:,ie_st_idx)', lon_query, lat_tc(it), 'nearest');
                    vs_val = interp2(LonSt, LatSt, vst_all(:,:,ie_st_idx)', lon_query, lat_tc(it), 'nearest');
                else
                    us_val = interp2(LonSt, LatSt, ust_all', lon_query, lat_tc(it), 'nearest');
                    vs_val = interp2(LonSt, LatSt, vst_all', lon_query, lat_tc(it), 'nearest');
                end
                Stokes_ms(it) = hypot(us_val, vs_val);
            end

            % Thermodynamics + GSW
            t_box = ncread(fn,'water_temp',[ix(1) iy(1) 1],[length(ix) length(iy) nz]);
            s_box = ncread(fn,'salinity'  ,[ix(1) iy(1) 1],[length(ix) length(iy) nz]);
            T_p   = squeeze(mean(t_box, [1,2], 'omitnan'));
            S_p   = squeeze(mean(s_box, [1,2], 'omitnan'));
            if size(T_p,1)==1, T_p=T_p'; S_p=S_p'; end

            P  = gsw_p_from_z(-z,lat_tc(it));
            SA = gsw_SA_from_SP(S_p,P,lon_tc(it),lat_tc(it));
            CT = gsw_CT_from_pt(SA,T_p);
            rho = gsw_rho(SA,CT,P);
            rho_save(it,:) = rho;
            [N2,~] = gsw_Nsquared(SA,CT,P,lat_tc(it));
            N2_save(it,:) = max(N2',1e-9);

            % MLD via dT = 0.2 C from shallowest valid level (de Boyer Montegut 2004)
            i_ref_mld = find(isfinite(T_p) & isfinite(rho), 1, 'first');
            if isempty(i_ref_mld)
                MLD_m(it) = NaN;
            else
                T_ref_mld = T_p(i_ref_mld);
                idx_thr = find(z > z(i_ref_mld) & T_p < T_ref_mld - 0.2, 1, 'first');
                if isempty(idx_thr)
                    valid_TS = find(isfinite(T_p) & isfinite(rho), 1, 'last');
                    MLD_m(it) = z(valid_TS);
                else
                    MLD_m(it) = z(idx_thr);
                end
            end

            [~, idx_mld_b] = min(abs(z - MLD_m(it)));
            ML_Buoyancy(it) = g*(rho(idx_mld_b) - rho(1)) / rho0;

            % Effective Coriolis (Doppler)
            u_surf = u_box(:,:,1); v_surf = v_box(:,:,1);
            if ~isequal(size(u_surf),size(mask)), u_surf=u_surf'; v_surf=v_surf'; end
            dy = 111e3 * mean(diff(lat_g(iy)));
            dx = 111e3 * cosd(lat_tc(it)) * mean(diff(lon_g(ix)));
            [du_dy,~] = gradient(u_surf,dy,dx);
            [~,dv_dx] = gradient(v_surf,dy,dx);
            zeta = mean(dv_dx-du_dy,'all','omitnan');
            f_eff(it) = f_val+0.5*zeta;
            WindStress_Nm2(it) = hypot(tau_ts(jt,1),tau_ts(jt,2));

            tx_R = mean(tx_s .* mask_e_right, 'all', 'omitnan');
            ty_R = mean(ty_s .* mask_e_right, 'all', 'omitnan');
            tx_L = mean(tx_s .* mask_e_left,  'all', 'omitnan');
            ty_L = mean(ty_s .* mask_e_left,  'all', 'omitnan');
            WindStress_right_Nm2(it) = hypot(tx_R, ty_R);
            WindStress_left_Nm2(it)  = hypot(tx_L, ty_L);

            ustar = sqrt(WindStress_Nm2(it)/rho0);
            if ~isnan(Stokes_ms(it)) && ~isnan(ustar)
                La_t(it) = sqrt(ustar / (Stokes_ms(it) + eps));
            end
        end
    end

    if all(isnan(u_ts(:))), continue; end

    arrays_to_fill = {'u_ts','v_ts','u_ts_right','v_ts_right','u_ts_left','v_ts_left'};
    for k = 1:nz
        for a = 1:length(arrays_to_fill)
            arr = eval(arrays_to_fill{a});
            good = isfinite(arr(:,k));
            if sum(good) >= 2
                arr(:,k) = fillmissing(arr(:,k),'linear','SamplePoints',twin_hours);
                eval([arrays_to_fill{a}, '(:,k) = arr(:,k);']);
            end
        end
    end

    u_surf_save(it,:) = u_ts(:,1)';
    v_surf_save(it,:) = v_ts(:,1)';

    % Single-frequency harmonic fit at f, depth-resolved
    [u_ni,       v_ni      ] = harmonic_fit_NI(u_ts,       v_ts,       f_val, tsec, nz);
    [u_ni_right, v_ni_right] = harmonic_fit_NI(u_ts_right, v_ts_right, f_val, tsec, nz);
    [u_ni_left,  v_ni_left ] = harmonic_fit_NI(u_ts_left,  v_ts_left,  f_val, tsec, nz);

    for k=1:nz
        NIKE_coherent_z_save(it,k)  = 0.5 * rho0 * (u_ni(k)^2       + v_ni(k)^2);
        NIKE_coherent_z_right(it,k) = 0.5 * rho0 * (u_ni_right(k)^2 + v_ni_right(k)^2);
        NIKE_coherent_z_left(it,k)  = 0.5 * rho0 * (u_ni_left(k)^2  + v_ni_left(k)^2);

        NIKE_broadband_z_save(it,k)  = 0.5 * rho_save(it,k) * mean(u_ts(:,k).^2       + v_ts(:,k).^2,       'omitnan');
        NIKE_broadband_z_right(it,k) = 0.5 * rho_save(it,k) * mean(u_ts_right(:,k).^2 + v_ts_right(:,k).^2, 'omitnan');
        NIKE_broadband_z_left(it,k)  = 0.5 * rho_save(it,k) * mean(u_ts_left(:,k).^2  + v_ts_left(:,k).^2,  'omitnan');
    end

    z_col = z(:);
    NIKE_coherent_Jm2(it)        = safe_trapz(z_col, NIKE_coherent_z_save(it,:)');
    NIKE_broadband_Jm2(it)       = safe_trapz(z_col, NIKE_broadband_z_save(it,:)');
    NIKE_coherent_right_Jm2(it)  = safe_trapz(z_col, NIKE_coherent_z_right(it,:)');
    NIKE_coherent_left_Jm2(it)   = safe_trapz(z_col, NIKE_coherent_z_left(it,:)');
    NIKE_broadband_right_Jm2(it) = safe_trapz(z_col, NIKE_broadband_z_right(it,:)');
    NIKE_broadband_left_Jm2(it)  = safe_trapz(z_col, NIKE_broadband_z_left(it,:)');

    % Asymmetry ratio with bilateral 10% guard + cap at 10
    if isfinite(NIKE_coherent_right_Jm2(it)) && ...
       isfinite(NIKE_coherent_left_Jm2(it))  && ...
       NIKE_coherent_left_Jm2(it)  > 0.10 * NIKE_coherent_Jm2(it) && ...
       NIKE_coherent_right_Jm2(it) > 0.10 * NIKE_coherent_Jm2(it)
        ar_raw = NIKE_coherent_right_Jm2(it) / NIKE_coherent_left_Jm2(it);
        Asymmetry_Ratio(it) = min(ar_raw, 10);
    end

    % ML / deep split
    nike_coher_prof = NIKE_coherent_z_save(it,:)';
    v_coh = isfinite(nike_coher_prof);
    idx_ml = (z_col <= MLD_m(it)) & v_coh;
    if sum(idx_ml) >= 2
        NIKE_above_MLD_Jm2(it) = trapz(z_col(idx_ml), nike_coher_prof(idx_ml));
    else
        NIKE_above_MLD_Jm2(it) = NaN;
    end
    NIKE_below_MLD_Jm2(it) = NIKE_coherent_Jm2(it) - NIKE_above_MLD_Jm2(it);

    % Energy-weighted centroid and upper-50m fraction
    prof_coh = NIKE_coherent_z_save(it,:)';
    prof_brd = NIKE_broadband_z_save(it,:)';
    v_c      = isfinite(prof_coh);
    v_brd    = isfinite(prof_brd);
    if sum(v_c)>=2
        NIKE_centroid_m(it) = trapz(z_col(v_c), z_col(v_c) .* prof_coh(v_c)) / (NIKE_coherent_Jm2(it) + eps);
    end
    idx50 = (z_col<=50);
    v50   = v_brd(idx50);
    if sum(v50)>=2
        NIKE_surface_frac(it) = trapz(z_col(idx50 & v_brd), prof_brd(idx50 & v_brd)) / (NIKE_broadband_Jm2(it) + eps);
    end

    % WKB-normalised energy and PE/KE
    N_prof   = sqrt(max([N2_save(it,:), N2_save(it,end)], 1e-9))';
    wkb_prof = (prof_brd .* (abs(f_val)./(N_prof+eps)));
    pe_prof  = (prof_brd .* (f_val^2 ./ (N_prof.^2+eps)));
    NIKE_WKB_Jm2(it) = safe_trapz(z_col, wkb_prof);
    NIW_PE_Jm2(it)   = safe_trapz(z_col, pe_prof);
    PE_to_KE_Ratio(it) = NIW_PE_Jm2(it) / (NIKE_broadband_Jm2(it) + eps);

    % Richardson number: S2 floored at 1e-7 (resolution noise floor),
    % mask cells with weak stratification or extreme Ri (numerical artifacts)
    idx0 = find(twin_hours==0);
    u0 = u_ts(idx0,:)'; v0 = v_ts(idx0,:)';
    du = diff(u0)./diff(z);
    dv = diff(v0)./diff(z);

    S2_raw = du.^2 + dv.^2;
    S2     = max(S2_raw, 1e-7);
    S2_save(it,:) = S2';

    N2_row = N2_save(it,:)';
    Ri_row = N2_row ./ S2;

    bad_mask = (N2_row < 1e-7) | (S2_raw < 1e-8) | ...
               ~isfinite(Ri_row) | (Ri_row < 1e-3) | (Ri_row > 1e4);
    Ri_row(bad_mask) = NaN;
    Ri_save(it,:) = Ri_row;

    [~,idx_mld] = min(abs(z-MLD_m(it)));
    idx_mld     = min(nz-1, max(1, idx_mld));
    if isfinite(Ri_row(idx_mld))
        Ri_Bulk_MLD(it) = Ri_row(idx_mld);
    else
        Ri_Bulk_MLD(it) = NaN;
    end

    if it > 1
        Centroid_Trend_m_day(it) = (NIKE_centroid_m(it) - NIKE_centroid_m(it-1)) / days(time_tc(it)-time_tc(it-1));
        if abs(Centroid_Trend_m_day(it)) > 240
            Centroid_Trend_m_day(it) = NaN;
        end
        Entrainment_Velocity_m_day(it) = (MLD_m(it) - MLD_m(it-1)) / days(time_tc(it)-time_tc(it-1));
    end

    % WKB vertical group velocity: cgz = (omega^2 - f^2)/(omega*m), omega = 1.10*f
    nike_prof_cgz = NIKE_coherent_z_save(it,:)';
    valid_pc = isfinite(nike_prof_cgz) & nike_prof_cgz > 0;
    if sum(valid_pc) >= 4 && ~isnan(MLD_m(it))
        idx_below = z > MLD_m(it) & valid_pc;
        if sum(idx_below) >= 4
            [~, iz_local] = max(nike_prof_cgz(idx_below));
            z_below = z(idx_below);
            z_at_peak = z_below(iz_local);

            quarter_lam = max(z_at_peak - MLD_m(it), 30);
            lam_z = max(4 * quarter_lam, 80);
            lam_z = min(lam_z, 800);
            m_dom = 2*pi / lam_z;

            iz_n_lo = find(z >= z_at_peak - 50, 1, 'first');
            iz_n_hi = find(z <= z_at_peak + 50, 1, 'last');
            if isempty(iz_n_lo) || isempty(iz_n_hi) || iz_n_lo > iz_n_hi
                N_at_peak = sqrt(max(N2_save(it, min(iz_local, size(N2_save,2))), 1e-9));
            else
                iz_n_hi = min(iz_n_hi, size(N2_save,2));
                N_at_peak = sqrt(max(mean(N2_save(it, iz_n_lo:iz_n_hi),'omitnan'), 1e-9));
            end

            omega_ni = abs(f_val) * 1.10;
            cgz_wkb_ms = abs((omega_ni^2 - f_val^2) / (omega_ni * m_dom + eps));
            cgz_candidate = cgz_wkb_ms * 86400;

            if cgz_candidate >= 1 && cgz_candidate <= 240
                Inferred_Cgz_m_day(it) = cgz_candidate;
            end
        end
    end

    Coherence_Ratio(it) = NIKE_coherent_Jm2(it) / (NIKE_broadband_Jm2(it) + eps);

    % Wind power (instantaneous, t=0 snapshot)
    [~, ie_t0] = min(abs(t_e - time_tc(it)));
    tx_t0 = tau_x_all(:,:,ie_t0)';
    ty_t0 = tau_y_all(:,:,ie_t0)';

    distE_t0  = distdim(distance(lat_tc(it),lon_tc(it),LatE,LonE),'deg','km');
    dy_E0     = (LatE - lat_tc(it)) * 111;
    dx_E0     = (LonE - lon_tc(it)) * 111 * cosd(lat_tc(it));
    cross_E0  = dx_E0 * sin(h_rad) - dy_E0 * cos(h_rad);
    mask_e_t0       = double(distE_t0 <= R_NI); mask_e_t0(mask_e_t0 == 0) = NaN;
    mask_e_t0_right = mask_e_t0; mask_e_t0_right(cross_E0 <  0) = NaN;
    mask_e_t0_left  = mask_e_t0; mask_e_t0_left( cross_E0 >= 0) = NaN;

    tau_x0 = mean(tx_t0 .* mask_e_t0, 'all', 'omitnan');
    tau_y0 = mean(ty_t0 .* mask_e_t0, 'all', 'omitnan');

    idx0_wp = find(twin_hours == 0);
    if ~isempty(idx0_wp) && isfinite(u_ts(idx0_wp,1)) && ...
       isfinite(tau_x0) && isfinite(tau_y0)
        u_sfc_save(it)    = u_ts(idx0_wp, 1);
        v_sfc_save(it)    = v_ts(idx0_wp, 1);
        tau_x_save(it)    = tau_x0;
        tau_y_save(it)    = tau_y0;
        WindPower_Wm2(it) = tau_x0 * u_ts(idx0_wp,1) + tau_y0 * v_ts(idx0_wp,1);

        tx_R0 = mean(tx_t0 .* mask_e_t0_right, 'all', 'omitnan');
        ty_R0 = mean(ty_t0 .* mask_e_t0_right, 'all', 'omitnan');
        tx_L0 = mean(tx_t0 .* mask_e_t0_left,  'all', 'omitnan');
        ty_L0 = mean(ty_t0 .* mask_e_t0_left,  'all', 'omitnan');
        if isfinite(u_ts_right(idx0_wp,1)) && isfinite(tx_R0)
            WindPower_right_Wm2(it) = tx_R0 * u_ts_right(idx0_wp,1) + ty_R0 * v_ts_right(idx0_wp,1);
        end
        if isfinite(u_ts_left(idx0_wp,1)) && isfinite(tx_L0)
            WindPower_left_Wm2(it)  = tx_L0 * u_ts_left(idx0_wp,1)  + ty_L0 * v_ts_left(idx0_wp,1);
        end
    else
        WindPower_Wm2(it) = NaN;
        tau_x_save(it)    = tau_x0;
        tau_y_save(it)    = tau_y0;
        u_sfc_save(it)    = NaN;
        v_sfc_save(it)    = NaN;
    end

end

% Budget and cumulative wind work
v_nike    = NIKE_above_MLD_Jm2(isfinite(NIKE_above_MLD_Jm2));
nike_gain = max(v_nike) - prctile(v_nike,20);

t_sec = seconds(time_tc - time_tc(1));
if any(diff(t_sec) <= 0)
    warning('Duplicate or non-monotonic timestamps detected.');
    bad = find(diff(t_sec) <= 0) + 1;
    t_sec(bad) = t_sec(bad-1) + eps;
end

dNIKE_dt_Wm2        = gradient(NIKE_broadband_Jm2, t_sec);
dNIW_PE_dt_Wm2      = gradient(NIW_PE_Jm2,         t_sec);
Eff_Energy_Loss_Wm2 = -abs(dNIKE_dt_Wm2);

% Cumulative wind work with 3-day pre-storm offset (0.5x mean P)
WindPower_pos = WindPower_Wm2;
WindPower_pos(WindPower_pos < 0 | ~isfinite(WindPower_pos)) = 0;

pre_storm_hours  = 72;
pp_valid = WindPower_pos > 0;
if any(pp_valid)
    pre_storm_mean_P = 0.5 * mean(WindPower_pos(pp_valid));
else
    pre_storm_mean_P = 0;
end
pre_storm_offset = pre_storm_mean_P * pre_storm_hours * 3600;

WindWork_cum_Jm2 = cumtrapz(t_sec, WindPower_pos) + pre_storm_offset;

% Build output table
Tmain = table(time_tc, lat_tc, lon_tc, ...
    NIKE_coherent_Jm2, NIKE_broadband_Jm2, dNIKE_dt_Wm2, ...
    NIKE_surface_frac, NIKE_centroid_m, ...
    WindStress_Nm2, WindPower_Wm2, WindWork_cum_Jm2, ...
    SST_C, OHC_Jm2, MLD_m, Stokes_ms, La_t, ML_Buoyancy, Inertial_hr, ...
    Ri_Bulk_MLD, Entrainment_Velocity_m_day, Coherence_Ratio, NIKE_WKB_Jm2, ...
    NIKE_above_MLD_Jm2, NIKE_below_MLD_Jm2, NIW_PE_Jm2, PE_to_KE_Ratio, ...
    Inferred_Cgz_m_day, Centroid_Trend_m_day, Eff_Energy_Loss_Wm2, dNIW_PE_dt_Wm2, ...
    repmat(nike_gain/1000, nt, 1), ...
    NIKE_coherent_right_Jm2, NIKE_coherent_left_Jm2, ...
    NIKE_broadband_right_Jm2, NIKE_broadband_left_Jm2, ...
    Asymmetry_Ratio, ...
    WindStress_right_Nm2, WindStress_left_Nm2, ...
    WindPower_right_Wm2,  WindPower_left_Wm2, ...
    heading_deg_save, R_NI_used_km);

Tmain.Properties.VariableNames = { ...
    'time','lat','lon', 'NIKE_Coherent','NIKE_Broadband','dNIKE_dt', ...
    'SurfFrac','Centroid', 'Stress', 'WindPower_Wm2','WindWork_cum_Jm2', ...
    'SST','OHC','MLD','Stokes','La_t','Buoyancy','Inertial_hr', ...
    'Ri_Bulk','Entrainment','Coherence_Ratio','NIKE_WKB', 'NIKE_ML','NIKE_Deep', ...
    'NIW_PE','PE_KE_Ratio', 'Cgz_Proxy','Centroid_Trend_m_day','EnergyLoss','dPE_dt', ...
    'NetGain_kJ', ...
    'NIKE_Coherent_right','NIKE_Coherent_left', ...
    'NIKE_Broadband_right','NIKE_Broadband_left', ...
    'Asymmetry_Ratio', ...
    'WindStress_right','WindStress_left', ...
    'WindPower_right','WindPower_left', ...
    'TC_heading_deg','R_NI_used_km'};

save(out_mat,'z','NIKE_coherent_z_save','NIKE_broadband_z_save', ...
    'NIKE_coherent_z_right','NIKE_coherent_z_left', ...
    'NIKE_broadband_z_right','NIKE_broadband_z_left', ...
    'N2_save','Ri_save','S2_save','rho_save','u_surf_save','v_surf_save', ...
    'f_tc','f_eff','SST_C', ...
    'OHC_Jm2','Stokes_ms','La_t','MLD_m','NIKE_broadband_Jm2', ...
    'NIKE_coherent_Jm2','WindStress_Nm2', ...
    'WindPower_Wm2','WindWork_cum_Jm2','tau_x_save','tau_y_save','u_sfc_save','v_sfc_save', ...
    'dNIKE_dt_Wm2','NIKE_surface_frac','NIKE_centroid_m','ML_Buoyancy', ...
    'Inertial_hr','Ri_Bulk_MLD','Entrainment_Velocity_m_day','Coherence_Ratio', ...
    'NIKE_WKB_Jm2','NIKE_above_MLD_Jm2','NIKE_below_MLD_Jm2','NIW_PE_Jm2', ...
    'PE_to_KE_Ratio','Inferred_Cgz_m_day','Centroid_Trend_m_day', ...
    'Eff_Energy_Loss_Wm2','dNIW_PE_dt_Wm2', ...
    'NIKE_coherent_right_Jm2','NIKE_coherent_left_Jm2', ...
    'NIKE_broadband_right_Jm2','NIKE_broadband_left_Jm2', ...
    'Asymmetry_Ratio','WindStress_right_Nm2','WindStress_left_Nm2', ...
    'WindPower_right_Wm2','WindPower_left_Wm2', ...
    'heading_deg_save','R_NI_used_km', ...
    '-v7.3');

writetable(Tmain,out_xlsx);

fprintf('\nML / Total NIKE             : %.2f %%\n',100*mean(NIKE_above_MLD_Jm2./(NIKE_coherent_Jm2+eps),'omitnan'));
fprintf('Below-ML NIKE fraction      : %.2f %%\n',100*mean(1 - (NIKE_above_MLD_Jm2./(NIKE_coherent_Jm2+eps)),'omitnan'));
fprintf('NIKE centroid depth         : %.1f m\n',mean(NIKE_centroid_m,'omitnan'));
fprintf('Wind Power mean / max       : %.4f / %.4f W/m2\n', mean(WindPower_Wm2,'omitnan'), max(WindPower_Wm2));
fprintf('Wind Stress mean            : %.4f N/m2\n', mean(WindStress_Nm2,'omitnan'));
fprintf('Cumulative WindWork (max)   : %.2f kJ/m2\n', max(WindWork_cum_Jm2)/1000);
fprintf('Mean Cgz (WKB)              : %.1f m/day\n', mean(Inferred_Cgz_m_day,'omitnan'));
fprintf('Mean / median R/L           : %.2f / %.2f\n', ...
    mean(Asymmetry_Ratio,'omitnan'), median(Asymmetry_Ratio,'omitnan'));
fprintf('=== %s : DONE ===\n', cfg.name);


% Summary report
storm_tag = upper(extractBefore(out_xlsx, '_V16'));
summary_file = sprintf('%s_SUMMARY_v6.txt', storm_tag);
fid = fopen(summary_file, 'w');
write_line = @(varargin) [fprintf(varargin{:}), fprintf(fid, varargin{:})];

write_line('\n=== %s RESULTS SUMMARY ===\n', storm_tag);
write_line('Track timesteps    : %d\n', nt);
write_line('Mean lat / lon     : %.2f N / %.2f E\n', mean(lat_tc,'omitnan'), mean(lon_tc,'omitnan'));
write_line('R_NI used          : %.0f km\n', R_NI);
write_line('Heading smoothing  : %d-pt\n', heading_window);
write_line('Time window        : %d to %d hours\n', twin_hours(1), twin_hours(end));

write_line('\n--- ENERGETICS ---\n');
write_line('Peak NIKE coherent : %.2f kJ/m2\n', max(NIKE_coherent_Jm2)/1000);
write_line('Peak NIKE broadband: %.2f kJ/m2\n', max(NIKE_broadband_Jm2)/1000);
write_line('Mean NIKE_ML       : %.2f kJ/m2\n', mean(NIKE_above_MLD_Jm2,'omitnan')/1000);
write_line('Mean NIKE_deep     : %.2f kJ/m2\n', mean(NIKE_below_MLD_Jm2,'omitnan')/1000);
write_line('Below-ML fraction  : %.1f %%\n', ...
    100*mean(NIKE_below_MLD_Jm2./(NIKE_coherent_Jm2+eps),'omitnan'));
write_line('Mean MLD           : %.1f m\n', mean(MLD_m,'omitnan'));
write_line('NIKE centroid mean : %.1f m\n', mean(NIKE_centroid_m,'omitnan'));
write_line('Coherence ratio    : %.3f mean\n', mean(Coherence_Ratio,'omitnan'));

write_line('\n--- WIND FORCING ---\n');
write_line('Mean / max WindStress : %.4f / %.4f N/m2\n', mean(WindStress_Nm2,'omitnan'), max(WindStress_Nm2));
write_line('Mean / max WindPower  : %.4f / %.4f W/m2\n', mean(WindPower_Wm2,'omitnan'), max(WindPower_Wm2));
write_line('Pre-storm offset      : %.2f kJ/m2\n', pre_storm_offset/1000);
write_line('Cumulative WindWork   : %.2f kJ/m2\n', max(WindWork_cum_Jm2)/1000);

write_line('\n--- STRATIFICATION & MIXING ---\n');
write_line('MLD range  : %.1f to %.1f m\n', min(MLD_m), max(MLD_m));

% Ri summary on upper 100 m, using 5th percentile (Jampana 2018)
Ri_upper100 = Ri_save(:, z(1:size(Ri_save,2)) <= 100);
Ri_valid    = Ri_upper100(isfinite(Ri_upper100));
n_valid     = numel(Ri_valid);
if n_valid > 0
    Ri_p05  = prctile(Ri_valid, 5);
    Ri_p25  = prctile(Ri_valid, 25);
    Ri_med  = median(Ri_valid);
    pct_KH  = 100 * sum(Ri_valid < 0.25) / n_valid;
else
    Ri_p05 = NaN; Ri_p25 = NaN; Ri_med = NaN; pct_KH = NaN;
end
write_line('5th %%ile Ri (upper100): %.3f\n', Ri_p05);
write_line('25th %%ile Ri (upper)  : %.3f\n', Ri_p25);
write_line('Median Ri (upper100)  : %.3f\n', Ri_med);
write_line('Sub-critical %%        : %.1f %%\n', pct_KH);
write_line('N valid Ri cells      : %d / %d\n', n_valid, numel(Ri_upper100));
write_line('Mean / median Cgz     : %.1f / %.1f m/day\n', ...
    mean(Inferred_Cgz_m_day,'omitnan'), median(Inferred_Cgz_m_day,'omitnan'));
write_line('Centroid trend (mean) : %.1f m/day\n', mean(Centroid_Trend_m_day,'omitnan'));
write_line('Mean PE/KE ratio      : %.3f\n', mean(PE_to_KE_Ratio,'omitnan'));

write_line('\n--- LANGMUIR ---\n');
La_valid = La_t(isfinite(La_t));
if ~isempty(La_valid)
    write_line('Mean La_t          : %.3f\n', mean(La_valid));
    write_line('Strong LC fraction : %.1f %%\n', 100*sum(La_valid<0.35)/length(La_valid));
    write_line('Weak LC fraction   : %.1f %%\n', 100*sum(La_valid>0.7)/length(La_valid));
end

write_line('\n--- ASYMMETRY ---\n');
AR_v = Asymmetry_Ratio(isfinite(Asymmetry_Ratio));
if ~isempty(AR_v)
    write_line('Mean R/L ratio     : %.2f\n', mean(AR_v));
    write_line('Median R/L ratio   : %.2f\n', median(AR_v));
    n_R = sum(AR_v >= 1.5);
    n_S = sum(AR_v >= 0.85 & AR_v < 1.18);
    n_L = sum(AR_v <= 0.67);
    n_X = sum( (AR_v > 0.67 & AR_v < 0.85) | (AR_v >= 1.18 & AR_v < 1.5) );
    write_line('Right-bias steps   : %.1f %%\n', 100*n_R/nt);
    write_line('Symmetric steps    : %.1f %%\n', 100*n_S/nt);
    write_line('Left-bias steps    : %.1f %%\n', 100*n_L/nt);
    write_line('Transitional steps : %.1f %%\n', 100*n_X/nt);
    write_line('Mean NIKE_right    : %.2f kJ/m2\n', mean(NIKE_coherent_right_Jm2,'omitnan')/1000);
    write_line('Mean NIKE_left     : %.2f kJ/m2\n', mean(NIKE_coherent_left_Jm2,'omitnan')/1000);
    write_line('Valid R/L steps    : %d / %d\n', length(AR_v), nt);
end

fclose(fid);
fprintf('Summary saved: %s\n', summary_file);


% --- Local helpers ---
function [u_ni, v_ni] = harmonic_fit_NI(u_ts, v_ts, f_val, tsec, nz)
% Fit u(t) = A*cos(f*t) + B*sin(f*t) at each depth.
% Amplitude capped at sqrt(2)*std to prevent leakage overshoot.
    u_ni = nan(nz,1);
    v_ni = nan(nz,1);
    for k = 1:nz
        g_t = isfinite(u_ts(:,k)) & isfinite(v_ts(:,k));
        if sum(g_t) < 11, continue; end
        X_fit = [cos(f_val*tsec(g_t)) sin(f_val*tsec(g_t))];
        cu = X_fit \ u_ts(g_t,k);
        cv = X_fit \ v_ts(g_t,k);

        amp_u = hypot(cu(1),cu(2));
        amp_v = hypot(cv(1),cv(2));
        std_u = std(u_ts(g_t,k), 'omitnan');
        std_v = std(v_ts(g_t,k), 'omitnan');

        u_ni(k) = min(amp_u, std_u * sqrt(2));
        v_ni(k) = min(amp_v, std_v * sqrt(2));
    end
end

function I = safe_trapz(z, prof)
% Trapezoidal integration ignoring NaNs; returns NaN if fewer than 2 valid points.
    v = isfinite(prof);
    if sum(v) >= 2
        I = trapz(z(v), prof(v));
    else
        I = NaN;
    end
end

function cfg = get_storm_config(storm_id)
% Per-storm configuration. To add a new storm, add a case block.
% Data is expected under ./data/ — adjust base_data if you keep data elsewhere.
    base_data   = fullfile(pwd, 'data');
    base_track  = fullfile(base_data, 'IMD',        filesep);
    base_hycom  = fullfile(base_data, 'HYCOM',      filesep);
    base_stress = fullfile(base_data, 'Stress',     filesep);
    base_wind   = fullfile(base_data, 'Wind',       filesep);
    base_curr   = fullfile(base_data, 'Currents',   filesep);
    base_ohc    = fullfile(base_data, 'INCOIS_OHC', filesep);

    switch upper(storm_id)
        case 'KYARR'
            cfg.name        = 'KYARR';
            cfg.year_ref    = 2019;
            cfg.heading_smooth_window = 3;
            cfg.hycom_prefix = 'KYARR';
            cfg.track_path  = [base_track  'Kyarr.csv'];
            cfg.hycom_path  = [base_hycom  'Kyarr' filesep];
            cfg.stress_file = [base_stress 'Kyarr' filesep 'era5_hourly_stress_KYARR_FINAL.nc'];
            cfg.wind_file   = [base_wind   'Kyarr' filesep 'hourly_wind_kyarr.nc'];
            cfg.stokes_file = [base_curr   'Kyarr' filesep 'ocean_stokes_kyarr.nc'];
            cfg.sst_file    = [base_curr   'Kyarr' filesep 'ocean_sst_kyarr.nc'];
            cfg.ohc_file    = [base_ohc    'Kyarr' filesep '2019_20_ohc700.nc'];
            cfg.out_mat     = 'KYARR_V16_SUPERCHARGED.mat';
            cfg.out_xlsx    = 'KYARR_V16_SUPERCHARGED.xlsx';

        case 'AMPHAN'
            cfg.name        = 'AMPHAN';
            cfg.year_ref    = 2020;
            cfg.heading_smooth_window = 3;
            cfg.hycom_prefix = 'AMPHAN';
            cfg.track_path  = [base_track  'Amphan.csv'];
            cfg.hycom_path  = [base_hycom  'Amphan' filesep];
            cfg.stress_file = [base_stress 'Amphan' filesep 'era5_hourly_stress_FINAL.nc'];
            cfg.wind_file   = [base_wind   'Amphan' filesep 'era5_10m_wind_May_2020.nc'];
            cfg.stokes_file = [base_curr   'Amphan' filesep 'ocean_stokes.nc'];
            cfg.sst_file    = [base_curr   'Amphan' filesep 'ocean_sst.nc'];
            cfg.ohc_file    = [base_ohc    'Amphan' filesep 'ohc700_May2020.nc'];
            cfg.out_mat     = 'AMPHAN_V16_SUPERCHARGED.mat';
            cfg.out_xlsx    = 'AMPHAN_V16_SUPERCHARGED.xlsx';

        case 'FANI'
            cfg.name        = 'FANI';
            cfg.year_ref    = 2019;
            cfg.heading_smooth_window = 5;   % recurves near landfall
            cfg.hycom_prefix = 'FANI';
            cfg.track_path  = [base_track  'Fani.csv'];
            cfg.hycom_path  = [base_hycom  'Fani' filesep];
            cfg.stress_file = [base_stress 'Fani' filesep 'era5_hourly_stress_FANI_FINAL.nc'];
            cfg.wind_file   = [base_wind   'Fani' filesep 'hourly_wind_fani.nc'];
            cfg.stokes_file = [base_curr   'Fani' filesep 'ocean_stokes_fani.nc'];
            cfg.sst_file    = [base_curr   'Fani' filesep 'ocean_sst_fani.nc'];
            cfg.ohc_file    = [base_ohc    'Fani' filesep '2019_20_ohc700.nc'];
            cfg.out_mat     = 'FANI_V16_SUPERCHARGED.mat';
            cfg.out_xlsx    = 'FANI_V16_SUPERCHARGED.xlsx';

        case 'TAUKTAE'
            cfg.name        = 'TAUKTAE';
            cfg.year_ref    = 2021;
            cfg.heading_smooth_window = 5;   % sharp recurve along W coast
            cfg.hycom_prefix = 'TAUKTAE';
            cfg.track_path  = [base_track  'Tauktae.csv'];
            cfg.hycom_path  = [base_hycom  'Tauktae' filesep];
            cfg.stress_file = [base_stress 'Tauktae' filesep 'era5_hourly_stress_TAUKTAE_FINAL.nc'];
            cfg.wind_file   = [base_wind   'tauktae' filesep 'hourly_wind_tauktae.nc'];
            cfg.stokes_file = [base_curr   'tauktae' filesep 'ocean_stokes_tauktae.nc'];
            cfg.sst_file    = [base_curr   'tauktae' filesep 'ocean_sst_tauktae.nc'];
            cfg.ohc_file    = [base_ohc    'tauktae' filesep '2021_22_ohc700.nc'];
            cfg.out_mat     = 'TAUKTAE_V16_SUPERCHARGED.mat';
            cfg.out_xlsx    = 'TAUKTAE_V16_SUPERCHARGED.xlsx';

        otherwise
            error('Unknown STORM_ID: %s. Valid: KYARR | AMPHAN | FANI | TAUKTAE', storm_id);
    end
end
