% Main-text figures 2-7, plus S10 and S11. Reads *_V17.mat and *_POST_V17.mat.

clc; clear; close all;

DATA_DIR = 'path/to/data';        % folder with HYCOM, IMD, Stress, ARGO etc.

%% USER CONFIG
BASE_DIR   = fullfile(DATA_DIR,'Final','All','Combined');
TRACK_DIR  = fullfile(DATA_DIR,'IMD');
VPDIR      = fullfile(DATA_DIR,'Final','All','Combined','tests');  % *_vertprop.mat
OUT_DIR    = fullfile(pwd, 'MAIN_FIGURES');

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

STORM_NAMES = {'KYARR','AMPHAN','FANI','TAUKTAE'};
STORM_BASIN = {'AS',   'BoB',   'BoB', 'AS'};
TRACK_CSVS  = {'Kyarr.csv','Amphan.csv','Fani.csv','Tauktae.csv'};

FIGURES_TO_MAKE = [2 3 4 5 6 7 8 9];

n_cyc = numel(STORM_NAMES);
Omega = 7.2921e-5;

%% GLOBAL STYLE
% [F12] The interactive toolbar was being captured inside the exported PNG
% (visible in the Fig 6b Fani panel). Turn the toolbar and menu bar off for
% every figure this script creates, before any figure is opened.
set(groot, ...
    'defaultFigureToolBar','none', ...
    'defaultFigureMenuBar','none', ...
    'defaultFigureColor','w', ...
    'defaultAxesFontName','Helvetica', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesLineWidth',1.1, ...
    'defaultLineLineWidth',1.6, ...
    'defaultAxesTitleFontWeight','bold');

COL.stress  = [0.20 0.40 0.70];
COL.wpos    = [0.10 0.55 0.30];
COL.wneg    = [0.85 0.20 0.20];
COL.work    = [0.45 0.20 0.55];
COL.work2   = [0.75 0.55 0.20];
COL.bb      = [0.05 0.32 0.58];
COL.coh     = [0.85 0.45 0.10];
COL.S       = [0.45 0.20 0.55];
COL.cent    = [0.00 0.00 0.00];
COL.right   = [0.85 0.20 0.20];
COL.left    = [0.20 0.45 0.85];
COL.grid    = [0.85 0.85 0.85];
COL.zoneSub = [0.85 0.92 1.00];
COL.zoneRes = [0.85 1.00 0.85];
COL.zoneSup = [1.00 0.88 0.85];
COL.litband = [0.85 1.00 0.85];
COL.cw      = [0.05 0.32 0.58];
COL.ccw     = [0.85 0.45 0.10];

beautify = @() set(gca,'TickDir','out','Box','on','Layer','top', ...
    'FontSize',10,'GridAlpha',0.22,'XGrid','on','YGrid','on', ...
    'GridColor',COL.grid);

savepro = @(name) print(gcf, fullfile(OUT_DIR,name), '-dpng','-r300');

LIT.RL_lo=1.5; LIT.RL_hi=3.0;
LIT.S_lo=0.7;  LIT.S_hi=1.3;
% [F6] Penetration literature band disabled. Under the revised definition the
% penetration fraction is normalised by the TOTAL windowed kinetic energy, not
% by a near-inertial quantity, so the 10-45% literature range (which is
% near-inertial over near-inertial) is no longer a like-for-like comparison.
% Set to NaN => no shaded band, neutral bar colour, no pass/fail claim.
LIT.pen_lo=NaN; LIT.pen_hi=NaN;
LIT.RC_thr=0.5;
LIT.lamz_lo=200; LIT.lamz_hi=400;   % vertical wavelength lit range (m)

fprintf('\n================================================================\n');
fprintf(' PAPER 1 - COMBINED MAIN FIGURES v7 (Fig 2..Fig 9)\n');
fprintf(' Output : %s\n', OUT_DIR);
fprintf('================================================================\n');

%% LOAD ALL DATA ONCE (v7)
D = struct();
for i = 1:n_cyc
    storm = STORM_NAMES{i};
    fprintf('Loading %s ...\n', storm);
    mat_main = fullfile(BASE_DIR, [storm, '_V17.mat']);
    mat_post = fullfile(BASE_DIR, [storm, '_POST_V17.mat']);
    vp_file  = fullfile(VPDIR,    [storm, '_vertprop.mat']);
    trk_csv  = fullfile(TRACK_DIR, TRACK_CSVS{i});
    if ~isfile(mat_main), error('Missing %s', mat_main); end
    M = load(mat_main);
    time = M.time_tc(:); nT = numel(time);

    D(i).name=storm; D(i).basin=STORM_BASIN{i}; D(i).time=time; D(i).M=M;
    D(i).lat_mean = mean(M.lat_tc,'omitnan');

    gv = @(f) gof(M,f,nT);
    D(i).stress   = gv('WindStress_Nm2');
    D(i).wp       = gv('WindPower_NI_Wm2');
    D(i).ww_signed= gv('WindWork_signed_Jm2')/1000;
    D(i).ww_pos   = gv('WindWork_pos_Jm2')/1000;
    D(i).nike_bb  = gv('NIKE_broadband_Jm2')/1000;
    D(i).nike_coh = gv('NIKE_coherent_Jm2')/1000;
    D(i).dnike_dt = gv('dNIKE_dt_Wm2');
    D(i).MLD      = gv('MLD_m');
    D(i).ILD      = gv('ILD_m');
    D(i).BLT      = gv('BLT_m');
    D(i).cent     = gv('NIKE_centroid_m');
    D(i).coh_ratio= gv('Coherence_Ratio');
    D(i).nike_ml  = gv('NIKE_above_MLD_Jm2')/1000;
    D(i).nike_dp  = gv('NIKE_below_MLD_Jm2')/1000;
    D(i).nike_R   = gv('NIKE_coherent_right_Jm2')/1000;
    D(i).nike_L   = gv('NIKE_coherent_left_Jm2')/1000;
    D(i).AR       = gv('Asymmetry_Ratio');
    D(i).stress_R = gv('WindStress_right_Nm2');
    D(i).stress_L = gv('WindStress_left_Nm2');
    D(i).phase    = gv('phase_tau_uI_deg');
    D(i).phase_sd = gv('phase_sd_deg');
    D(i).uI_amp   = gv('uI_amp_surf');
    D(i).zeta     = gv('zeta_save');
    D(i).strain   = gv('strain_rate');
    D(i).f_tc     = gv('f_tc');
    D(i).f_eff    = gv('f_eff');
    D(i).Ri_min   = gv('Ri_min_100m');
    D(i).La_t     = gv('La_t');
    D(i).SST      = gv('SST_C');
    D(i).OHC      = gv('OHC_Jm2')/1e7;
    D(i).peke     = gv('PE_to_KE');
    D(i).zeta_f   = abs(D(i).zeta)./abs(D(i).f_tc);

    if isfield(M,'z'), D(i).z=M.z(:); else, D(i).z=[]; end
    D(i).NIKE_coh_z = foe(M,'NIKE_coherent_z_save')/1000;
    D(i).NIKE_brd_z = foe(M,'NIKE_broadband_z_save')/1000;
    if ~isempty(D(i).NIKE_coh_z) && ~isempty(D(i).NIKE_brd_z)
        cz = D(i).NIKE_coh_z ./ (D(i).NIKE_brd_z + eps);
        cz(~isfinite(cz))=NaN; cz(cz<0)=NaN; D(i).Coh_z=cz;
    else, D(i).Coh_z=[]; end

    % surface velocity time series (rotary spectra)
    D(i).u_surf = []; D(i).v_surf = [];
    if isfield(M,'u_surf_save')
        i0 = ceil(size(M.u_surf_save,2)/2);   % window centre column
        D(i).u_surf = M.u_surf_save(:,i0);
        D(i).v_surf = M.v_surf_save(:,i0);
    end

    % vertical propagation file (Fig 5 phase tilt)
    D(i).has_vp=false;
    if isfile(vp_file)
        try
            VP = load(vp_file);
            D(i).vp = VP; D(i).has_vp = true;
        catch, end
    end

    % post file (penetration)
    D(i).pen=nan(nT,1); D(i).pen_R=nan(nT,1); D(i).pen_L=nan(nT,1);
    D(i).has_post=false;
    if isfile(mat_post)
        try
            P=load(mat_post); Tp=P.Tpost; D(i).Tp=Tp; D(i).stat=P.st;
            if ismember('NIW_Penetration_Frac',Tp.Properties.VariableNames)
                D(i).pen = Tp.NIW_Penetration_Frac; end
            if ismember('Penetration_right',Tp.Properties.VariableNames)
                D(i).pen_R=Tp.Penetration_right; D(i).pen_L=Tp.Penetration_left; end
            D(i).has_post=true;
        catch, end
    end

    % S parameter from track CSV
    D(i).S_t=nan(nT,1); D(i).S_median_trk=NaN;
    if isfile(trk_csv)
        try
            T=readtable(trk_csv);
            tt=datetime(T.year,T.month,T.day,T.hour,0,0);
            ft=2*Omega*sind(T.lat); ft=(sign(ft)+(ft==0)).*max(abs(ft),1e-5);
            if ismember('Speed',T.Properties.VariableNames), Uh=T.Speed*0.514444;
            else, Uh=nan(height(T),1); end
            if ismember('rmax',T.Properties.VariableNames), Rm=T.rmax*1852;
            elseif ismember('Rmax',T.Properties.VariableNames), Rm=T.Rmax*1852;
            else, Rm=repmat(25*1852,height(T),1); end
            St=movmean(Uh./(abs(ft).*Rm+eps),3,'omitnan');
            D(i).S_median_trk=median(St,'omitnan');
            D(i).S_t=interp1(datenum(tt),St,datenum(time),'linear',NaN);
        catch, end
    end

    % storm-active range
    va = isfinite(D(i).stress) | isfinite(D(i).nike_bb);
    if any(va), a=find(va,1,'first'); b=find(va,1,'last');
        keep=false(size(va)); keep(a:b)=true;
    else, keep=false(size(va)); end
    D(i).keep=keep;
    if any(keep), D(i).xmin=min(time(keep)); D(i).xmax=max(time(keep));
    else, D(i).xmin=min(time); D(i).xmax=max(time); end

    % [F7] 'phase', 'nike_ml' and 'nike_dp' deliberately EXCLUDED from the fill
    % list. Gap-filling is fine for continuous curves but must not feed a
    % statistic: interpolating them made the figure annotations disagree with
    % the post summaries (below-ML 52.1 vs 53.5 for Amphan, align 16.7 vs 16.9
    % for Fani). Linear interpolation of a phase ANGLE is also wrong in
    % principle, since angles wrap at +/-180 deg.
    fl={'stress','wp','ww_signed','ww_pos','nike_bb','nike_coh','dnike_dt', ...
        'MLD','ILD','BLT','cent','coh_ratio','nike_R','nike_L', ...
        'stress_R','stress_L','S_t','pen','peke','SST','OHC'};
    for j=1:numel(fl)
        if ~isfield(D(i),fl{j}), continue; end
        v=D(i).(fl{j});
        if isnumeric(v)&&numel(v)==numel(keep)&&any(keep)
            vf=v; vf(keep)=fillmissing(v(keep),'linear','EndValues','none');
            D(i).(fl{j})=vf;
        end
    end
end

%% DRIVER
if ismember(2,FIGURES_TO_MAKE), make_fig2(D,n_cyc,OUT_DIR,COL,beautify,savepro); end
if ismember(3,FIGURES_TO_MAKE), make_fig3(D,n_cyc,OUT_DIR,COL,LIT,beautify,savepro); end
if ismember(4,FIGURES_TO_MAKE), make_fig4(D,n_cyc,OUT_DIR,COL,beautify,savepro); end
if ismember(5,FIGURES_TO_MAKE), make_fig5(D,n_cyc,OUT_DIR,COL,LIT,beautify,savepro); end
if ismember(6,FIGURES_TO_MAKE), make_fig6(D,n_cyc,OUT_DIR,COL,LIT,beautify,savepro,Omega); end
if ismember(7,FIGURES_TO_MAKE), make_fig7(D,n_cyc,OUT_DIR,COL,LIT,beautify,savepro); end
% Figs 8 and 9 were moved to the Supporting Information as S10 and S11;
% they are still built here but written under the FigS10_/FigS11_ names.
if ismember(8,FIGURES_TO_MAKE), make_fig8(D,n_cyc,OUT_DIR,COL,LIT,beautify,savepro); end
if ismember(9,FIGURES_TO_MAKE), make_fig9(D,n_cyc,OUT_DIR,COL,beautify,savepro); end

fprintf('\n=== ALL FIGURES DONE. Output in %s ===\n', OUT_DIR);

%% SHARED HELPERS
function cb = shrink_cb(cb,lab)
    % [F14] Four-panel rows leave little room at the right edge, so the long
    % colorbar labels were being clipped. Shorten the label, put the units in
    % the tick area instead, and reduce the font so it fits.
    cb.Label.String = '';           % [F16] a rotated side label overlapped the ticks
    cb.Title.String  = lab;         % put it above the bar instead
    cb.Title.FontSize = 8;
    cb.Title.FontWeight = 'normal';
    cb.FontSize = 7;
end

function v = gof(M,f,n)
    if isfield(M,f)&&numel(M.(f))==n, v=M.(f)(:); else, v=nan(n,1); end
end
function A = foe(M,f)
    if isfield(M,f), A=M.(f); else, A=[]; end
end
function [zu, Mu] = unifdepth(zp, Mp, ZMAX, method)
    % [F2] imagesc spaces rows uniformly between y(1) and y(end), so a
    % non-uniform depth vector is drawn distorted. Resample onto a uniform
    % grid first. Mp is [ntime x ndepth]; returns Mu as [nzu x ntime],
    % which is the orientation imagesc expects.
    % zu starts at zp(1) rather than 0 so that a depth vector which does not
    % begin at the surface (e.g. VP.z) does not produce NaN rows at the top.
    if nargin < 4, method = 'linear'; end
    zp = zp(:);
    z0 = min(zp);
    if ~isfinite(z0), z0 = 0; end
    zu = linspace(z0, ZMAX, 200).';
    Mu = interp1(zp, Mp.', zu, method);
end
function v = lastfin(x)
    % [F1] Final finite value of a cumulative series. Cumulative wind work is
    % read at the end of the record, not at its maximum: a curve that ends
    % below its own peak (net negative transfer) is otherwise reported at the
    % peak, which returned 0.00 for Fani instead of its true net value.
    i = find(isfinite(x),1,'last');
    if isempty(i), v = NaN; else, v = x(i); end
end
function set_xax_date(ax)
    try, datetick(ax,'x','dd','keeplimits'); catch, end
end
function Mf = fill_nan_cols(M)
    % Interpolate across all-NaN timestep rows (dropped/coastal) so imagesc
    % shows no dark vertical stripes. M is [nt x nz]; rows = timesteps.
    % NOTE: this is cosmetic gap-filling and should be stated in the caption.
    Mf = M;
    if isempty(M), return; end
    bad = all(~isfinite(M),2);              % all-NaN timesteps
    good = ~bad;
    if sum(good) < 2 || ~any(bad), return; end
    idx = (1:size(M,1))';
    for c = 1:size(M,2)
        col = M(:,c);
        g = isfinite(col);
        if sum(g) >= 2
            Mf(:,c) = interp1(idx(g), col(g), idx, 'linear', NaN);
        end
    end
end
function make_fig2(D, n_cyc, OUT_DIR, C, beautify, savepro)
    fprintf('\n--- Fig 2 : Storm forcing ---\n');
    figure('Name','Fig 2 - Forcing','Units','normalized','Position',[0.04 0.05 0.92 0.88]);
    for i = 1:n_cyc
        d = D(i); k = d.keep;
        % (a) wind stress
        ax = subplot(3,n_cyc,i);
        if any(k)
            t=d.time(k);
            fill([t;flipud(t)],[d.stress(k);zeros(numel(t),1)],C.stress, ...
                'FaceAlpha',0.25,'EdgeColor','none'); hold on;
            plot(t,d.stress(k),'Color',C.stress,'LineWidth',1.8);
            xlim([d.xmin d.xmax]);
            ym=max(d.stress(k),[],'omitnan')*1.15; if ~isfinite(ym)||ym<=0,ym=0.2;end
            ylim([0 ym]);
            text(0.97,0.97, ...
                sprintf('mean=%.3f\npeak=%.3f',mean(d.stress,'omitnan'),max(d.stress,[],'omitnan')), ...
                'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
                'FontSize',8,'BackgroundColor','w','Margin',2,'EdgeColor',[0.7 0.7 0.7]);
        end
        title(sprintf('%s (%s)',d.name,d.basin),'FontSize',11);
        if i==1, ylabel({'\bf(a)';'\tau (N m^{-2})'},'FontSize',10); else, ylabel('\tau (N m^{-2})'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
        % (b) wind power tau.u_I (signed, pos green / neg red)
        ax = subplot(3,n_cyc,n_cyc+i);
        if any(k)
            t=d.time(k); wpk=d.wp(k);
            pos=wpk; pos(pos<0)=0; neg=wpk; neg(neg>0)=0;
            area(t,pos,'FaceColor',C.wpos,'FaceAlpha',0.5,'EdgeColor','none'); hold on;
            area(t,neg,'FaceColor',C.wneg,'FaceAlpha',0.5,'EdgeColor','none');
            plot(t,wpk,'Color',[0.2 0.2 0.2],'LineWidth',1.2);
            yline(0,'k-','LineWidth',0.6);
            xlim([d.xmin d.xmax]);
            % [F3] '\cdot u' with a space; MATLAB TeX cannot parse '\cdotu'.
            if i==1, legend({'\tau\cdot u_I > 0','\tau\cdot u_I < 0'},'Box','off','Location','best','FontSize',8); end
        end
        if i==1, ylabel({'\bf(b)';'\tau\cdot u_I (W m^{-2})'},'FontSize',10); else, ylabel('\tau\cdot u_I (W m^{-2})'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
        % (c) cumulative wind work: signed and positive
        ax = subplot(3,n_cyc,2*n_cyc+i);
        if any(k)
            t=d.time(k);
            plot(t,d.ww_pos(k),'Color',C.work,'LineWidth',2.0); hold on;
            plot(t,d.ww_signed(k),'--','Color',C.work2,'LineWidth',1.8);
            yline(0,'k:','LineWidth',0.6,'HandleVisibility','off');
            xlim([d.xmin d.xmax]);
            % annotation reads the END of the cumulative curve, not its peak
            text(0.97,0.03,sprintf('pos=%.2f\nsigned=%+.2f',lastfin(d.ww_pos),lastfin(d.ww_signed)), ...
                'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
                'FontSize',8,'BackgroundColor','w','Margin',2,'EdgeColor',[0.7 0.7 0.7]);
            if i==1, legend({'positive','signed'},'Box','off','Location','best','FontSize',8); end
        end
        if i==1, ylabel({'\bf(c)';'Cum WW (kJ m^{-2})'},'FontSize',10); else, ylabel('Cum WW (kJ m^{-2})'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
    end
    savepro('Fig02_Storm_Forcing.png'); fprintf(' Saved Fig 2\n');
    fprintf('\n[FIG 2 - FORCING]\n%-8s | %-9s | %-9s | %-9s | %-9s\n', ...
        'Storm','MeanTau','WWpos','WWsign','Retain%');
    for i=1:n_cyc, d=D(i);
        wp_=lastfin(d.ww_pos); ws_=lastfin(d.ww_signed);
        fprintf('%-8s | %-9.4f | %-+9.2f | %-+9.2f | %-9.0f\n',d.name, ...
            mean(d.stress,'omitnan'),wp_,ws_,100*ws_/max(abs(wp_),eps));
    end
end
function make_fig3(D, n_cyc, OUT_DIR, C, LIT, beautify, savepro)
    fprintf('\n--- Fig 3 : NIKE / S / phase ---\n');
    figure('Name','Fig 3','Units','normalized','Position',[0.04 0.05 0.92 0.88]);
    for i = 1:n_cyc
        d = D(i); k = d.keep;
        % (a) NIKE coherent + broadband
        ax = subplot(3,n_cyc,i);
        if any(k)
            t=d.time(k);
            fill([t;flipud(t)],[d.nike_bb(k);zeros(numel(t),1)],C.bb,'FaceAlpha',0.25,'EdgeColor','none'); hold on;
            fill([t;flipud(t)],[d.nike_coh(k);zeros(numel(t),1)],C.coh,'FaceAlpha',0.30,'EdgeColor','none');
            plot(t,d.nike_bb(k),'Color',C.bb,'LineWidth',1.8);
            plot(t,d.nike_coh(k),'Color',C.coh,'LineWidth',1.8);
            yl=max([d.nike_bb(k);d.nike_coh(k)],[],'omitnan')*1.15; if ~isfinite(yl)||yl<=0,yl=1;end
            xlim([d.xmin d.xmax]); ylim([0 yl]);
            % [F11] normalized units, top-RIGHT, so the box cannot collide with the
            % top-left legend or ride on the data (R2 figure comments).
            text(0.97,0.97, ...
                sprintf('Coh=%.1f\nWin=%.1f',max(d.nike_coh,[],'omitnan'),max(d.nike_bb,[],'omitnan')), ...
                'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
                'FontSize',8,'BackgroundColor','w','Margin',2,'EdgeColor',[0.7 0.7 0.7]);
            % [F10] 'Win' = windowed kinetic energy (Eq. 6): the detrended windowed
            % variance, NOT a bandpassed quantity. Label must match the manuscript.
            if i==1, legend({'Windowed KE','Coherent NIKE'},'Box','off','Location','northwest','FontSize',8); end
        end
        title(sprintf('%s (%s)',d.name,d.basin),'FontSize',11);
        if i==1, ylabel({'\bf(a)';'NIKE (kJ m^{-2})'}); else, ylabel('NIKE (kJ m^{-2})'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
        % (b) S parameter
        ax = subplot(3,n_cyc,n_cyc+i);
        if any(k)&&any(isfinite(d.S_t))
            Sk=d.S_t(k); Sk(Sk<=0|~isfinite(Sk))=NaN;
            xb=[d.xmin d.xmax d.xmax d.xmin];
            patch(xb,[0.1 0.1 LIT.S_lo LIT.S_lo],C.zoneSub,'EdgeColor','none','FaceAlpha',0.7,'HandleVisibility','off'); hold on;
            patch(xb,[LIT.S_lo LIT.S_lo LIT.S_hi LIT.S_hi],C.zoneRes,'EdgeColor','none','FaceAlpha',0.7,'HandleVisibility','off');
            patch(xb,[LIT.S_hi LIT.S_hi 100 100],C.zoneSup,'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');
            plot(d.time(k),Sk,'Color',C.S,'LineWidth',2.0);
            xlim([d.xmin d.xmax]); set(gca,'YScale','log');
            Sv=Sk(isfinite(Sk)&Sk>0);
            if isempty(Sv), ylim([0.2 10]); else, ylim([max(0.15,min(Sv)*0.7) max([max(Sv)*1.3 LIT.S_hi*1.1])]); end
            if isfinite(d.S_median_trk)
                % [F11] was pinned to y=0.25 on a log axis whose limits vary per storm,
                % so it fell outside the axes in three of four panels.
                text(0.03,0.06,sprintf('Med S=%.2f',d.S_median_trk),'Units','normalized', ...
                    'HorizontalAlignment','left','VerticalAlignment','bottom', ...
                    'FontSize',9,'FontWeight','bold','BackgroundColor','w','Margin',2,'EdgeColor',[0.5 0.5 0.5]);
            end
        end
        if i==1, ylabel({'\bf(b)';'S = U_h/(|f|R_{max})'}); else, ylabel('S'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
        % (c) tau-u_I phase
        ax = subplot(3,n_cyc,2*n_cyc+i);
        if any(k)&&any(isfinite(d.phase))
            t=d.time(k); ph=d.phase(k);
            xb=[d.xmin d.xmax d.xmax d.xmin];
            patch(xb,[-45 -45 45 45],C.zoneRes,'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off'); hold on;
            yline(0,'k-','LineWidth',0.8,'HandleVisibility','off');
            yline(45,'k--','LineWidth',0.7,'HandleVisibility','off');
            yline(-45,'k--','LineWidth',0.7,'HandleVisibility','off');
            plot(t,ph,'Color',C.S,'LineWidth',1.6);
            xlim([d.xmin d.xmax]); ylim([-180 180]); set(gca,'YTick',-180:90:180);
            mp=median(d.phase,'omitnan'); fr=100*sum(abs(d.phase)<45,'omitnan')/max(sum(isfinite(d.phase)),1);
            text(0.97,0.97,sprintf('Med=%.0f\\circ\nAlign=%.0f%%',mp,fr),'Units','normalized', ...
                'HorizontalAlignment','right','VerticalAlignment','top', ...
                'FontSize',9,'FontWeight','bold','BackgroundColor','w','Margin',2,'EdgeColor',[0.5 0.5 0.5]);
            if i==1
                hb=patch(NaN,NaN,C.zoneRes,'EdgeColor','none','FaceAlpha',0.5);
                legend(hb,{'Resonant |\Delta\phi|<45\circ'},'Box','off','Location','south','FontSize',8);
            end
        end
        if i==1, ylabel({'\bf(c)';'\Delta\phi(\tau,u_I) (deg)'}); else, ylabel('\Delta\phi (deg)'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
    end
    savepro('Fig03_NIKE_S_Phase.png'); fprintf(' Saved Fig 3\n');
    fprintf('\n[FIG 3 - NIKE/S/PHASE]\n%-8s | %-8s | %-8s | %-8s | %-8s | %-7s | %-6s\n', ...
        'Storm','PkWin','PkCoh','MedS','MedPh','Align%','Nvalid');
    for i=1:n_cyc, d=D(i);
        nv=sum(isfinite(d.phase));
        fprintf('%-8s | %-8.2f | %-8.2f | %-8.2f | %-8.0f | %-7.1f | %-6d\n',d.name, ...
            max(d.nike_bb,[],'omitnan'),max(d.nike_coh,[],'omitnan'),d.S_median_trk, ...
            median(d.phase,'omitnan'),100*sum(abs(d.phase)<45,'omitnan')/max(nv,1),nv);
    end
end
function make_fig4(D, n_cyc, OUT_DIR, C, beautify, savepro)
    fprintf('\n--- Fig 4 : Vertical NIKE ---\n');
    figure('Name','Fig 4','Units','normalized','Position',[0.04 0.05 0.92 0.85]);
    ZMAX=400; use_cmo=exist('cmocean.m','file');
    for i=1:n_cyc
        d=D(i); k=d.keep;
        if isempty(d.NIKE_coh_z), continue; end
        iz=d.z<=ZMAX; zp=d.z(iz);
        Mc=fill_nan_cols(d.NIKE_coh_z(:,iz)); Mb=fill_nan_cols(d.NIKE_brd_z(:,iz)); Mr=fill_nan_cols(d.Coh_z(:,iz));
        % (a) coherent
        ax=subplot(3,n_cyc,i);
        if any(k)
            t=d.time(k); Mp=Mc(k,:); cl=prctile(Mp(:),98); if ~isfinite(cl)||cl<=0,cl=0.3;end
            [zu,Mu]=unifdepth(zp,Mp,ZMAX);            % [F2] uniform depth axis
            % [F19] log colour scale. On a linear scale the sub-200 m signal is
            % invisible, which is what prompted R1's objection. Two and a half decades.
            Mu(Mu<=0)=NaN;
            imagesc(datenum(t),zu,log10(Mu)); set(gca,'YDir','reverse');
            caxis([log10(cl)-2.5 log10(cl)]);
            if use_cmo, colormap(gca,cmocean('thermal')); else, colormap(gca,hot); end
            hold on; plot(datenum(t),d.MLD(k),'w-','LineWidth',1.5);
            plot(datenum(t),d.cent(k),'k--','LineWidth',1.5);
            xlim(datenum([d.xmin d.xmax])); ylim([0 ZMAX]); set_xax_date(ax);
        end
        title(sprintf('%s (%s)',d.name,d.basin),'FontSize',11);
        if i==1, ylabel({'\bf(a)';'Depth (m)'}); else, ylabel('Depth (m)'); end
                xlabel('Date (UTC)'); cb=colorbar; cb=shrink_cb(cb,'log_{10} KE_{win} (kJ m^{-3})'); beautify();
        % (b) broadband
        ax=subplot(3,n_cyc,n_cyc+i);
        if any(k)
            t=d.time(k); Mp=Mb(k,:); cl=prctile(Mp(:),98); if ~isfinite(cl)||cl<=0,cl=0.3;end
            [zu,Mu]=unifdepth(zp,Mp,ZMAX);
            % [F19] log colour scale. On a linear scale the sub-200 m signal is
            % invisible, which is what prompted R1's objection. Two and a half decades.
            Mu(Mu<=0)=NaN;
            imagesc(datenum(t),zu,log10(Mu)); set(gca,'YDir','reverse');
            caxis([log10(cl)-2.5 log10(cl)]);
            if use_cmo, colormap(gca,cmocean('thermal')); else, colormap(gca,hot); end
            hold on; plot(datenum(t),d.MLD(k),'w-','LineWidth',1.5);
            xlim(datenum([d.xmin d.xmax])); ylim([0 ZMAX]); set_xax_date(ax);
        end
        if i==1, ylabel({'\bf(b)';'Depth (m)'}); else, ylabel('Depth (m)'); end
                xlabel('Date (UTC)'); cb=colorbar; cb=shrink_cb(cb,'log_{10} KE_{win} (kJ m^{-3})'); beautify();
        % (c) coherence ratio (uncapped)
        ax=subplot(3,n_cyc,2*n_cyc+i);
        if any(k)&&~isempty(Mr)
            t=d.time(k); Mp=Mr(k,:);
            [zu,Mu]=unifdepth(zp,Mp,ZMAX);
            imagesc(datenum(t),zu,Mu); set(gca,'YDir','reverse'); caxis([0 1]); colormap(gca,parula);
            hold on; plot(datenum(t),d.MLD(k),'w-','LineWidth',1.5);
            xlim(datenum([d.xmin d.xmax])); ylim([0 ZMAX]); set_xax_date(ax);
        end
        if i==1, ylabel({'\bf(c)';'Depth (m)'}); else, ylabel('Depth (m)'); end
        xlabel('Date (UTC)'); cb=colorbar; cb=shrink_cb(cb,'R_C'); beautify();
    end
    savepro('Fig04_Vertical_NIKE.png'); fprintf(' Saved Fig 4\n');
    fprintf('\n[FIG 4 - VERTICAL NIKE]\n%-8s | %-7s | %-9s | %-9s | %-7s\n', ...
        'Storm','MLDmn','PkCohSurf','PkCohDeep','Pen(m)');
    for i=1:n_cyc, d=D(i);
        if isempty(d.NIKE_coh_z), continue; end
        ML=mean(d.MLD,'omitnan'); is=d.z<=ML; id=d.z>ML&d.z<=400;
        ps=max(d.NIKE_coh_z(:,is),[],'all','omitnan'); pd=max(d.NIKE_coh_z(:,id),[],'all','omitnan');
        cm=mean(d.NIKE_coh_z,1,'omitnan'); thr=0.10*cm(1);
        ip=find(cm<thr&d.z(:)'<=400,1,'first'); if isempty(ip),pen=400;else,pen=d.z(ip);end
        fprintf('%-8s | %-7.1f | %-9.3f | %-9.3f | %-7.1f\n',d.name,ML,ps,pd,pen);
    end
end
%  FIG 5 - VERTICAL PROPAGATION (replaces WKB Cgz; R1#4)
%  Row (a) NIKE_coh depth-time to 400 m with 200 m line (penetration)
%  Row (b) phase-tilt panel from *_vertprop.mat (mechanism), fallback text
%  Row (c) single: per-storm below-200m penetration fraction bar chart
function make_fig5(D, n_cyc, OUT_DIR, C, LIT, beautify, savepro)
    fprintf('\n--- Fig 5 : Vertical propagation ---\n');
    figure('Name','Fig 5','Units','normalized','Position',[0.04 0.04 0.92 0.92]);
    use_cmo=exist('cmocean.m','file'); ZMAX=400;
    frac200=nan(1,n_cyc); lamz=nan(1,n_cyc);
    for i=1:n_cyc
        d=D(i); k=d.keep;
        % (a) NIKE penetration depth-time
        ax=subplot(3,n_cyc,i);
        if any(k)&&~isempty(d.NIKE_coh_z)
            iz=d.z<=ZMAX; zp=d.z(iz); t=d.time(k); Mp=fill_nan_cols(d.NIKE_coh_z(k,iz));
            cl=prctile(Mp(:),98); if ~isfinite(cl)||cl<=0,cl=0.3;end
            [zu,Mu]=unifdepth(zp,Mp,ZMAX);            % [F2] uniform depth axis
            % [F19] log colour scale. On a linear scale the sub-200 m signal is
            % invisible, which is what prompted R1's objection. Two and a half decades.
            Mu(Mu<=0)=NaN;
            imagesc(datenum(t),zu,log10(Mu)); set(gca,'YDir','reverse');
            caxis([log10(cl)-2.5 log10(cl)]);
            if use_cmo, colormap(gca,cmocean('thermal')); else, colormap(gca,hot); end
            hold on; plot(datenum(t),d.MLD(k),'w-','LineWidth',1.5);
            yline(200,'c--','LineWidth',1.3,'HandleVisibility','off');
            xlim(datenum([d.xmin d.xmax])); ylim([0 ZMAX]); set_xax_date(ax);
            % fraction of coherent energy below 200 m
            % NOTE: this is the ratio of TIME-MEAN profiles (trapz of the mean
            % profile), which is a different convention from the Fig 6 below-ML
            % percentage (time-mean of the instantaneous ratio). State both.
            cmn=mean(d.NIKE_coh_z,1,'omitnan'); z=d.z(:);
            below=z>200&z<=700; tot=safe_trapz2(z,cmn); bel=safe_trapz2(z(below),cmn(below));
            frac200(i)=bel/(tot+eps);
        end
        title(sprintf('%s (%s)',d.name,d.basin),'FontSize',11);
        if i==1, ylabel({'\bf(a)';'Depth (m)'}); else, ylabel('Depth (m)'); end
                xlabel('Date (UTC)'); cb=colorbar; cb=shrink_cb(cb,'log_{10} NIKE_{coh} (kJ m^{-3})'); beautify();
        % (b) phase tilt from vertprop
        ax=subplot(3,n_cyc,n_cyc+i);
        if d.has_vp && isfield(d.vp,'phase')
            VP=d.vp; zc=VP.z(:); izc=zc<=200; ph=VP.phase(:,izc);
            % 'nearest' so the +/-pi wrap is not interpolated into false colours
            [zub,Pu]=unifdepth(zc(izc),ph,200,'nearest');
            imagesc(1:size(ph,1),zub,Pu); set(gca,'YDir','reverse');
            colormap(gca,hsv); caxis([-pi pi]);
            xlabel('time index'); ylim([0 200]);
            % [F13] A rejected fit leaves an unannotated phase field, which reads
            % as a broken panel rather than as a reported non-result. Label it.
            mm=median(VP.mslopes,'omitnan');
            % [F8] Explicit acceptance threshold. The phase slope is fitted over
            % the upper 200 m, so a returned wavelength far exceeding that
            % domain is not resolved by the fit. LAMZ_MAX = 500 m is 2.5x the
            % fit domain: it accepts Kyarr (321), Amphan (401) and Fani (329)
            % and rejects Tauktae (1117), which previously printed as if it were
            % a measurement while Table 2 reported it as not resolved.
            LAMZ_MAX = 500;
            if isfinite(mm) && abs(mm) > 2*pi/LAMZ_MAX, lamz(i)=2*pi/abs(mm); end
            if isfinite(lamz(i))
                text(0.03,0.9,sprintf('\\lambda_z=%.0f m',lamz(i)),'Units','normalized', ...
                    'FontSize',9,'FontWeight','bold','BackgroundColor','w','Margin',2,'EdgeColor',[0.5 0.5 0.5]);
            else
                text(0.5,0.5,{'\lambda_z not resolved';'(fit fails criterion)'}, ...
                    'Units','normalized','HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontSize',10,'FontWeight','bold','Color',[0.15 0.15 0.15], ...
                    'BackgroundColor','w','Margin',4,'EdgeColor',[0.35 0.35 0.35]);
            end
            cb=colorbar; cb.Label.String='phase (rad)';
        else
            text(0.5,0.5,'vertprop n/a','Units','normalized','HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
            axis off;
        end
        if i==1, ylabel({'\bf(b)';'Depth (m)'}); end
        beautify();
    end
    % (c) below-200m penetration fraction bar chart
    ax=subplot(3,1,3);
    bar(1:n_cyc,frac200*100,0.55,'FaceColor',[0.30 0.55 0.75],'EdgeColor','k'); hold on;
    for i=1:n_cyc
        if isfinite(frac200(i))
            text(i,frac200(i)*100+0.5,sprintf('%.1f%%',frac200(i)*100), ...
                'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
        end
        % [F17] the lambda_z annotations sat on top of the storm tick labels;
        % fold them into a two-line tick label instead.
    end
    % [F18] Built by concatenation, not sprintf. In a sprintf format string
    % MATLAB reads '\newline' as the escape '\n' followed by "ewline" and then
    % warns on the invalid '\l', which mangled these labels. A plain string
    % literal passes the backslashes through to the TeX interpreter untouched.
    xl = cell(1,n_cyc);
    for i = 1:n_cyc
        if isfinite(lamz(i))
            xl{i} = [D(i).name '\newline\lambda_z=' num2str(round(lamz(i))) ' m'];
        else
            xl{i} = [D(i).name '\newline\lambda_z n.r.'];
        end
    end
    set(ax,'XTick',1:n_cyc,'XTickLabel',xl,'TickLabelInterpreter','tex');
    ylabel({'\bf(c)';'NIKE below 200 m (%)'}); xlabel('Storm');
        title('Deep penetration fraction (coherent NIKE below 200 m) with vertical wavelength', ...
        'FontSize',11,'Units','normalized','Position',[0.5 1.06 0]);
    ym=max(frac200*100,[],'omitnan')*1.3; if ~isfinite(ym)||ym<=0,ym=10;end
    ylim([0 ym]); beautify();
    savepro('Fig05_Vertical_Propagation.png'); fprintf(' Saved Fig 5\n');
    fprintf('\n[FIG 5 - VERTICAL PROPAGATION]\n%-8s | %-11s | %-10s\n','Storm','Below200m%','LambdaZ(m)');
    for i=1:n_cyc
        fprintf('%-8s | %-11.1f | %-10.0f\n',D(i).name,frac200(i)*100,lamz(i));
    end
end
function I=safe_trapz2(z,p)
    z=z(:); p=p(:); v=isfinite(p)&isfinite(z);
    if sum(v)>=2, I=trapz(z(v),p(v)); else, I=NaN; end
end
%  FIG 6 - ENERGY PARTITION + ROTARY SPECTRUM (PE/KE dropped)
%  Row (a) NIKE_ML vs NIKE_deep stacked
%  Row (b) rotary spectrum CW vs CCW at surface
function make_fig6(D, n_cyc, OUT_DIR, C, LIT, beautify, savepro, Omega)
    fprintf('\n--- Fig 6 : Partition + rotary ---\n');
    c_ml=[0.20 0.55 0.80]; c_deep=[0.85 0.45 0.10]; c_line=[0.05 0.20 0.40];
    figure('Name','Fig 6','Units','normalized','Position',[0.04 0.06 0.92 0.82]);
    for i=1:n_cyc
        d=D(i); k=d.keep;
        % (a) partition
        ax=subplot(2,n_cyc,i);
        if any(k)
            t=d.time(k); ml=d.nike_ml(k); dp=d.nike_dp(k);
            ml(~isfinite(ml))=0; dp(~isfinite(dp))=0;
            ha=area(t,[dp(:) ml(:)]);
            ha(1).FaceColor=c_deep; ha(1).FaceAlpha=0.55; ha(1).EdgeColor='none';
            ha(2).FaceColor=c_ml;   ha(2).FaceAlpha=0.55; ha(2).EdgeColor='none';
            hold on; plot(t,d.nike_coh(k),'Color',c_line,'LineWidth',2.0);
            xlim([d.xmin d.xmax]);
            yl=max(d.nike_coh(k),[],'omitnan')*1.15; if ~isfinite(yl)||yl<=0,yl=max(ml+dp)*1.1;end
            if ~isfinite(yl)||yl<=0,yl=1; end
            ylim([0 yl]);
            % time-mean of the instantaneous ratio (NOT the ratio of time-means;
            % Fig 5 frac200 uses the other convention -- state both in captions)
            % [F9] Defer to the post statistic so the panel annotation matches
            % Table 2. Post applies its own validity guards, so recomputing the
            % ratio here over every finite pair gives a slightly different
            % answer (52.1 vs 53.5 for Amphan). Post is authoritative.
            if d.has_post && isfield(d.stat,'below_ML_frac')
                pct = d.stat.below_ML_frac;
            else
                pct = 100*mean(d.nike_dp./(d.nike_ml+d.nike_dp+eps),'omitnan');
            end
            % [F15] moved to the bottom right; the legend occupies the top of
            % this panel and the two boxes were overlapping in the Kyarr panel.
            text(0.97,0.04,sprintf('Below-ML=%.1f%%',pct),'Units','normalized', ...
                'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',8, ...
                'BackgroundColor','w','Margin',2,'EdgeColor',[0.7 0.7 0.7]);
            if i==1, legend({'deep','ML','coh total'},'Box','off','Location','northwest','FontSize',7); end
        end
        title(sprintf('%s (%s)',d.name,d.basin),'FontSize',11);
        if i==1, ylabel({'\bf(a)';'NIKE (kJ m^{-2})'}); else, ylabel('NIKE (kJ m^{-2})'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
        % (b) rotary spectrum
        ax=subplot(2,n_cyc,n_cyc+i);
        f_in=abs(2*Omega*sind(d.lat_mean))*86400/(2*pi);
        if ~isempty(d.u_surf)&&~isempty(d.v_surf)
            u=fillmissing(d.u_surf,'linear'); v=fillmissing(d.v_surf,'linear');
            ok=isfinite(u)&isfinite(v); u=u(ok); v=v(ok); N=numel(u);
            if N>=12
                Fs=24/3; fn=Fs/2;
                flo=0.5*f_in/fn; fhi=min(1.5*f_in/fn,0.95);
                if fhi>flo&&flo>0
                    [b,a]=butter(2,[flo fhi],'bandpass'); u=filtfilt(b,a,u); v=filtfilt(b,a,v);
                else
                    u=u-movmean(u,8,'omitnan'); v=v-movmean(v,8,'omitnan');
                end
                Nw=numel(u); w=hann(Nw)'; W=(u(:)+1i*v(:)).*w(:);
                Y=fftshift(fft(W)); df=Fs/Nw; fv=(-Nw/2:Nw/2-1)*df;
                P=abs(Y).^2/(sum(w.^2)*df); P(P<=0)=NaN;
                fcw=abs(fv(fv<0)); Pcw=P(fv<0); fccw=fv(fv>0); Pccw=P(fv>0);
                [fcw,ic]=sort(fcw); Pcw=Pcw(ic); [fccw,icc]=sort(fccw); Pccw=Pccw(icc);
                semilogy(fcw,Pcw,'-','Color',C.cw,'LineWidth',2.2); hold on;
                semilogy(fccw,Pccw,'--','Color',C.ccw,'LineWidth',1.5);
                xline(f_in,'k-','LineWidth',1.2,'HandleVisibility','off');
                % [F5] M2 shown as a reference frequency only. HYCOM GLBy0.08
                % is non-tidal, so no tidal energy is present at this line;
                % it marks that the retained band does not reach M2.
                xline(1.9323,'--','Color',[0.5 0.28 0.68],'LineWidth',0.9,'HandleVisibility','off');
                text(f_in*1.03,max(Pcw)*0.7,'f','FontSize',9);
                text(1.9323*1.02,max(Pcw)*0.15,'M_2','FontSize',8,'Color',[0.5 0.28 0.68]);
                xlim([0.05 4]);
                if any(isfinite(Pcw)&Pcw>0)
                    yl=[max(min(Pcw(Pcw>0))/2,1e-6) max(Pcw)*5];
                    if all(isfinite(yl))&&yl(2)>yl(1), ylim(yl); end
                end
                if i==1, legend({'CW (NIW)','CCW'},'Box','off','Location','best','FontSize',8); end
            end
        end
        xlabel('Frequency (cpd)');
        if i==1, ylabel({'\bf(b)';'Rotary power'}); else, ylabel('Rotary power'); end
        beautify();
    end
    savepro('Fig06_Partition_Rotary.png'); fprintf(' Saved Fig 6\n');
    fprintf('\n[FIG 6 - PARTITION+ROTARY]\n%-8s | %-8s | %-8s | %-9s\n','Storm','MeanML','MeanDp','BelowML%');
    for i=1:n_cyc, d=D(i);
        if d.has_post && isfield(d.stat,'below_ML_frac')
            pct = d.stat.below_ML_frac;
        else
            pct = 100*mean(d.nike_dp./(d.nike_ml+d.nike_dp+eps),'omitnan');
        end
        fprintf('%-8s | %-8.2f | %-8.2f | %-9.1f\n',d.name,mean(d.nike_ml,'omitnan'),mean(d.nike_dp,'omitnan'),pct);
    end
end
function make_fig7(D, n_cyc, OUT_DIR, C, LIT, beautify, savepro)
    fprintf('\n--- Fig 7 : R/L asymmetry ---\n');
    figure('Name','Fig 7','Units','normalized','Position',[0.04 0.02 0.92 0.95]);
    for i=1:n_cyc
        d=D(i); k=d.keep;
        % (a) NIKE R vs L
        ax=subplot(4,n_cyc,i);
        if any(k)
            t=d.time(k);
            fill([t;flipud(t)],[d.nike_R(k);zeros(numel(t),1)],C.right,'FaceAlpha',0.30,'EdgeColor','none'); hold on;
            fill([t;flipud(t)],[d.nike_L(k);zeros(numel(t),1)],C.left,'FaceAlpha',0.30,'EdgeColor','none');
            plot(t,d.nike_R(k),'Color',C.right,'LineWidth',1.8);
            plot(t,d.nike_L(k),'Color',C.left,'LineWidth',1.8);
            ym=max([d.nike_R(k);d.nike_L(k)],[],'omitnan')*1.15; if ~isfinite(ym)||ym<=0,ym=1;end
            xlim([d.xmin d.xmax]); ylim([0 ym]);
            if i==1, legend({'R','L'},'Box','off','Location','best','FontSize',8); end
        end
        title(sprintf('%s (%s)',d.name,d.basin),'FontSize',11);
        if i==1, ylabel({'\bf(a)';'NIKE (kJ m^{-2})'}); else, ylabel('NIKE (kJ m^{-2})'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
        % (b) R/L ratio
        ax=subplot(4,n_cyc,n_cyc+i);
        if any(k)
            t=d.time(k); ar=d.AR(k); ar(ar<=0|~isfinite(ar))=NaN;
            xb=[d.xmin d.xmax d.xmax d.xmin];
            patch(xb,[0.1 0.1 0.67 0.67],[0.85 0.92 1.00],'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off'); hold on;
            patch(xb,[0.85 0.85 1.18 1.18],[0.95 0.95 0.95],'EdgeColor','none','FaceAlpha',0.6,'HandleVisibility','off');
            patch(xb,[1.5 1.5 10 10],[1.00 0.88 0.85],'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
            yline(1,'k-','LineWidth',1.0,'HandleVisibility','off');
            plot(t,ar,'Color',[0.30 0.10 0.50],'LineWidth',2.0);
            xlim([d.xmin d.xmax]); ylim([0.2 10]); set(gca,'YScale','log');
            text(0.97,0.97,sprintf('Med R/L=%.2f',median(d.AR,'omitnan')),'Units','normalized', ...
                'HorizontalAlignment','right','VerticalAlignment','top', ...
                'FontSize',9,'FontWeight','bold','BackgroundColor','w','Margin',2,'EdgeColor',[0.5 0.5 0.5]);
        end
        if i==1, ylabel({'\bf(b)';'R/L (log)'}); else, ylabel('R/L'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
        % (c) penetration per side
        ax=subplot(4,n_cyc,2*n_cyc+i);
        pR=mean(d.pen_R,'omitnan'); pL=mean(d.pen_L,'omitnan');
        sR=std(d.pen_R,'omitnan'); sL=std(d.pen_L,'omitnan');
        if isfinite(pR)||isfinite(pL)
            b=bar([1 2],[pR pL],0.5,'FaceColor','flat','EdgeColor','k'); b.CData=[C.right;C.left];
            hold on; errorbar([1 2],[pR pL],[sR sL],'k','LineStyle','none','CapSize',6);
            set(gca,'XTick',[1 2],'XTickLabel',{'R','L'}); xlim([0.5 2.5]);
            ylim([0 max([pR+sR pL+sL 0.45])*1.2]);
        end
        if i==1, ylabel({'\bf(c)';'Penetration'}); else, ylabel('Penetration'); end
        xlabel('Side'); beautify();
        % (d) stress R vs L
        ax=subplot(4,n_cyc,3*n_cyc+i);
        if any(k)
            t=d.time(k); sR_=d.stress_R(k); sL_=d.stress_L(k);
            ym=max([sR_;sL_],[],'omitnan')*1.2; if ~isfinite(ym)||ym<=0,ym=0.2;end
            fill([t;flipud(t)],[sR_;zeros(numel(t),1)],C.right,'FaceAlpha',0.28,'EdgeColor','none'); hold on;
            fill([t;flipud(t)],[sL_;zeros(numel(t),1)],C.left,'FaceAlpha',0.28,'EdgeColor','none');
            plot(t,sR_,'Color',C.right,'LineWidth',2.0); plot(t,sL_,'Color',C.left,'LineWidth',2.0);
            xlim([d.xmin d.xmax]); ylim([0 ym]);
            tR=mean(d.stress_R,'omitnan'); tL=mean(d.stress_L,'omitnan');
            text(0.97,0.97,sprintf('\\tau_R/\\tau_L=%.2f',tR/(tL+eps)),'Units','normalized', ...
                'HorizontalAlignment','right','VerticalAlignment','top', ...
                'FontSize',8,'FontWeight','bold','BackgroundColor','w','Margin',2,'EdgeColor',[0.5 0.5 0.5]);
            if i==1, legend({'\tau_R','\tau_L'},'Box','off','Location','best','FontSize',8); end
        end
        if i==1, ylabel({'\bf(d)';'\tau (N m^{-2})'}); else, ylabel('\tau (N m^{-2})'); end
        xlabel('Date (UTC)'); set_xax_date(ax); beautify();
    end
    savepro('Fig07_RightLeft_Asymmetry.png'); fprintf(' Saved Fig 7\n');
    fprintf('\n[FIG 7 - R/L ASYMMETRY]\n%-8s | %-8s | %-8s | %-8s | %-8s\n','Storm','MeanR','MeanL','MedR/L','TauR/L');
    for i=1:n_cyc, d=D(i);
        fprintf('%-8s | %-8.2f | %-8.2f | %-8.2f | %-8.2f\n',d.name, ...
            mean(d.nike_R,'omitnan'),mean(d.nike_L,'omitnan'),median(d.AR,'omitnan'), ...
            mean(d.stress_R,'omitnan')/(mean(d.stress_L,'omitnan')+eps));
    end
end
%  FIG 8 - CROSS-STORM SCORECARD (conversion efficiency removed; R2#6)
%  Metrics: R/L, |RC|@f, S, penetration, mean phase-alignment
function make_fig8(D, n_cyc, OUT_DIR, C, LIT, beautify, savepro)
    fprintf('\n--- Fig 8 : Cross-storm scorecard ---\n');
    figure('Name','Fig 8','Units','normalized','Position',[0.04 0.10 0.92 0.70]);
    M.RL=nan(1,n_cyc); M.RC=nan(1,n_cyc); M.S=nan(1,n_cyc);
    M.pen=nan(1,n_cyc); M.align=nan(1,n_cyc);
    for i=1:n_cyc
        d=D(i);
        M.RL(i)=median(d.AR,'omitnan');
        M.RC(i)=compute_RC_at_f(d);
        M.S(i)=d.S_median_trk;
        M.pen(i)=mean(d.pen,'omitnan');
        M.align(i)=sum(abs(d.phase)<45,'omitnan')/max(sum(isfinite(d.phase)),1);
    end
    metrics={ 'R/L',M.RL,[LIT.RL_lo LIT.RL_hi]; ...
              '|RC|@f',M.RC,[LIT.RC_thr 1.0]; ...
              'S',M.S,[LIT.S_lo LIT.S_hi]; ...
              'Penetration',M.pen,[LIT.pen_lo LIT.pen_hi]; ...
              'Phase align',M.align,[0.4 1.0]};
    letters={'(a)','(b)','(c)','(d)','(e)'};
    for m=1:size(metrics,1)
        ax=subplot(1,size(metrics,1),m);
        vals=metrics{m,2}; lit=metrics{m,3};
        b=bar(1:n_cyc,vals,0.55,'FaceColor','flat','EdgeColor','k');
        for i=1:n_cyc
            if ~all(isfinite(lit))
                b.CData(i,:)=[0.55 0.60 0.65];        % no band defined -> neutral
            elseif isfinite(vals(i))&&vals(i)>=lit(1)&&vals(i)<=lit(2)
                b.CData(i,:)=[0.4 0.7 0.4];
            else
                b.CData(i,:)=[0.85 0.45 0.40];
            end
        end
        hold on; xl=[0.5 n_cyc+0.5];
        if all(isfinite(lit))
            patch([xl fliplr(xl)],[lit(1) lit(1) lit(2) lit(2)],C.litband,'EdgeColor','none','FaceAlpha',0.35);
            yline(lit(1),'k--','LineWidth',0.9); yline(lit(2),'k--','LineWidth',0.9);
        end
        % 'omitnan' already handles lit(2)=NaN, so the axis limit falls back to the data
        ym=max([vals lit(2)],[],'omitnan')*1.25; if ~isfinite(ym)||ym<=0,ym=1;end
        ylim([0 ym]); set(ax,'XTick',1:n_cyc,'XTickLabel',{D.name}); ax.XTickLabelRotation=35;
        title(sprintf('%s %s',letters{m},metrics{m,1}),'FontSize',12);
        for i=1:n_cyc
            if isfinite(vals(i))
                text(i,vals(i)*1.02,sprintf('%.2f',vals(i)),'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
            end
        end
        ylabel(metrics{m,1}); xlabel('Storm'); beautify();
    end
    savepro('FigS10_Cross_Storm_Summary.png'); fprintf(' Saved Fig 8\n');
    fprintf('\n[FIG 8 - SCORECARD]\n%-8s | %-7s | %-7s | %-7s | %-7s | %-7s\n','Storm','R/L','RC@f','S','Pen','Align');
    for i=1:n_cyc
        fprintf('%-8s | %-7.2f | %-7.3f | %-7.2f | %-7.3f | %-7.2f\n', ...
            D(i).name,M.RL(i),M.RC(i),M.S(i),M.pen(i),M.align(i));
    end
end
function rc=compute_RC_at_f(d)
    rc=NaN;
    if isempty(d.u_surf)||isempty(d.v_surf), return; end
    u=fillmissing(d.u_surf,'linear'); v=fillmissing(d.v_surf,'linear');
    ok=isfinite(u)&isfinite(v); u=u(ok); v=v(ok); N=numel(u); if N<12, return; end
    Omega=7.2921e-5; f_in=abs(2*Omega*sind(d.lat_mean))*86400/(2*pi);
    Fs=24/3; fn=Fs/2; flo=0.5*f_in/fn; fhi=min(1.5*f_in/fn,0.95);
    if fhi>flo&&flo>0, [b,a]=butter(2,[flo fhi],'bandpass'); u=filtfilt(b,a,u); v=filtfilt(b,a,v);
    else, u=u-movmean(u,8,'omitnan'); v=v-movmean(v,8,'omitnan'); end
    Nw=numel(u); w=hann(Nw)'; W=(u(:)+1i*v(:)).*w(:);
    Y=fftshift(fft(W)); df=Fs/Nw; fv=(-Nw/2:Nw/2-1)*df;
    P=abs(Y).^2/(sum(w.^2)*df); P(P<=0)=NaN;
    fcw=abs(fv(fv<0)); Pcw=P(fv<0); fccw=fv(fv>0); Pccw=P(fv>0);
    [fcw,ic]=sort(fcw); Pcw=Pcw(ic); [fccw,icc]=sort(fccw); Pccw=Pccw(icc);
    fc=linspace(0.05,4,250);
    Pci=interp1(fcw,Pcw,fc,'linear',0); Pcci=interp1(fccw,Pccw,fc,'linear',0);
    RC=(Pci-Pcci)./(Pci+Pcci+eps); [~,idf]=min(abs(fc-f_in)); rc=abs(RC(idf));
end
function make_fig9(D, n_cyc, OUT_DIR, C, beautify, savepro)
    fprintf('\n--- Fig 9 : Hysteresis ---\n');
    figure('Name','Fig 9','Units','normalized','Position',[0.04 0.10 0.92 0.65]);
    letters={'(a)','(b)','(c)','(d)'};
    for i=1:n_cyc
        d=D(i); k=d.keep;
        ax=subplot(1,n_cyc,i);
        if any(k)
            tau=d.stress(k); dndt=d.dnike_dt(k); t=d.time(k);
            ok=isfinite(tau)&isfinite(dndt); tau=tau(ok); dndt=dndt(ok); t=t(ok);
            if numel(tau)>=3
                hrs=hours(t-t(1));
                scatter(tau,dndt,28,hrs,'filled','MarkerEdgeColor','k'); hold on;
                plot(tau,dndt,'-','Color',[0.3 0.3 0.3],'LineWidth',0.8);
                plot(tau(1),dndt(1),'o','MarkerSize',10,'MarkerFaceColor',[0.2 0.6 0.3],'MarkerEdgeColor','k');
                plot(tau(end),dndt(end),'s','MarkerSize',10,'MarkerFaceColor',[0.85 0.2 0.2],'MarkerEdgeColor','k');
                A=0.5*sum(tau(1:end-1).*dndt(2:end)-tau(2:end).*dndt(1:end-1));
                r=corr(tau,dndt,'rows','complete');
                xline(0,'k:','HandleVisibility','off'); yline(0,'k:','HandleVisibility','off');
                xlim([min(tau)*1.05 max(tau)*1.05]); ym=max(abs(dndt))*1.15; ylim([-ym ym]);
                text(0.02,0.95,sprintf('Area=%+.2e\nCorr=%+.2f',A,r),'Units','normalized', ...
                    'FontSize',9,'FontWeight','bold','BackgroundColor','w','Margin',2,'EdgeColor',[0.5 0.5 0.5],'VerticalAlignment','top');
                cb=colorbar; cb.Label.String='hours from start';
            end
        end
        title(sprintf('%s %s (%s)',letters{i},d.name,d.basin),'FontSize',11);
        xlabel('\tau (N m^{-2})'); ylabel('dNIKE/dt (W m^{-2})'); beautify();
    end
    savepro('FigS11_Hysteresis.png'); fprintf(' Saved Fig 9\n');
    fprintf('\n[FIG 9 - HYSTERESIS]\n%-8s | %-12s | %-8s\n','Storm','LoopArea','Corr');
    for i=1:n_cyc, d=D(i);
        tau=d.stress(d.keep); dndt=d.dnike_dt(d.keep);
        ok=isfinite(tau)&isfinite(dndt); tau=tau(ok); dndt=dndt(ok);
        if numel(tau)>=3
            A=0.5*sum(tau(1:end-1).*dndt(2:end)-tau(2:end).*dndt(1:end-1));
            r=corr(tau,dndt,'rows','complete');
            fprintf('%-8s | %-12.3e | %-+8.3f\n',d.name,A,r);
        else
            fprintf('%-8s | %-12s | %-8s\n',d.name,'NA','NA');
        end
    end
end