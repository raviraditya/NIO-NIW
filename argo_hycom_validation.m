% Argo-HYCOM co-location and evaluation. Produces Table S1 and Figure S2,
% and the report that the supplement figure script reads.

clc; clear; close all; warning('off','all');

DATA_DIR = 'path/to/data';        % folder with HYCOM and ARGO subfolders

%% config: folder holding the daily *_prof.nc files per storm
ARGODIR = struct( ...
  'KYARR',  fullfile(DATA_DIR,'ARGO','Kyarr'), ...
  'AMPHAN', fullfile(DATA_DIR,'ARGO','Amp'), ...
  'FANI',   fullfile(DATA_DIR,'ARGO','Fani'), ...
  'TAUKTAE',fullfile(DATA_DIR,'ARGO','Tauk'));

TSDIR = struct( ...
  'KYARR',  fullfile(DATA_DIR,'HYCOM','Kyarr_TS3Z'), ...
  'AMPHAN', fullfile(DATA_DIR,'HYCOM','Amphan_TS3Z'), ...
  'FANI',   fullfile(DATA_DIR,'HYCOM','Fani_TS3Z'), ...
  'TAUKTAE',fullfile(DATA_DIR,'HYCOM','Tauktae_TS3Z'));

STORMS = {'KYARR','AMPHAN','FANI','TAUKTAE'};
OUTDIR = fullfile(DATA_DIR,'HYCOM','Validation');
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

MAX_KM=50; MAX_HR=3; REF_Z=5; DT_MLD=0.2;

fid_all = fopen(fullfile(OUTDIR,'argo_hycom_report.txt'),'w');

for si = 1:numel(STORMS)
    storm = STORMS{si};
    adir = ARGODIR.(storm);
    tsdir = TSDIR.(storm);
    fprintf('\n=== %s (individual) ===\n', storm);

    %% gather daily prof files
    pf = dir(fullfile(adir,'*_prof.nc'));
    if isempty(pf), warning('No *_prof.nc in %s', adir); continue; end

    % concatenate all profiles across days
    juld=[]; lat=[]; lon=[]; Pc={}; Tc={}; Sc={};
    for k=1:numel(pf)
        fn = fullfile(pf(k).folder, pf(k).name);
        try
            j = ncread(fn,'JULD'); la=ncread(fn,'LATITUDE'); lo=ncread(fn,'LONGITUDE');
            P=ncread(fn,'PRES'); T=ncread(fn,'TEMP'); S=ncread(fn,'PSAL');
        catch ME
            warning('skip %s : %s', pf(k).name, ME.message); continue;
        end
        try
            Pa=ncread(fn,'PRES_ADJUSTED'); Ta=ncread(fn,'TEMP_ADJUSTED'); Sa=ncread(fn,'PSAL_ADJUSTED');
            if ~all(isnan(Pa),'all')
                P(~isnan(Pa))=Pa(~isnan(Pa)); T(~isnan(Ta))=Ta(~isnan(Ta)); S(~isnan(Sa))=Sa(~isnan(Sa));
            end
        catch, end
        np = numel(j);
        for ip=1:np
            juld(end+1,1)=j(ip); lat(end+1,1)=la(ip); lon(end+1,1)=lo(ip); %#ok<AGROW>
            Pc{end+1}=P(:,ip); Tc{end+1}=T(:,ip); Sc{end+1}=S(:,ip); %#ok<AGROW>
        end
    end
    argo_time = datetime(1950,1,1)+days(juld);
    nprof = numel(juld);
    fprintf('  loaded %d profiles from %d daily files\n', nprof, numel(pf));
    fprintf('  Argo time range: %s to %s\n', datestr(min(argo_time)), datestr(max(argo_time)));

    %% index HYCOM TS3Z
    [tsf,tst,tsi] = index_ts_dir(tsdir);
    if isempty(tsf), warning('No TS3Z in %s', tsdir); continue; end
    fn0=tsf{1}; lon_g=ncread(fn0,'lon'); lat_g=ncread(fn0,'lat'); depth=ncread(fn0,'depth');
    idz=find(depth<=700); zc=depth(idz); nz=numel(zc);
    hmin=min(tst); hmax=max(tst);

    %% match
    rows=[];
    for ip=1:nprof
        tp=argo_time(ip); la=lat(ip); lo=lon(ip);
        if isnat(tp)||~isfinite(la)||~isfinite(lo), continue; end
        if tp<hmin-hours(MAX_HR)||tp>hmax+hours(MAX_HR), continue; end
        [dth,ki]=min(abs(tst-tp)); if dth>hours(MAX_HR), continue; end
        fn=tsf{ki}; tt=tsi(ki);
        [~,ilon]=min(abs(lon_g-lo)); [~,ilat]=min(abs(lat_g-la));
        if haversine(la,lo,lat_g(ilat),lon_g(ilon))>MAX_KM, continue; end

        pa=Pc{ip}; ta=Tc{ip}; sa=Sc{ip};
        good=isfinite(pa)&isfinite(ta)&isfinite(sa);
        if sum(good)<5, continue; end
        pa=pa(good); ta=ta(good); sa=sa(good);
        [pa,o]=sort(pa); ta=ta(o); sa=sa(o);
        [pa,iu]=unique(pa,'stable'); ta=ta(iu); sa=sa(iu);
        if numel(pa)<5 || min(pa)>15, continue; end
        T5a=interp1(pa,ta,REF_Z,'linear',NaN); S5a=interp1(pa,sa,REF_Z,'linear',NaN);
        MLDa=mld_temp(pa,ta,REF_Z,DT_MLD);

        try
            th=squeeze(ncread(fn,'water_temp',[ilon ilat idz(1) tt],[1 1 nz 1]));
            sh=squeeze(ncread(fn,'salinity',  [ilon ilat idz(1) tt],[1 1 nz 1]));
        catch, continue; end
        gh=isfinite(th)&isfinite(sh); if sum(gh)<5, continue; end
        zh=zc(gh); th=th(gh); sh=sh(gh);
        [zh,iu]=unique(zh,'stable'); th=th(iu); sh=sh(iu);
        T5h=interp1(zh,th,REF_Z,'linear',NaN); S5h=interp1(zh,sh,REF_Z,'linear',NaN);
        MLDh=mld_temp(zh,th,REF_Z,DT_MLD);

        if ~isfinite(T5a)||~isfinite(T5h), continue; end
        % reject nonphysical surface values (tropical NIO: T 10-35 C, S 25-40 psu)
        if T5a<10||T5a>35||T5h<10||T5h>35, continue; end
        if isfinite(S5a) && (S5a<25||S5a>40), S5a=NaN; end
        if isfinite(S5h) && (S5h<25||S5h>40), S5h=NaN; end

        % --- deep-profile correlation over 50-500 m (T and S) ---
        zg = (50:10:500)';
        Ta_i = interp1(pa,ta,zg,'linear',NaN);   Sa_i = interp1(pa,sa,zg,'linear',NaN);
        Th_i = interp1(zh,th,zg,'linear',NaN);    Sh_i = interp1(zh,sh,zg,'linear',NaN);
        rT_deep = pair_corr(Ta_i,Th_i);
        rS_deep = pair_corr(Sa_i,Sh_i);

        rows=[rows; T5a T5h S5a S5h MLDa MLDh la lo rT_deep rS_deep]; %#ok<AGROW>
    end
    n=size(rows,1);
    fprintf('  matched %d profiles\n', n);

    %% write block
    fprintf(fid_all,'HYCOM vs ARGO FULL VALIDATION REPORT\n');
    fprintf(fid_all,'Cyclone %s\n\n', storm);
    fprintf(fid_all,'PER-PROFILE COLOCATION TABLE\n');
    fprintf(fid_all,'%s\n', repmat('-',1,110));
    fprintf(fid_all,'%-4s %-8s %-9s %-8s %-9s %-5s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
        'idx','lat','lon','dt_h','dist','qc','T5a','T5h','S5a','S5h','MLDa','MLDh');
    for r=1:n
        fprintf(fid_all,'%-4d %-8.3f %-9.3f %-8.2f %-9.2f %-5d %-8.3f %-8.3f %-8.3f %-8.3f %-8.2f %-8.2f\n', ...
            r, rows(r,7),rows(r,8),0,0, 1, rows(r,1),rows(r,2),rows(r,3),rows(r,4),rows(r,5),rows(r,6));
    end
    fprintf(fid_all,'%s\n', repmat('-',1,110));
    fprintf(fid_all,'OUTPUT FILES\n\n');

    %% stats
    if n>0
        okT=isfinite(rows(:,1))&isfinite(rows(:,2));
        okS=isfinite(rows(:,3))&isfinite(rows(:,4));
        okM=isfinite(rows(:,5))&isfinite(rows(:,6));
        Tb=mean(rows(okT,2)-rows(okT,1)); Tr=sqrt(mean((rows(okT,2)-rows(okT,1)).^2));
        if sum(okT)>2, ccT=corrcoef(rows(okT,1),rows(okT,2)); TR=ccT(1,2); else, TR=NaN; end
        Sb=mean(rows(okS,4)-rows(okS,3)); Sr=sqrt(mean((rows(okS,4)-rows(okS,3)).^2));
        if sum(okS)>2, ccS=corrcoef(rows(okS,3),rows(okS,4)); SR=ccS(1,2); else, SR=NaN; end
        Mb=mean(rows(okM,6)-rows(okM,5)); Mr_rmse=sqrt(mean((rows(okM,6)-rows(okM,5)).^2));
        if sum(okM)>2, cc=corrcoef(rows(okM,5),rows(okM,6)); Mr=cc(1,2); else, Mr=NaN; end
        % deep 50-500 m correlations (columns 9,10)
        rTd = mean(rows(isfinite(rows(:,9)),9),'omitnan');
        rSd = mean(rows(isfinite(rows(:,10)),10),'omitnan');
        fprintf('  T: bias %+.2f RMSE %.2f R %.2f | S: bias %+.2f RMSE %.2f R %.2f | MLD: bias %+.1f R %.2f | deep rT %.2f rS %.2f\n', ...
            Tb,Tr,TR, Sb,Sr,SR, Mb,Mr, rTd,rSd);
    else
    end
    save(fullfile(OUTDIR,[storm '_coloc_indiv.mat']),'rows');
end

fprintf('\n=== DONE (individual) ===\n');
fprintf('Compare TableS1_stats_indiv.txt against TableS1_stats.txt (merged).\n');
function [files,times,tidx]=index_ts_dir(tsdir)
    files={}; times=datetime.empty(0,1); tidx=[];
    d=dir(fullfile(tsdir,'*.nc'));
    for i=1:numel(d)
        fn=fullfile(d(i).folder,d(i).name);
        try, tv=ncread(fn,'time'); catch, continue; end
        info=ncinfo(fn,'time'); un='';
        for a=1:numel(info.Attributes), if strcmpi(info.Attributes(a).Name,'units'), un=info.Attributes(a).Value; end, end
        tok=regexp(un,'hours since\s+([\d\-]+)[ T]([\d:]+)','tokens','once');
        if isempty(tok)
            tok=regexp(un,'hours since\s+([\d\-]+)','tokens','once'); if isempty(tok), continue; end
            t0=datetime(tok{1},'InputFormat','yyyy-MM-dd');
        else
            t0=datetime([tok{1} ' ' tok{2}],'InputFormat','yyyy-MM-dd HH:mm:ss');
        end
        tt=t0+hours(double(tv(:))); n=numel(tt);
        times=[times;tt]; files=[files;repmat({fn},n,1)]; tidx=[tidx;(1:n)'];
    end
    [times,o]=sort(times); files=files(o); tidx=tidx(o);
end
function km=haversine(lat1,lon1,lat2,lon2)
    R=6371; dlat=deg2rad(lat2-lat1); dlon=deg2rad(lon2-lon1);
    a=sin(dlat/2).^2+cosd(lat1).*cosd(lat2).*sin(dlon/2).^2;
    km=R*2*atan2(sqrt(a),sqrt(1-a));
end
function mld=mld_temp(z,T,zref,dT)
    mld=NaN; z=z(:); T=T(:); [z,o]=sort(z); T=T(o);
    [z,iu]=unique(z,'stable'); T=T(iu);
    if numel(z)<3, return; end
    Tref=interp1(z,T,zref,'linear',NaN); if ~isfinite(Tref), return; end
    below=z>=zref; zz=z(below); TT=T(below);
    idx=find(abs(TT-Tref)>=dT,1,'first');
    if isempty(idx), mld=max(zz); return; end
    if idx==1, mld=zz(1); return; end
    d1=abs(TT(idx-1)-Tref); d2=abs(TT(idx)-Tref);
    if d2==d1, mld=zz(idx); else, mld=zz(idx-1)+(zz(idx)-zz(idx-1))*(dT-d1)/(d2-d1); end
end

function r = pair_corr(a,b)
    % correlation of two profiles over their common finite depths
    r = NaN;
    ok = isfinite(a) & isfinite(b);
    if sum(ok) >= 4
        c = corrcoef(a(ok), b(ok));
        r = c(1,2);
    end
end