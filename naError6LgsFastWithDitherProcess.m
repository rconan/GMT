%%
telLowRes = giantMagellanTelescope('resolution',50,'samplingTime',1/500);
nZernMode = 28;
zernTTPFree = zernike(telLowRes,5:nZernMode,'unitNorm',true);
%%
atm = gmtAtmosphere(1);
gs = source('asterism',{[6,arcsec(30),0]},'wavelength',photometry.Na,'magnitude',6,'height',90e3);
ss = source('zenith',arcsec(0),'azimuth',zeros(1,1),'wavelength',photometry.K,'magnitude',6);
ltaoMmse = linearMMSE(telLowRes.resolution,telLowRes.D,atm,gs,ss,'pupil',telLowRes.pupilLogical,'unit',-9);

%%
pixelPerSubap = 15;
nPixelPerLenslet = ceil(30/pixelPerSubap)*pixelPerSubap;
tel = giantMagellanTelescope('resolution',50*nPixelPerLenslet,'samplingTime',1/500);

wfs15 = shackHartmann(50,pixelPerSubap*50,0.85);

wfs15.lenslets.nyquistSampling = 0.5;
wfs15.lenslets.fieldStopSize = nPixelPerLenslet/2;

u = (0:wfs15.lenslets.nLensletImagePx-1)-(wfs15.lenslets.nLensletImagePx)/2;
[x,y] = meshgrid(u);
r = hypot(x,y);
f = exp(-r.^2/(5./(2.*sqrt(2*log(2)))));
ngs15 = source;
ngs15.extent = f;

ngs15 = ngs15.*tel*wfs15;

subaps = fitsread('/home/rconan/matlab/GMT/mcode/lgsAberrations/subaps.fits');
wfs15.validLenslet = logical(subaps);
wfs15.referenceSlopes = wfs15.slopes;
% wfs15.INIT

zern = zernike(tel,2:nZernMode,'unitNorm',true);
stroke = ngs15.wavelength/2;
zern.c = eye(zern.nMode)*stroke;
ngs15 = ngs15.*zern*wfs15;
D15 = wfs15.slopes/stroke;
figure,imagesc(wfs15.slopes)

%%
pixelPerSubap = 10;
nPixelPerLenslet = ceil(30/pixelPerSubap)*pixelPerSubap;
tel = giantMagellanTelescope('resolution',50*nPixelPerLenslet,'samplingTime',1/500);

wfs10 = shackHartmann(50,pixelPerSubap*50,0.85);

wfs10.lenslets.nyquistSampling = 0.5;
wfs10.lenslets.fieldStopSize = nPixelPerLenslet/2;

u = (0:wfs10.lenslets.nLensletImagePx-1)-(wfs10.lenslets.nLensletImagePx)/2;
[x,y] = meshgrid(u);
r = hypot(x,y);
f = exp(-r.^2/(5./(2.*sqrt(2*log(2)))));
ngs10 = source;
ngs10.extent = f;

ngs10 = ngs10.*tel*wfs10;

subaps = fitsread('/home/rconan/matlab/GMT/mcode/lgsAberrations/subaps.fits');
wfs10.validLenslet = logical(subaps);
wfs10.referenceSlopes = wfs10.slopes;

zern = zernike(tel,2:nZernMode,'unitNorm',true);
stroke = ngs10.wavelength;
zern.c = eye(zern.nMode)*stroke;
ngs10 = ngs10.*zern*wfs10;
D10 = wfs10.slopes/stroke;
figure,imagesc(wfs10.slopes)

%%
pixelPerSubap = 4;
nPixelPerLenslet = ceil(30/pixelPerSubap)*pixelPerSubap;
tel = giantMagellanTelescope('resolution',50*nPixelPerLenslet,'samplingTime',1/500);

wfs04 = shackHartmann(50,pixelPerSubap*50,0.85);

wfs04.lenslets.nyquistSampling = 0.5;
wfs04.lenslets.fieldStopSize = nPixelPerLenslet/2;

u = (0:wfs04.lenslets.nLensletImagePx-1)-(wfs04.lenslets.nLensletImagePx)/2;
[x,y] = meshgrid(u);
r = hypot(x,y);
f = exp(-r.^2/(5./(2.*sqrt(2*log(2)))));
ngs04 = source;
ngs04.extent = f;

ngs04 = ngs04.*tel*wfs04;

subaps = fitsread('/home/rconan/matlab/GMT/mcode/lgsAberrations/subaps.fits');
wfs04.validLenslet = logical(subaps);
wfs04.referenceSlopes = wfs04.slopes;

zern = zernike(tel,2:nZernMode,'unitNorm',true);
stroke = ngs04.wavelength;
zern.c = eye(zern.nMode)*stroke;
ngs04 = ngs04.*zern*wfs04;
D04 = wfs04.slopes/stroke;
figure,imagesc(wfs04.slopes)
%%
fitsPath = '/home/rconan/public_html/share/adaptiveOptics/lgsAberrations/noiseless/dither';
ngs15 = ngs15.*zernTTPFree;
ngs10 = ngs10.*zernTTPFree;
ngs04 = ngs04.*zernTTPFree;

slopes15 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes15x15Na1-1000.fits'));
slopes10 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes10x10Na1-1000.fits'));
slopes04 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes4x4Na1-1000.fits'));
zernC = zeros(nZernMode-1,1000,6);
opds = cell(3,6);
h = waitbar(0,'Patience...');
for k=1:6
    
    waitbar(k/6)
    
    zernC(:,:,k) = D15\squeeze(slopes15(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs15;
    opds{1,k} = ngs15.opd(:,:,1)*1e9;
    
    zernC(:,:,k) = D10\squeeze(slopes10(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs10;
    opds{2,k} = ngs10.opd(:,:,1)*1e9;
    
    zernC(:,:,k) = D04\squeeze(slopes04(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs04;
    opds{3,k} = ngs04.opd(:,:,1)*1e9;
%     imagesc(ngs.opd(:,:,1)*1e9);
%     title(sprintf('%d',k))
%     pause
end
close(h)

%%
tomoOpd = zeros(telLowRes.resolution);
count = 1;
figure
for k=1:3
    
    buf = ltaoMmse.mmseBuilder{1}*cell2mat( cellfun( @(x) x(telLowRes.pupilLogical) ,  opds(k,:)' , 'uniformOutput' , false ) );
    tomoOpd(telLowRes.pupilLogical) = buf;

    subplot(3,1,count);count=count+1;
    imagesc( [ cell2mat(opds(k,:)) tomoOpd] )
    axis equal tight,ylabel(colorbar('location','westOutside'),'wfe [nm]')
    title(sprintf('Tomo. opd rms = %4.2fnm',std(buf)))
    
%     subplot(3,2,count);count=count+1;
%     imagesc(tomoOpd),axis equal tight,ylabel(colorbar('location','eastOutside'),'wfe [nm]')
%     set(gca,'visible','off')

end
%%
fitsPath = '/home/rconan/public_html/share/adaptiveOptics/lgsAberrations/noise/noThreshold';
ngs15 = ngs15.*zernTTPFree;
ngs10 = ngs10.*zernTTPFree;
ngs04 = ngs04.*zernTTPFree;

slopes04 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes4x4Na1-1000.fits'));
slopes15 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes15x15Na1-1000.fits'));
slopes10 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes10x10Na1-1000.fits'));
fitsPath = '/home/rconan/public_html/share/adaptiveOptics/lgsAberrations/noise/';
slopes15t = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes15x15Na1-1000.fits'));
slopes10t = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes10x10Na1-1000.fits'));

zernC = zeros(27,1000,6);
opdsNoise = cell(3,6);
h = waitbar(0,'Patience...');
for k=1:6
    
    waitbar(k/6)
    
    zernC(:,:,k) = D15\squeeze(slopes15t(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs15;
    opds{1,k} = ngs15.opd(:,:,1)*1e9;
    
    zernC(:,:,k) = D10\squeeze(slopes10t(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs10;
    opds{2,k} = ngs10.opd(:,:,1)*1e9;
    
    zernC(:,:,k) = D04\squeeze(slopes04(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs04;
    opds{3,k} = ngs04.opd(:,:,1)*1e9;

end
close(h)

%%
tomoOpd = zeros(telLowRes.resolution);
count = 1;
figure
for k=1:3
    
    buf = ltaoMmse.mmseBuilder{1}*cell2mat( cellfun( @(x) x(telLowRes.pupilLogical) ,  opds(k,:)' , 'uniformOutput' , false ) );
    tomoOpd(telLowRes.pupilLogical) = buf;

    subplot(3,1,count);count=count+1;
    imagesc( [ cell2mat(opds(k,:)) tomoOpd] )
    axis equal tight,ylabel(colorbar('location','westOutside'),'wfe [nm]')
    title(sprintf('Tomo. opd rms = %4.2fnm',std(buf)))
    
%     subplot(3,2,count);count=count+1;
%     imagesc(tomoOpd),axis equal tight,ylabel(colorbar('location','eastOutside'),'wfe [nm]')
%     set(gca,'visible','off')

end

%%
% tomoOpdNoise = zeros(telLowRes.resolution);
% buf = ltaoMmse.mmseBuilder{1}*cell2mat( cellfun( @(x) x(telLowRes.pupilLogical) ,  opdsNoise(2,:)' , 'uniformOutput' , false ) );
% std(buf)
% tomoOpdNoise(telLowRes.pupilLogical) = buf;
% figure,imagesc([tomoOpd,tomoOpdNoise]);axis equal tight; colorbar

%%
fitsPath = '/home/rconan/public_html/share/adaptiveOptics/lgsAberrations/noise/';
slopes = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes15x15Na1-1000.fits'));
wfs10.slopes = squeeze(slopes(:,1,:))*zernCoef*2*1e9/2;
zern = zern.\wfs10;
zc1515 = zern.c';
mean(sqrt(sum(zc1515(:,5:end).^2,2)))
%%
fitsPath = '/home/rconan/public_html/share/adaptiveOptics/lgsAberrations/noise/';
slopes = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes10x10Na1-1000.fits'));
wfs10.slopes = squeeze(slopes(:,1,:))*zernCoef*2*1e9/2;
zern = zern.\wfs10;
zc1010 = zern.c';
mean(sqrt(sum(zc1010(:,5:end).^2,2)))
%%
fitsPath = '/home/rconan/public_html/share/adaptiveOptics/lgsAberrations/noise/';
slopes = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes15x15Na1-1000radThresh.fits'));
wfs10.slopes = squeeze(slopes(:,1,:))*zernCoef*2*1e9/2;
zern = zern.\wfs10;
zc1515 = zern.c';
mean(sqrt(sum(zc1515(:,5:end).^2,2)))

