matPath = '/priv/meggs2/rconan/mat';
naData = load(fullfile(matPath,'naUbc0.mat'));
[nH,nT] = size(naData.naUbc0);
resT = 1;  %[s]
resH = 24; %[m]

% naAltBin = 0.5e3;
% wfsNaFov = 10e3;
% wfsNaMeanFocus = 90e3;
% lgsAlt = linspace(wfsNaMeanFocus-(naAltBin+wfsNaFov)/2,wfsNaMeanFocus+(naAltBin+wfsNaFov)/2,840);
% lgsAltBin = wfsNaMeanFocus-wfsNaFov/2:naAltBin:wfsNaMeanFocus+wfsNaFov/2;
% nT = size(naData.naUbc0,2);
% ;
% [x,y] = meshgrid(u,naData.haug);
% naSubProfile = interp2(x,y,naData.naUbc0,u,lgsAlt'*1e-3);
% naBinnedSubProfile = utilities.binning(naSubProfile,[length(lgsAltBin),nT]);
% naBinnedSubProfile = 500*50^2*1e3*naBinnedSubProfile./max(naBinnedSubProfile(:));

hBinFactor = 12;
naBinnedSubProfile = utilities.binning(naData.naUbc0,[nH/hBinFactor,1000]);
lgsAltBin = resH*hBinFactor;
hBin = interp1(1:nH,naData.haug,0.5*(hBinFactor+1) + (0:nH/hBinFactor-1)*hBinFactor);
u = 1:nT;
timeSubSamp = (0:nT-1)/10 + 770 - (1000/10)*75/100;
[x,y] = meshgrid( u , hBin);
naBinnedSubProfileSubSamp = interp2(x,y,naBinnedSubProfile,timeSubSamp,hBin');
figure(314)
subplot(4,1,1)
imagesc(u,hBin,naBinnedSubProfile)
axis xy
ylabel('Altitude [km]')
xlabel('Time [s]')
%%
telLowRes = giantMagellanTelescope('resolution',50,'samplingTime',1/500);
nZernMode = 28;
zernTTPFree = zernike(telLowRes,5:nZernMode,'unitNorm',true);
%%
atm = gmtAtmosphere(1);
gs = source('asterism',{[6,arcsec(30),0]},'wavelength',photometry.Na,'magnitude',6,'height',90e3);
ss = source('zenith',arcsec(60),'azimuth',zeros(1,1),'wavelength',photometry.K,'magnitude',6);
ltaoMmse = linearMMSE(telLowRes.resolution,telLowRes.D,atm,gs,ss,'pupil',telLowRes.pupilLogical,'unit',-9);
ltaoMmse.rms
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
zernC15 = zeros(nZernMode-1,1000,6);
zernC10 = zeros(nZernMode-1,1000,6);
zernC04 = zeros(nZernMode-1,1000,6);
opds = cell(3,6);
h = waitbar(0,'Patience...');
for k=1:6
    
    waitbar(k/6)
    
    zernC15(:,:,k) = D15\squeeze(slopes15(:,k,:));
    zernTTPFree.c = zernC15(4:end,:,k);
    +ngs15;
    opds{1,k} = ngs15.opdVector*1e9;
    
    zernC10(:,:,k) = D10\squeeze(slopes10(:,k,:));
    zernTTPFree.c = zernC10(4:end,:,k);
    +ngs10;
    opds{2,k} = ngs10.opdVector*1e9;
    
    zernC04(:,:,k) = D04\squeeze(slopes04(:,k,:));
    zernTTPFree.c = zernC04(4:end,:,k);
    +ngs04;
    opds{3,k} = ngs04.opdVector*1e9;
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
    
    buf = ltaoMmse.mmseBuilder{1}*cell2mat( opds(k,:)'  );
    tomoOpd(telLowRes.pupilLogical) = buf(:,1);
    stdTomoOpd(:,k) = std(buf)';

    subplot(3,1,count);count=count+1;
    map = zeros(telLowRes.resolution^2,6);
    map( telLowRes.pupilLogical,: ) = cell2mat( cellfun( @(x) x(:,1) , opds(k,:), 'uniformOutput', false));
    map = reshape(map,telLowRes.resolution,[]);
    imagesc( [ map tomoOpd] )
    axis equal tight,ylabel(colorbar('location','westOutside'),'wfe [nm]')
    title(sprintf('Tomo. opd rms = %4.2fnm',std(buf(:,1))))
    
%     subplot(3,2,count);count=count+1;
%     imagesc(tomoOpd),axis equal tight,ylabel(colorbar('location','eastOutside'),'wfe [nm]')
%     set(gca,'visible','off')

end
figure(314)
subplot(4,1,2)
plot(stdTomoOpd)
grid
ylabel('Altitude [km]')
xlabel('Time [s]')
title('Noiseless detector!')
legend('E2V 800','ANDOR Neo','CCD 220',0)
%%
fitsPath = '/home/rconan/public_html/share/adaptiveOptics/lgsAberrations/noise/noThreshold';
ngs15 = ngs15.*zernTTPFree;
ngs10 = ngs10.*zernTTPFree;
ngs04 = ngs04.*zernTTPFree;

slopes04 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes4x4Na1-1000.fits'));
slopes15 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes15x15Na1-1000.fits'));
slopes10 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes10x10Na1-1000.fits'));

zernC = zeros(27,1000,6);
opdsNoise = cell(3,6);
h = waitbar(0,'Patience...');
for k=1:6
    
    waitbar(k/6)
    
    zernC(:,:,k) = D15\squeeze(slopes15(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs15;
    opds{1,k} = ngs15.opdVector*1e9;
    
    zernC(:,:,k) = D10\squeeze(slopes10(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs10;
    opds{2,k} = ngs10.opdVector*1e9;
    
    zernC(:,:,k) = D04\squeeze(slopes04(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs04;
    opds{3,k} = ngs04.opdVector*1e9;

end
close(h)

%%
tomoOpd = zeros(telLowRes.resolution);
count = 1;
figure
% for k=1:3
%     
%     buf = ltaoMmse.mmseBuilder{1}*cell2mat( cellfun( @(x) x(telLowRes.pupilLogical) ,  opds(k,:)' , 'uniformOutput' , false ) );
%     tomoOpd(telLowRes.pupilLogical) = buf;
% 
%     subplot(3,1,count);count=count+1;
%     imagesc( [ cell2mat(opds(k,:)) tomoOpd] )
%     axis equal tight,ylabel(colorbar('location','westOutside'),'wfe [nm]')
%     title(sprintf('Tomo. opd rms = %4.2fnm',std(buf)))
%     
% %     subplot(3,2,count);count=count+1;
% %     imagesc(tomoOpd),axis equal tight,ylabel(colorbar('location','eastOutside'),'wfe [nm]')
% %     set(gca,'visible','off')
% 
% end
for k=1:3
    
    buf = ltaoMmse.mmseBuilder{1}*cell2mat( opds(k,:)'  );
    tomoOpd(telLowRes.pupilLogical) = buf(:,1);
    stdTomoOpd(:,k) = std(buf)';

    subplot(3,1,count);count=count+1;
    map = zeros(telLowRes.resolution^2,6);
    map( telLowRes.pupilLogical,: ) = cell2mat( cellfun( @(x) x(:,1) , opds(k,:), 'uniformOutput', false));
    map = reshape(map,telLowRes.resolution,[]);
    imagesc( [ map tomoOpd] )
    axis equal tight,ylabel(colorbar('location','westOutside'),'wfe [nm]')
    title(sprintf('Tomo. opd rms = %4.2fnm',std(buf(:,1))))
    
%     subplot(3,2,count);count=count+1;
%     imagesc(tomoOpd),axis equal tight,ylabel(colorbar('location','eastOutside'),'wfe [nm]')
%     set(gca,'visible','off')

end
e2v800Noise   = sqrt((pi^2/3)*(2.5^2/(450*0.9/0.448)^2)*(15^4/2^2))*589/(2*pi);
andorNeoNoise = sqrt((pi^2/3)*(1.4^2/(450*0.6/0.448)^2)*(10^4/2^2))*589/(2*pi);
figure(314)
subplot(4,1,3)
plot(stdTomoOpd)
h1 = line(get(gca,'xlim'),ones(1,2)*e2v800Noise,'color','b');
h2 = line(get(gca,'xlim'),ones(1,2)*andorNeoNoise,'color',[0 0.5 0]);
grid
ylabel('wfe [nm]')
xlabel('Time [s]')
title('Noisy detector!')
legend('E2V 800','ANDOR Neo','CCD 220',0)

%%
ngs15 = ngs15.*zernTTPFree;
ngs10 = ngs10.*zernTTPFree;
ngs04 = ngs04.*zernTTPFree;

fitsPath = '/home/rconan/public_html/share/adaptiveOptics/lgsAberrations/noise/';
slopes15 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes15x15Na1-1000.fits'));
slopes10 = fitsread(fullfile(fitsPath,'edgeLaunchWfsSlopes10x10Na1-1000.fits'));

zernC = zeros(27,1000,6);
opdsNoise = cell(3,6);
h = waitbar(0,'Patience...');
for k=1:6
    
    waitbar(k/6)
    
    zernC(:,:,k) = D15\squeeze(slopes15(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs15;
    opds{1,k} = ngs15.opdVector*1e9;
    
    zernC(:,:,k) = D10\squeeze(slopes10(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs10;
    opds{2,k} = ngs10.opdVector*1e9;
    
    zernC(:,:,k) = D04\squeeze(slopes04(:,k,:));
    zernTTPFree.c = zernC(4:end,:,k);
    +ngs04;
    opds{3,k} = ngs04.opdVector*1e9;

end
close(h)

%%
tomoOpd = zeros(telLowRes.resolution);
count = 1;
figure
% for k=1:3
%     
%     buf = ltaoMmse.mmseBuilder{1}*cell2mat( cellfun( @(x) x(telLowRes.pupilLogical) ,  opds(k,:)' , 'uniformOutput' , false ) );
%     tomoOpd(telLowRes.pupilLogical) = buf;
% 
%     subplot(3,1,count);count=count+1;
%     imagesc( [ cell2mat(opds(k,:)) tomoOpd] )
%     axis equal tight,ylabel(colorbar('location','westOutside'),'wfe [nm]')
%     title(sprintf('Tomo. opd rms = %4.2fnm',std(buf)))
%     
% %     subplot(3,2,count);count=count+1;
% %     imagesc(tomoOpd),axis equal tight,ylabel(colorbar('location','eastOutside'),'wfe [nm]')
% %     set(gca,'visible','off')
% 
% end
for k=1:3
    
    buf = ltaoMmse.mmseBuilder{1}*cell2mat( opds(k,:)'  );
    tomoOpd(telLowRes.pupilLogical) = buf(:,1);
    stdTomoOpd(:,k) = std(buf)';

    subplot(3,1,count);count=count+1;
    map = zeros(telLowRes.resolution^2,6);
    map( telLowRes.pupilLogical,: ) = cell2mat( cellfun( @(x) x(:,1) , opds(k,:), 'uniformOutput', false));
    map = reshape(map,telLowRes.resolution,[]);
    imagesc( [ map tomoOpd] )
    axis equal tight,ylabel(colorbar('location','westOutside'),'wfe [nm]')
    title(sprintf('Tomo. opd rms = %4.2fnm',std(buf(:,1))))
    
%     subplot(3,2,count);count=count+1;
%     imagesc(tomoOpd),axis equal tight,ylabel(colorbar('location','eastOutside'),'wfe [nm]')
%     set(gca,'visible','off')

end
figure(314)
subplot(4,1,4)
plot(stdTomoOpd)
h1 = line(get(gca,'xlim'),ones(1,2)*e2v800Noise,'color','b');
h2 = line(get(gca,'xlim'),ones(1,2)*andorNeoNoise,'color',[0 0.5 0]);
grid
ylabel('Altitude [km]')
xlabel('Time [s]')
title('Noisy thresholded detector!')
legend('E2V 800','ANDOR Neo','CCD 220',0)


