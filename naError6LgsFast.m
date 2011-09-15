close all
clear
addpath('/home/rconan/matlab/GMT/mcode/gmt')

%%
pixelPerSubap = 15;
nPixelPerLenslet = ceil(30/pixelPerSubap)*pixelPerSubap;

tel = giantMagellanTelescope('resolution',50*nPixelPerLenslet,'samplingTime',1/500);

wfs = shackHartmann(50,pixelPerSubap*50,0.85);

wfs.lenslets.nyquistSampling = 0.5;
wfs.lenslets.fieldStopSize = nPixelPerLenslet;

ngs = source.*tel*wfs;

subaps = fitsread('/home/rconan/matlab/GMT/mcode/lgsAberrations/subaps.fits');
wfs.validLenslet = logical(subaps);
wfs.referenceSlopes = wfs.slopes;

+wfs;
figure
imagesc(wfs.camera)
%%
% [x,y] = meshgrid(linspace(-1,1,wfs.lenslets.nLenslet)*tel.R*(wfs.lenslets.nLenslet-1)/wfs.lenslets.nLenslet);
% lensletsCoord = [x(wfs.validLenslet) y(wfs.validLenslet)];
% fits_write('subapertureCoordinates.fits',lensletsCoord)


%%
% bif = influenceFunction('monotonic',0.5);
% dm = deformableMirror(wfs.lenslets.nLenslet+1,'modes',bif,'resolution',tel.resolution,'validActuator',wfs.validActuator);
% ngs = ngs.*tel;
% dmWfsCalib = calibration(dm,wfs,ngs,ngs.wavelength*8,dm.nValidActuator);
%  dmWfsCalib.threshold = 2e5;
% % load('../mat/dmData','bif','dm','dmWfsCalib')
wfsId = sprintf('%dx%d',wfs.lenslets.nLensletImagePx,wfs.lenslets.nLensletImagePx);
dmDataFile = ['../mat/dmData',wfsId];
fprintf('Loading %s ...',dmDataFile)
load(dmDataFile,'bif','dm','dmWfsCalib')
fprintf('\b\b\b!\n')
%% 
bifLowRes = influenceFunction('monotonic',0.5);
dmLowRes = deformableMirror(wfs.lenslets.nLenslet+1,'modes',bifLowRes,'resolution',wfs.lenslets.nLenslet+1,'validActuator',wfs.validActuator);
F = dmLowRes.modes.modes(dmLowRes.validActuator,:);
% P = D/F;
iP = F*dmWfsCalib.M;
iP = repmat( {iP} , 1 ,6 );
iP = blkdiag(iP{:});

%%
naData = load('../mat/naUbc0.mat');
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
figure
imagesc(u,hBin,naBinnedSubProfile)
axis xy
ylabel('Altitude [km]')
colorbar
figure
imagesc(u,hBin,naBinnedSubProfileSubSamp)
axis xy
ylabel('Altitude [km]')
colorbar
% naBinnedSubProfile = naBinnedSubProfileSubSamp;
% nT = 20;
% naBinnedSubProfile = repmat( naBinnedSubProfile(:,1) , 1, nT);


%%
C = tel.petalRingRadius;

% Side launch 
xL = 2*C/2;
yL = 2*sqrt(3)*C/6;
[oL,rL] = cart2pol(xL,yL);
rL = rL*1.1;
[xL,yL] = pol2cart(oL+2*(0:2)*pi/3,rL*ones(1,3));
xL = reshape([xL;xL],[],1);
yL = reshape([yL;yL],[],1);

% Central launch
% xL = zeros(6,1);
% yL = zeros(6,1);

if all(xL) || all(yL)
    launchType = 'Edge';
else
    launchType = 'Central';
end

%  index = hBin >= 85 & hBin<=95;
a = 90e3:lgsAltBin:95e3;
b = -fliplr(-90e3:lgsAltBin:-85e3);
lgsHeight0 = [b(1:end-1) a];

nLgs = 6;
% clear lgs
% for kLgs = 1:6
%     lgs{kLgs} = source('zenith',arcsec(35),'azimuth',(kLgs-1)*pi/3,...
%         'height',lgsHeight0,'wavelength',photometry.Na,'viewPoint',[xL(kLgs),yL(kLgs)]);
%     set(lgs{kLgs},'objectiveFocalLength',90e3)
%     lgsWfs{kLgs} = shackHartmann(50,tel.resolution/2,0.85);
% %     lgsWfs{kLgs}.validLenslet = logical(subaps);
%     lgsWfs{kLgs}.lenslets.nyquistSampling = wfsNyquistSampling;
%     lgsWfs{kLgs}.lenslets.fieldStopSize = 30;
%     ngs = source.*tel*lgsWfs{kLgs};
%     setValidLenslet(lgsWfs{kLgs})
%     lgsWfs{kLgs}.referenceSlopes = lgsWfs{kLgs}.slopes;
%     +lgsWfs{kLgs};
% %     lgs{kLgs} = lgs{kLgs}.*tel*lgsWfs{kLgs};
% %     lgsWfs{kLgs}.referenceSlopes = lgsWfs{kLgs}.slopes + lgsWfs{kLgs}.referenceSlopes;
% end
% tel.focalDistance = 90e3;

lgs = source('asterism',{[6,arcsec(35),0]},...
        'height',lgsHeight0,'wavelength',photometry.Na);
for kLgs=1:nLgs
    set(lgs(1,kLgs,:),...
        'objectiveFocalLength',90e3,...
        'viewPoint',[xL(kLgs),yL(kLgs)],...
        'nPhoton',285e4/length(lgsHeight0));
end

[x,y] = meshgrid((0:14)-7);
r = hypot(x,y);
f = exp(-r.^2/(5./(2.*sqrt(2*log(2)))));
set(lgs,'extent',f)

%%
% imagelets = zeros(wfs.lenslets.nLensletsImagePx, wfs.lenslets.nLensletsImagePx*length(lgsHeight0) , nLgs);
% for kLgs=1:nLgs
%      fprintf(' >> LGS #%d/%d\n',kLgs,nLgs)
%      lgs_k = lgs(1,kLgs,:);
%      lgs_k = lgs_k.*tel;
%      propagateThrough(wfs.lenslets,lgs_k(:))
%      imagelets(:,:,kLgs) = wfs.lenslets.imagelets;
% end
imageletsFile = ['../mat/imagelets',launchType,'Launch_50x50lenslets_',wfsId,'pixels_35na_marcos.mat'];
% fprintf('Saving imagelets to %s ...',imageletsFile)
% save(imageletsFile,'imagelets')
% fprintf('\b\b\b!\n')
fprintf('Loading imagelets to %s ...',imageletsFile)
load(imageletsFile,'imagelets')
fprintf('\b\b\b!\n')
%%
imagelets = reshape(imagelets,[wfs.lenslets.nLensletsImagePx,wfs.lenslets.nLensletsImagePx,35,6]);
% telFlux = tel.samplingTime*lgs(1,1,1).nPhoton.*tel.area;
% imagelets = imagelets./length(lgsHeight0);
%%
% lgs = lgs.*tel;

lpfN = 1;
lpfCount = 0;
lpfSlopes = zeros(wfs.nSlope,1);

slopes = zeros(wfs.nSlope,nLgs,nT);
a4 = zeros(1,nT);
hNaMean = zeros(1,nT);
naBinnedSubProfile = naBinnedSubProfile./max(naBinnedSubProfile(:));
lgsHeight = lgsHeight0'*ones(1,nT+lpfN);
naHeight  = lgsHeight0'*ones(1,nT+lpfN);
deltaHeightTrig = 0;

figure(101)
naProfile = ones(size(lgsHeight0));
subplot(2,2,3)
h1 = plot(lgsHeight0*1e-3,naProfile,'r.-');
ax1 = get(h1,'Parent');
h1a = line(ones(1,2)*90,get(ax1,'ylim'),'color','k','Parent',ax1);
h1b = line(ones(1,2)*90,get(ax1,'ylim'),'color','m','Parent',ax1);
grid
xlabel('Height [km]')
ylabel('Na Profile')
subplot(2,2,2)
h2 = plot(0,0,'.-');
grid
ax2 = get(h2,'Parent');
xlabel('Time [s]')
ylabel('Focus [nm]')
subplot(2,2,1)
h3 = plot(0,90,'.-');
grid
ax3 = get(h3,'Parent');
xlabel('Time [s]')
ylabel('Mean Na height [km]')
drawnow

kT0 = 0;%770;
gain = 0.0;
focus.c = 0;

naProfile3D = zeros(1,1,length(lgsHeight0));

%%
% h = waitbar(0,'patience!');
wfs.lenslets.nArray = 6;
nT = 500:700;
slopes = zeros(wfs.nSlope,nLgs,length(nT));

wfs.camera.photonNoise  = true;
wfs.camera.readOutNoise = 2.5;

tStart = tic;
for kT=nT
    
    fprintf(' Profile #%3d\n',kT)
    naProfile = interp1(hBin,naBinnedSubProfile(:,kT+kT0),lgsHeight0*1e-3);
    
    set(h1,'xdata',lgsHeight0*1e-3,'ydata',naProfile)
    
    naProfile3D(1,1,:) = naProfile;
    wfs.lenslets.imagelets = squeeze( sum( bsxfun( @times , imagelets, naProfile3D ) , 3) );
    wfs.lenslets.imagelets = reshape(  wfs.lenslets.imagelets , wfs.lenslets.nLensletsImagePx , [] );
    
    %         lgs = lgs.*tel*wfs;
    spotsSrcKernelConvolution(wfs,lgs)
    grab(wfs.camera)
    dataProcessing(wfs);
    
    slopes(:,:,kT-nT(1)+1) = wfs.slopes;
    
end
toc(tStart)
%% INTENSITY CENTER
% centroid = zeros(wfs.nSlope,nLgs);
% for kLgs=1:nLgs
%     obj = lgsWfs{kLgs};
%     [nPx,mPx,nFrame]  = size(obj.camera.frame);
%     nLensletArray = obj.lenslets.nArray;
%     nPxLenslet = nPx/obj.lenslets.nLenslet;
%     mPxLenslet = mPx/obj.lenslets.nLenslet/nLensletArray;
%     indexRasterLenslet ...
%         = utilities.rearrange([nPx,mPx/nLensletArray,nLensletArray*nFrame],[nPxLenslet,mPxLenslet]);
%     % remove index from non-valid lenslets
%     v = ~obj.validLenslet(:);
%     v = repmat(v,nLensletArray,1);
%     v = repmat(v,nFrame,1);
%     indexRasterLenslet(:,v) = [];
%     buffer     = obj.camera.frame(indexRasterLenslet);
%     centroid(:,kLgs) = intensityCenter(buffer);
% end
% cogSlopes = slopes;
% slopes = centroid;
%%
S = zeros(50,50*6);
s2 = hypot(slopes(1:end/2,:,1),slopes(1+end/2:end,:,1));
S(repmat(wfs.validLenslet,1,6)) = s2;

figure('Name','Slopes map')
imagesc(S)
axis xy equal tight
colorbar('location','northOutside')
%%
Sx = zeros(50,50*6);
Sy = zeros(50,50*6);
p = repmat(wfs.validLenslet,1,6);
Sx(p) = slopes(1:end/2,:,1);
Sy(p) = slopes(1+end/2:end,:,1);
figure('Name','X/Y slopes map')
imagesc([Sx;Sy])
axis xy equal tight
colorbar('location','northOutside')
%%
wft = iP*reshape(slopes(:,:,1),[],1); 
buf = reshape(wft,[dm.nValidActuator,6]);
buf = squeeze(std(buf,[],1));
fprintf('pupil rms mean %4.2fnm and std %4.2fnm\n',...
    [1e9*mean(buf,2),1e9*std(buf,[],2)]')
maps = zeros(dm.nActuator^2*6,1);
maps(repmat(dm.validActuator(:),6,1),:) = wft;
maps = reshape(maps,[dm.nActuator,dm.nActuator*6]);
figure('Name','LGS Aberrations')
imagesc(maps*1e9)
axis xy equal tight
xlabel(colorbar('location','northOutside'),'nm')
drawnow
%%
atm = gmtAtmosphere(1);
gs = source('asterism',{[6,arcsec(30),0]},'wavelength',photometry.Na,'magnitude',6,'height',90e3);
ss = source('zenith',arcsec(0:30:60),'azimuth',zeros(1,3),'wavelength',photometry.K,'magnitude',6);
ltaoMmse = linearMMSE(dm.nActuator,tel.D,atm,gs,ss,'pupil',dm.validActuator,'unit',-9);
%%
nSs = length(ss);
buf = cell2mat( cellfun(@(x) x*wft,ltaoMmse.mmseBuilder,'UniformOutput',false)');
fprintf('asm rms mean %4.2fnm and std %4.2fnm\n',...
    [mean(std(buf,[],1))*1e9,std(std(buf,[],1))*1e9]')
asm = zeros(dm.nActuator^2*nSs,1);
asm(repmat(dm.validActuator(:),nSs,1)) = bsxfun(@minus,buf,mean(buf));
asm = reshape(asm,[51,51*nSs]);
figure('Name','Tomographic LGS Aberrations')
imagesc(asm*1e9)
axis xy equal tight
xlabel(colorbar('location','southOutside'),'nm')
set(gca,'xtickLabel',[],'ytickLabel',[])
title(sprintf('wavefront rms, @%d": %.2fnm - @%d": %.2fnm - @%d": %.2fnm',[(0:30:60)',std(buf)'*1e9]') )
%%
%  save('/home/rconan/matlab/GMT/mcode/mat/naErrorClosedLoop6LgsEdgeLaunch4x4_na500-700_marcos.mat','slopes')
%%
fitsPath = '~/public_html/share/adaptiveOptics/lgsAberrations/noise/noThreshold/';
ccdId = sprintf('%dx%d',wfs.camera.resolution/wfs.lenslets.nLenslet);
fitsFile = fullfile(fitsPath,[lower(launchType),'LaunchWfsSlopes',ccdId,'.fits']);
fits_write(fitsFile,slopes(:,:,1)*589e-9/0.5*constants.radian2arcsec*2)
% fits_write('wfsPhase15x15_conv_na580.fits',maps)
% fits_write('wfsTomoPhase15x15_conv_na580.fits',asm)
