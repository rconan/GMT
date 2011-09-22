function naError6LgsFastWithDitherFun(nTrange)


%%
pixelPerSubap = 15;
nPixelPerLenslet = ceil(30/pixelPerSubap)*pixelPerSubap;

tel = giantMagellanTelescope('resolution',50*nPixelPerLenslet,'samplingTime',1/500);

wfs = shackHartmann(50,pixelPerSubap*50,0.85);

wfs.lenslets.nyquistSampling = 0.5;
wfs.lenslets.fieldStopSize = nPixelPerLenslet/2;

ngs = source.*tel*wfs;

subaps = fitsread('/home/rconan/matlab/GMT/mcode/lgsAberrations/subaps.fits');
wfs.validLenslet = logical(subaps);
wfs.referenceSlopes = wfs.slopes;

+wfs;
figure
imagesc(wfs.camera)

% %%
% wfs.camera.frameListener.Enabled = true;
% figure
% slopesDisplay(wfs)
% wfs.slopesListener.Enabled = true;

% %% Centroid gain estimation
% zern = zernike(tel,2);
% fwhmPerPixel = length(wfs.lenslets.imagelets)/length(wfs.camera.frame);
% pixelSize = skyAngle(fwhmPerPixel*ngs.wavelength*50/tel.D);
% zc = 0:0.1:1;
% c = zeros(size(zc));
% nZc = length(zc);
% a = skyAngle(4*zc*1e-6/tel.D);
% a.ARCSEC
% ngs = ngs.*tel*zern*wfs;
% for kZc = 1:nZc
%     zern.c = zc(kZc)*1e-6;
%     +ngs
%     c(kZc) = pixelSize.ARCSEC*median(wfs.slopes(1:end/2));
% %     gain(kZc) = a.ARCSEC/c;
% end
% figure
% plot(a.ARCSEC,c)
% grid
% %%
% % [x,y] = meshgrid(linspace(-1,1,wfs.lenslets.nLenslet)*tel.R*(wfs.lenslets.nLenslet-1)/wfs.lenslets.nLenslet);
% % lensletsCoord = [x(wfs.validLenslet) y(wfs.validLenslet)];
% % fits_write('subapertureCoordinates.fits',lensletsCoord)
% 
% 
% %%
% % bif = influenceFunction('monotonic',0.5);
% % dm = deformableMirror(wfs.lenslets.nLenslet+1,'modes',bif,'resolution',tel.resolution,'validActuator',wfs.validActuator);
% % ngs = ngs.*tel;
% % dmWfsCalib = calibration(dm,wfs,ngs,ngs.wavelength/2,dm.nValidActuator);
% %  dmWfsCalib.threshold = 2e5;
% % save('../mat/dmData32x32','bif','dm','dmWfsCalib')
% %%
wfsId = sprintf('%dx%d',wfs.lenslets.nLensletImagePx,wfs.lenslets.nLensletImagePx);
% dmDataFile = ['../mat/dmData',wfsId];
% fprintf('Loading %s ...',dmDataFile)
% load(dmDataFile,'bif','dm','dmWfsCalib')
% fprintf('\b\b\b!\n')
% %% 
% bifLowRes = influenceFunction('monotonic',0.5);
% dmLowRes = deformableMirror(wfs.lenslets.nLenslet+1,'modes',bifLowRes,'resolution',wfs.lenslets.nLenslet+1,'validActuator',wfs.validActuator);
% F = dmLowRes.modes.modes(dmLowRes.validActuator,:);
% % P = D/F;
% iP = F*dmWfsCalib.M;
% iP = repmat( {iP} , 1 ,6 );
% iP = blkdiag(iP{:});

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

%%
nLgs = 6;

u = (0:wfs.lenslets.nLensletImagePx-1)-(wfs.lenslets.nLensletImagePx)/2;
[x,y] = meshgrid(u);
r = hypot(x,y);
f = exp(-r.^2/(5./(2.*sqrt(2*log(2)))));

lgs = source('asterism',{[6,arcsec(35),0]},...
        'height',lgsHeight0,'wavelength',photometry.Na);
for kLgs=1:nLgs
    set(lgs(1,kLgs,:),...
        'objectiveFocalLength',90e3,...
        'viewPoint',[xL(kLgs),yL(kLgs)],...
        'nPhoton',285e4/length(lgsHeight0));
end
set(lgs,'extent',f)

%%
imageletsFile = ['../mat/imageletsDithered',launchType,'Launch_50x50lenslets_',wfsId,'pixels_35na_marcos.mat'];
matObj = matfile(imageletsFile);
if matObj.Properties.Writable % check file existence
    lgs = lgs.*tel;
    zern = zernike(tel,2:3);
    tiltAngle = skyAngle(1,'mas');
    zernCoef = tiltAngle*tel.D/4;
    matObj.imagelets = zeros(wfs.lenslets.nLensletsImagePx, wfs.lenslets.nLensletsImagePx*length(lgsHeight0) , nLgs);
    matObj.imageletsTip  = zeros(wfs.lenslets.nLensletsImagePx, wfs.lenslets.nLensletsImagePx*length(lgsHeight0) , nLgs);
    matObj.imageletsTilt = zeros(wfs.lenslets.nLensletsImagePx, wfs.lenslets.nLensletsImagePx*length(lgsHeight0) , nLgs);
    for kLgs=1:nLgs
        fprintf(' >> LGS #%d/%d\n',kLgs,nLgs)
        lgs_k = lgs(1,kLgs,:);
        fprintf('    -> straight!\n')
        lgs_k = lgs_k.*tel;
        propagateThrough(wfs.lenslets,lgs_k(:))
        matObj.imagelets(:,:,kLgs) = wfs.lenslets.imagelets;
        % Tip
        fprintf('    -> tip dithered!\n')
        zern.c = [zernCoef;0];
        lgs_k = lgs_k.*tel*zern;
        propagateThrough(wfs.lenslets,lgs_k(:))
        matObj.imageletsTip(:,:,kLgs) = wfs.lenslets.imagelets;
        % Tilt
        fprintf('    -> tilt dithered!\n')
        zern.c = [0;zernCoef];
        lgs_k = lgs_k.*tel*zern;
        propagateThrough(wfs.lenslets,lgs_k(:))
        matObj.imageletsTilt(:,:,kLgs) = wfs.lenslets.imagelets;
    end
    % imageletsFile = ['../mat/imageletsDithered',launchType,'Launch_50x50lenslets_',wfsId,'pixels_35na_marcos.mat'];
    %  fprintf('Saving imagelets to %s ...',imageletsFile)
    %  save(imageletsFile,'imagelets','imageletsTip','imageletsTilt')
    %  fprintf('\b\b\b!\n')
    %fprintf('Loading imagelets to %s ...',imageletsFile)
    %load(imageletsFile,'imagelets')
    %fprintf('\b\b\b!\n')
    %%
    matObj.imagelets = reshape(matObj.imagelets,[wfs.lenslets.nLensletsImagePx,wfs.lenslets.nLensletsImagePx,35,6]);
    matObj.imageletsTip = reshape(matObj.imageletsTip,[wfs.lenslets.nLensletsImagePx,wfs.lenslets.nLensletsImagePx,35,6]);
    matObj.imageletsTilt = reshape(matObj.imageletsTilt,[wfs.lenslets.nLensletsImagePx,wfs.lenslets.nLensletsImagePx,35,6]);
    % telFlux = tel.samplingTime*lgs(1,1,1).nPhoton.*tel.area;
    % imagelets = imagelets./length(lgsHeight0);
end
%%
% lgs = lgs.*tel;


naProfile3D = zeros(1,1,length(lgsHeight0));

%%
wfs.lenslets.nArray = 6;
nT = nTrange;
slopes = zeros(wfs.nSlope,nLgs,length(nT));

% wfs.camera.photonNoise  = true;
% wfs.camera.readOutNoise = 2.5;

imagelets = matObj.imagelets;

tStart = tic;
for kT=nT
    
    fprintf(' Profile #%3d\n',kT)
    naProfile = interp1(hBin,naBinnedSubProfile(:,kT),lgsHeight0*1e-3);
    
    
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
fitsPath = '~/public_html/share/adaptiveOptics/lgsAberrations/noiseless/dither';
ccdId = sprintf('%dx%d',wfs.camera.resolution/wfs.lenslets.nLenslet);
naId  = sprintf('Na%d-%d',nT(1),nT(end));
fitsFile = fullfile(fitsPath,[lower(launchType),'LaunchWfsSlopes',ccdId,naId,'.fits']);
fits_write(fitsFile,slopes)

% %%
% wfs.lenslets.nArray = 6;
% nT = 500:700;
% slopes = zeros(wfs.nSlope,nLgs,length(nT));
% 
% % wfs.camera.photonNoise  = true;
% % wfs.camera.readOutNoise = 2.5;
% 
% imagelets = matObj.imageletsTip;
% 
% tStart = tic;
% for kT=nT
%     
%     fprintf(' Profile #%3d\n',kT)
%     naProfile = interp1(hBin,naBinnedSubProfile(:,kT),lgsHeight0*1e-3);
%     
%     
%     naProfile3D(1,1,:) = naProfile;
%     wfs.lenslets.imagelets = squeeze( sum( bsxfun( @times , imagelets, naProfile3D ) , 3) );
%     wfs.lenslets.imagelets = reshape(  wfs.lenslets.imagelets , wfs.lenslets.nLensletsImagePx , [] );
%     
%     %         lgs = lgs.*tel*wfs;
%     spotsSrcKernelConvolution(wfs,lgs)
%     grab(wfs.camera)
%     dataProcessing(wfs);
%     
%     slopes(:,:,kT-nT(1)+1) = wfs.slopes;
%     
% end
% toc(tStart)
% fitsPath = '~/public_html/share/adaptiveOptics/lgsAberrations/noiseless/dither';
% ccdId = sprintf('%dx%d',wfs.camera.resolution/wfs.lenslets.nLenslet);
% fitsFile = fullfile(fitsPath,[lower(launchType),'LaunchWfsSlopes',ccdId,'TipDither.fits']);
% fits_write(fitsFile,slopes)
% 
% %%
% wfs.lenslets.nArray = 6;
% nT = 500:700;
% slopes = zeros(wfs.nSlope,nLgs,length(nT));
% 
% % wfs.camera.photonNoise  = true;
% % wfs.camera.readOutNoise = 2.5;
% 
% imagelets = matObj.imageletsTilt;
% 
% tStart = tic;
% for kT=nT
%     
%     fprintf(' Profile #%3d\n',kT)
%     naProfile = interp1(hBin,naBinnedSubProfile(:,kT),lgsHeight0*1e-3);
%     
%     
%     naProfile3D(1,1,:) = naProfile;
%     wfs.lenslets.imagelets = squeeze( sum( bsxfun( @times , imagelets, naProfile3D ) , 3) );
%     wfs.lenslets.imagelets = reshape(  wfs.lenslets.imagelets , wfs.lenslets.nLensletsImagePx , [] );
%     
%     %         lgs = lgs.*tel*wfs;
%     spotsSrcKernelConvolution(wfs,lgs)
%     grab(wfs.camera)
%     dataProcessing(wfs);
%     
%     slopes(:,:,kT-nT(1)+1) = wfs.slopes;
%     
% end
% toc(tStart)
% fitsPath = '~/public_html/share/adaptiveOptics/lgsAberrations/noiseless/dither';
% ccdId = sprintf('%dx%d',wfs.camera.resolution/wfs.lenslets.nLenslet);
% fitsFile = fullfile(fitsPath,[lower(launchType),'LaunchWfsSlopes',ccdId,'TiltDither.fits']);
% fits_write(fitsFile,slopes)
