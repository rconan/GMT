function lgsAberrationsTest2(nTrange)

%% TESTING OF LGS ABERRATIONS

% keck_pupil = fitsread('keck_pupil.fits');
% keck_subaps = fitsread('keck_subaps.fits');

%% Parameters definition
telDiameter    = 10;  % [m]
launchLocation = 6.2;  % [m]
naHeight       = 90e3; % [m]
naWidth        = 10e3; % [m]
elongation     = skyAngle(naWidth*(launchLocation+telDiameter/2)/naHeight^2);
nLenslet       = 20;
sizeSubap      = telDiameter/nLenslet; % [m]
fwhmSubap      = skyAngle(589e-9/sizeSubap);
pixelSize      = skyAngle(2.4,'arcsec');
nPixelSubAp    = ceil(2*pixelSize/fwhmSubap);
resolution     = nLenslet*nPixelSubAp;

%% Wavefront sensor
wfs = shackHartmann(nLenslet,nLenslet*2,0.85);
wfs.quadCell = true;
wfs.centroiding = false;
wfs.camera.frameListener.Enabled = true;
wfs.lenslets.nyquistSampling = 0.5;
wfs.lenslets.fieldStopSize   = nPixelSubAp;
figure

%% Telescope
tel = telescope(telDiameter,'resolution',resolution,'samplingTime',1/500);
%%
[x,y] = meshgrid(linspace(-1,1,wfs.lenslets.nLenslet)*tel.R*(wfs.lenslets.nLenslet-1)/wfs.lenslets.nLenslet);
lensletsCoord = [x(wfs.validLenslet) y(wfs.validLenslet)];
fits_write('subapertureCoordinates.fits',lensletsCoord)

%%
ngs = source;

u = (0:wfs.lenslets.nLensletImagePx-1)-(wfs.lenslets.nLensletImagePx)/2;
[x,y] = meshgrid(u);
r = hypot(x,y);
f = exp(-r.^2/(3./(2.*sqrt(2*log(2)))));
ngs.extent = f;

ngs = ngs.*tel*wfs;
wfs.INIT;
+wfs;

% %%
% bif = influenceFunction('monotonic',0.5);
% dm = deformableMirror(nLenslet+1,'modes',bif,...
%     'resolution',tel.resolution,'validActuator',wfs.validActuator);
% ngs = ngs.*tel;
% dmWfsCalib = calibration(dm,wfs,ngs,ngs.wavelength/2);
% dmWfsCalib.threshold = 2e5;
% 
% %%
lgs = source('wavelength',photometry.Na,...
    'height',linspace(85,95,27)*1e3,...
    'viewpoint',[0,launchLocation]);
set(lgs,'objectiveFocalLength',90e3,'extent',f)
% dmWfsCalib = calibration(dm,wfs,lgs,lgs(1).wavelength/4,dm.nValidActuator);
% 
% wfs.camera.frameListener.Enabled = false;
% % dmWfsCalib = calibration(dm,wfs,lgs,lgs(1).wavelength);

%%
naData = load('../mat/naUbc0.mat');
[nH,nT] = size(naData.naUbc0);
resT = 1;  %[s]
resH = 24; %[m]

hBinFactor = 16;
naBinnedSubProfile = utilities.binning(naData.naUbc0,[nH/hBinFactor,1000]);
lgsAltBin = resH*hBinFactor;
hBin = interp1(1:nH,naData.haug,0.5*(hBinFactor+1) + (0:nH/hBinFactor-1)*hBinFactor);
u = 1:nT;
figure
imagesc(u,hBin,naBinnedSubProfile)
axis xy
ylabel('Altitude [km]')
colorbar

a = 90e3:lgsAltBin:95e3;
b = -fliplr(-90e3:lgsAltBin:-85e3);
lgsHeight0 = [b(1:end-1) a];
naBinnedSubProfile = 1e4*naBinnedSubProfile./max(naBinnedSubProfile(:));
nLgsHeight0 = length(lgsHeight0);

%%
nLgs = 1;
% imagelets = zeros(wfs.lenslets.nLensletsImagePx,...
%     wfs.lenslets.nLensletsImagePx*nLgsHeight0 , nLgs);
% for kLgs=1:nLgs
%      fprintf(' >> LGS #%d/%d\n',kLgs,nLgs)
%      lgs_k = lgs(1,kLgs,:);
%      lgs_k = lgs_k.*tel;
%      propagateThrough(wfs.lenslets,lgs_k(:))
%      imagelets(:,:,kLgs) = wfs.lenslets.imagelets;
% end
% imagelets = reshape(imagelets,[wfs.lenslets.nLensletsImagePx,...
%     wfs.lenslets.nLensletsImagePx,nLgsHeight0]);
%%
wfsId = sprintf('%dx%d',wfs.lenslets.nLensletImagePx,wfs.lenslets.nLensletImagePx);
imageletsFile = ['../mat/keckImageletsDithered_20x20lenslets_',wfsId,'pixels_',sprintf('%d',nLgsHeight0),'na_marcos.mat'];
matObj = matfile(imageletsFile);
if matObj.Properties.Writable % check file existence
    lgs = lgs.*tel;
    zern = zernike(tel,2:3);
    tiltAngle = skyAngle(1,'mas');
    zernCoef = tiltAngle*tel.D/4;
%         fprintf(' >> LGS #%d/%d\n',kLgs,nLgs)
        lgs_k = lgs;
        fprintf('    -> straight!\n')
        lgs_k = lgs_k.*tel;
        propagateThrough(wfs.lenslets,lgs_k(:))
        matObj.imagelets = wfs.lenslets.imagelets;
        % Tip
        fprintf('    -> tip dithered!\n')
        zern.c = [zernCoef;0];
        lgs_k = lgs_k.*tel*zern;
        propagateThrough(wfs.lenslets,lgs_k(:))
        matObj.imageletsTip = wfs.lenslets.imagelets;
        % Tilt
        fprintf('    -> tilt dithered!\n')
        zern.c = [0;zernCoef];
        lgs_k = lgs_k.*tel*zern;
        propagateThrough(wfs.lenslets,lgs_k(:))
        matObj.imageletsTilt = wfs.lenslets.imagelets;
    % imageletsFile = ['../mat/imageletsDithered',launchType,'Launch_50x50lenslets_',wfsId,'pixels_35na_marcos.mat'];
    %  fprintf('Saving imagelets to %s ...',imageletsFile)
    %  save(imageletsFile,'imagelets','imageletsTip','imageletsTilt')
    %  fprintf('\b\b\b!\n')
    %fprintf('Loading imagelets to %s ...',imageletsFile)
    %load(imageletsFile,'imagelets')
    %fprintf('\b\b\b!\n')
    %%
    matObj.imagelets = reshape(matObj.imagelets,[wfs.lenslets.nLensletsImagePx,wfs.lenslets.nLensletsImagePx,nLgsHeight0]);
    matObj.imageletsTip = reshape(matObj.imageletsTip,[wfs.lenslets.nLensletsImagePx,wfs.lenslets.nLensletsImagePx,nLgsHeight0]);
    matObj.imageletsTilt = reshape(matObj.imageletsTilt,[wfs.lenslets.nLensletsImagePx,wfs.lenslets.nLensletsImagePx,nLgsHeight0]);
    % telFlux = tel.samplingTime*lgs(1,1,1).nPhoton.*tel.area;
    % imagelets = imagelets./nLgsHeight0;
end

%%
wfs.lenslets.nArray = 1;
nT = nTrange;
slopes = zeros(wfs.nSlope,length(nT));

imagelets = matObj.imagelets;
parfor kT =1:length(nT)
    slopes(:,kT) = computeLgsSlopes(wfs,lgs,hBin,naBinnedSubProfile,lgsHeight0,imagelets,nT(kT));
end
fitsPath = '~/public_html/share/adaptiveOptics/lgsAberrations/noiseless/dither';
ccdId = sprintf('%dx%d',wfs.camera.resolution/wfs.lenslets.nLenslet);
naId  = sprintf('Na%d-%d',nT(1),nT(end));
fitsFile = fullfile(fitsPath,['keckLike_WfsSlopes',ccdId,naId,'.fits']);
fits_write(fitsFile,slopes)

imagelets = matObj.imageletsTip;
parfor kT =1:length(nT)
    slopes(:,kT) = computeLgsSlopes(wfs,lgs,hBin,naBinnedSubProfile,lgsHeight0,imagelets,nT(kT));
end
fitsPath = '~/public_html/share/adaptiveOptics/lgsAberrations/noiseless/dither';
ccdId = sprintf('%dx%d',wfs.camera.resolution/wfs.lenslets.nLenslet);
naId  = sprintf('Na%d-%d',nT(1),nT(end));
fitsFile = fullfile(fitsPath,['keckLike_WfsSlopes',ccdId,naId,'TipDither.fits']);
fits_write(fitsFile,slopes)

imagelets = matObj.imageletsTilt;
parfor kT =1:length(nT)
    slopes(:,kT) = computeLgsSlopes(wfs,lgs,hBin,naBinnedSubProfile,lgsHeight0,imagelets,nT(kT));
end
fitsPath = '~/public_html/share/adaptiveOptics/lgsAberrations/noiseless/dither';
ccdId = sprintf('%dx%d',wfs.camera.resolution/wfs.lenslets.nLenslet);
naId  = sprintf('Na%d-%d',nT(1),nT(end));
fitsFile = fullfile(fitsPath,['keckLike_WfsSlopes',ccdId,naId,'TiltDither.fits']);
fits_write(fitsFile,slopes)

end

function slopes = computeLgsSlopes(wfs,lgs,hBin,naBinnedSubProfile,lgsHeight0,imagelets,kT)

% slopes = zeros(wfs.nSlope,length(nT));
naProfile3D = zeros(1,1,length(lgsHeight0));

% tStart = tic;
% for kT=nT
%     
%     waitbar((kT-nT(1)+1)/length(nT))
%     fprintf(' Profile #%3d\n',kT)
    naProfile = interp1(hBin,naBinnedSubProfile(:,kT),lgsHeight0*1e-3);
    
%     set(h1,'xdata',lgsHeight0*1e-3,'ydata',naProfile)
    
    naProfile3D(1,1,:) = naProfile;
    wfs.lenslets.imagelets = squeeze( sum( bsxfun( @times , imagelets, naProfile3D ) , 3) );
    wfs.lenslets.imagelets = reshape(  wfs.lenslets.imagelets , wfs.lenslets.nLensletsImagePx , [] );
    
    %         lgs = lgs.*tel*wfs;
    spotsSrcKernelConvolution(wfs,lgs)
    grab(wfs.camera)
    dataProcessing(wfs);
    
%     set(lgs,{'nPhoton'},num2cell(naProfile'))
%     lgs = lgs.*tel*wfs;
    
    slopes = wfs.slopes;
    
% end
% toc(tStart)

end
