%% TESTING OF LGS ABERRATIONS

%% Parameters definition
telDiameter    = 8;  % [m]
launchLocation = 4;  % [m]
naHeight       = 90e3; % [m]
naWidth        = 10e3; % [m]
elongation     = skyAngle(2*naWidth*launchLocation/naHeight^2);
nLenslet       = 17;
sizeSubap      = telDiameter/nLenslet; % [m]
fwhmSubap      = skyAngle(589e-9/sizeSubap);
nPixelSubAp    = ceil(1+elongation/(0.5*fwhmSubap));
resolution     = nLenslet*nPixelSubAp;

%% Wavefront sensor
wfs = shackHartmann(nLenslet,resolution,0.85);
wfs.camera.frameListener.Enabled = true;
figure

%% Telescope
tel = telescope(8,'resolution',resolution);

%%
ngs = source.*tel*wfs;
wfs.INIT;
+wfs;

%%
lgs = source('wavelength',photometry.Na,...
    'height',linspace(85,95,31)*1e3,...
    'viewpoint',[launchLocation,0]);
set(lgs,'objectiveFocalLength',90e3)
lgs = lgs.*tel*wfs;
