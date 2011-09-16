%% ADAPTIVE OPTICS MODELING WITH OOMAO
% Demonstrate how to build a simple closed-loop single conjugated adaptive
% optics system

%% Definition of the atmosphere 
% atm = atmosphere(photometry.V,0.15,30,...
%     'altitude',[0,4,10]*1e3,...
%     'fractionnalR0',[0.7,0.25,0.05],...
%     'windSpeed',[5,10,20],...
%     'windDirection',[0,pi/4,pi]);
atm = gmtAtmosphere(1);

%% Definition of the telescope
nPx = 300;
% tel = telescope(25,...
%     'fieldOfViewInArcMin',2.5,...
%     'resolution',nPx,...
%     'samplingTime',1/500);
tel = giantMagellanTelescope(...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/500);

%%
gmtPups = zeros([size(tel.pupil),7]);
for k=1:7;gmtPups(:,:,k)=tel.segment{k}.pupil;end
gmtPups = reshape(gmtPups,[],7);
sumGmtPups = sum(gmtPups)';

%% Definition of a calibration source
ngs = source('wavelength',photometry.Na);

%% Definition of the wavefront sensor
nLenslet = 50;
wfs = shackHartmann(nLenslet,nPx,0.75);
ngs = ngs.*tel*wfs;
wfs.INIT;
+wfs;
% The WFS camera display:
figure
subplot(1,2,1)
imagesc(wfs.camera)
% The WFS slopes display:
subplot(1,2,2)
slopesDisplay(wfs)

%% Tip-Tilt Sensor
% Tip-Tilt source
tt = source('wavelength',photometry.K);
% GMT Tip-Tilt IR sensor
tipTiltWfs = gmtInfraredQuadCellDetector(tel,1/tel.samplingTime,40e-3,tt);
tipTiltWfs.camera.readOutNoise = 0;
tipTiltWfs.camera.photonNoise = false;
tt = tt.*tel*tipTiltWfs;
figure
imagesc(tipTiltWfs.camera)

%% Definition of the deformable mirror
bif = influenceFunction('monotonic',50/100);
% % Cut of the influence function
% figure
% show(bif)
nActuator = nLenslet + 1;
tic
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',wfs.validActuator);
toc
%%% Interaction matrix: DM/WFS calibration
ngs = ngs.*tel;
dmWfsCalib = calibration(dm,wfs,ngs,ngs.wavelength/2);
dmWfsCalib.threshold = 3e5;
commandMatrix = dmWfsCalib.M;
%%% Interaction matrix: DM/TT sensor calibration
dmTipTiltWfsCalib = calibration(dm,tipTiltWfs,tt,tt.wavelength/4);
dmTipTiltWfsCalib.nThresholded = 0;
commandTipTilt = dmTipTiltWfsCalib.M;

%% Adaptive Secondary Mirror
% bifAsm = influenceFunction('monotonic',50/100);
% zActuator = gmtAsm(0.3);
% zActuator(1:72) = [];
% bifAsm.actuatorCoord = zActuator;
% tic
% asm = deformableMirror(length(zActuator),...
%     'modes',bifAsm,...
%     'resolution',nPx,...
%     'validActuator',true(length(zActuator),1));
% toc
% %%% Interaction matrix: ASM/WFS calibration
% ngs = ngs.*tel;
% asmWfsCalib = calibration(asm,wfs,ngs,ngs.wavelength/2);
% asmWfsCalib.threshold = 10.^6;
% asmWfsCalib.nThresholded = length(asmWfsCalib.eigenValues)-2000+1;
% commandMatrix = asmWfsCalib.M;
% %%% Interaction matrix: ASM/TT sensor calibration
% asmTipTiltWfsCalib = calibration(asm,tipTiltWfs,tt,tt.wavelength/4);
% asmTipTiltWfsCalib.nThresholded = 0;
% commandTipTilt = asmTipTiltWfsCalib.M;

%% The closed loop
% Combining the atmosphere and the telescope
% tel = tel+atm;
% figure
% imagesc(tel)
%%
% Resetting the DM command
dm.coefs = 0;
%%
% Propagation throught the atmosphere to the telescope
ngs=ngs.*tel;
tt = tt.*tel;
%%
% Saving the turbulence aberrated phase
turbPhase = ngs.meanRmPhase;
%%
% Propagation to the WFS
ngs=ngs*dm*wfs;
%%
% Display of turbulence and residual phase
rad2mic = 1e6/ngs.waveNumber;
rad2nm  = 1e9/ngs.waveNumber;
figure(11)
subplot(1,2,1)
h1 = imagesc(turbPhase*rad2mic);
axis equal tight
ylabel(colorbar,'\mum')
subplot(1,2,2)
h2 = imagesc(ngs.meanRmPhase*rad2nm);
axis equal tight
ylabel(colorbar,'nm')
%%
% Closed loop integrator gain:
loopGain = 0.5;
%%
% closing the loop
nIteration = 500;
srcorma = sourceorama('gmtSingleDmGmt.h5',ngs,nIteration*tel.samplingTime,tel);
total  = zeros(1,nIteration);
residue = zeros(1,nIteration);
diffPistonCoefs = 0;
pistonActuator = {tel.pupil,zeros(tel.resolution)};
zern = zernike(tel,1);
tic
for kIteration=1:nIteration
    % Propagation throught the atmosphere to the telescope, +tel means that
    % all the layers move of one step based on the sampling time and the
    % wind vectors of the layers
%     ngs=ngs.*+tel; 
%     ngs.resetPhase = data(:,:,kIteration);
    +srcorma;
    tt.resetPhase = ngs.opd*tt.waveNumber;
%     srcorma.srcs=srcorma.srcs.*tel;
    % Saving the turbulence aberrated phase
%     ngs = ngs*-(zern.\ngs);
    turbPhase = ngs.phase;
    % Variance of the atmospheric wavefront
    total(kIteration) = var(ngs);
    % Propagation to the WFS
    ngs=ngs*pistonActuator*dm*wfs; 
    % Variance of the residual wavefront
    residue(kIteration) = var(ngs);
    % Computing the DM residual coefficients
    residualDmCoefs = commandMatrix*( wfs.slopes - mean(wfs.slopes) );
    % -.- TT sensing 
%     tt = tt.*srcorma.tel*dm*tipTiltWfs;
    tt = tt*dm*tipTiltWfs;
    residualDmTtCoefs = commandTipTilt*tipTiltWfs.slopes;
    % TT sensing -.-
    % Integrating the DM coefficients
    dm.coefs = dm.coefs - loopGain*(residualDmCoefs + residualDmTtCoefs);
    % piston circus
    diffPistonCoefs = diffPistonCoefs - loopGain*(gmtPups'*ngs.phase(:))./sumGmtPups;
    diffPiston = gmtPups*diffPistonCoefs;
    pistonActuator{2} = reshape(diffPiston,tel.resolution,tel.resolution);
    % Display of turbulence and residual phase
    set(h1,'Cdata',turbPhase*rad2mic)
    set(h2,'Cdata',ngs.meanRmPhase*rad2nm)
    drawnow
end
toc
clear srcorma
u = (0:nIteration-1).*tel.samplingTime;
atm.wavelength = ngs.wavelength;
%%
% Piston removed phase variance
totalTheory = phaseStats.zernikeResidualVariance(1,atm,tel);
atm.wavelength = photometry.V;
%%
% Phase variance to micron rms converter 
rmsMicron = @(x) 1e6*sqrt(x).*ngs.wavelength/2/pi;
figure(12)
plot(u,rmsMicron(total),u([1,end]),rmsMicron(totalTheory)*ones(1,2),u,rmsMicron(residue))
grid
legend('Full','Full (theory)','Residue',0)
xlabel('Time [s]')
ylabel('Wavefront rms [\mum]')
