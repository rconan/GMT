%% GMT LTAO MODELING WITH OOMAO
% Demonstrate how to build the GMT LTAO system

forceSettings = false;
filename  = '../mat/gmtNonSegSingleDmLtaoSettings.mat';
if exist(filename,'file') && ~forceSettings
    
    fprintf('   >>> LOAD SETTINGS FROM %s ....',upper(filename))
    load(filename)
    fprintf('\b\b\b\b!\n')
    
else
    %% Definition of the atmosphere
    % atm = atmosphere(photometry.V,0.15,30,...
    %     'altitude',[0,4,10]*1e3,...
    %     'fractionnalR0',[0.7,0.25,0.05],...
    %     'windSpeed',[5,10,20],...
    %     'windDirection',[0,pi/4,pi]);
    atm = gmtAtmosphere(1);
    
    %% Definition of the telescope
    nLenslet = 50;
    nPx = nLenslet*8;
    tel = telescope(25,...
        'obstructionRatio',0.4,...
        'fieldOfViewInArcMin',2.5,...
        'resolution',nPx,...
        'samplingTime',1/500);
%     tel = giantMagellanTelescope(...
%         'fieldOfViewInArcMin',2.5,...
%         'resolution',nPx,...
%         'samplingTime',1/500);
%     
    %% Definition of a calibration source
    ngs = source('wavelength',photometry.Na);
    
    %% Definition of the wavefront sensor
    wfs = shackHartmann(nLenslet,nPx,0.85);
    %%
    % Propagation of the calibration source to the WFS through the telescope
    ngs = ngs.*tel*wfs;
    wfs.INIT;
    +wfs;
    %%
    % The WFS camera display:
    figure
    subplot(1,2,1)
    imagesc(wfs.camera)
    %%
    % The WFS slopes display:
    subplot(1,2,2)
    slopesDisplay(wfs)
    
    %% Definition of the deformable mirror
    bif = influenceFunction('monotonic',50/100);
    % Cut of the influence function
    % figure
    % show(bif)
    nActuator = nLenslet + 1;
    dm = deformableMirror(nActuator,...
        'modes',bif,...
        'resolution',nPx,...
        'validActuator',wfs.validActuator);
    
    
    %% Interaction matrix: DM/WFS calibration
    ngs = ngs.*tel;
    dmWfsCalib = calibration(dm,wfs,ngs,ngs.wavelength/2);
    dmWfsCalib.threshold = 4e5;
    commandMatrix = dmWfsCalib.M;
    
    %% Tip-Tilt Sensor
    % Tip-Tilt source
    tt = source('wavelength',photometry.K);
    % GMT Tip-Tilt IR sensor
    tipTiltWfs = gmtInfraredQuadCellDetector(tel,1/tel.samplingTime,40e-3,tt);
    tipTiltWfs.camera.readOutNoise = 0;
    tipTiltWfs.camera.photonNoise = false;
    tt = tt.*tel*tipTiltWfs;
    figure
    imagesc(tipTiltWfs.camera.frame)
    %% Interaction matrix: DM/TT sensor calibration
    dmTipTiltWfsCalib = calibration(dm,tipTiltWfs,tt,tt.wavelength/4);
    dmTipTiltWfsCalib.nThresholded = 0;
    commandTipTilt = dmTipTiltWfsCalib.M;
    %%
    lgs = source('asterism',{[6,arcsec(35),0]},'wavelength',photometry.Na,'height',90e3);
    ltaoMmse = linearMMSE(dm.nActuator,tel.D,atm,lgs,ngs,'pupil',dm.validActuator,'unit',-9);
    %%
    bifLowRes = influenceFunction('monotonic',0.5);
    dmLowRes = deformableMirror(wfs.lenslets.nLenslet+1,'modes',bifLowRes,...
        'resolution',wfs.lenslets.nLenslet+1,'validActuator',wfs.validActuator);
    F = dmLowRes.modes.modes(dmLowRes.validActuator,:);
    iP0 = F*commandMatrix;
    iP = repmat( {iP0} , 1 ,6 );
    iP = blkdiag(iP{:});
    M = ltaoMmse.mmseBuilder{1}*iP;
    iF = pinv(full(F));
    M = iF*M;
    
    % Combining the atmosphere and the telescope
    tel = tel+atm;
    figure
    imagesc(tel)
    
    fprintf('   >>> SAVING SETTINGS TO %s ....',upper(filename))
    save(filename)
    fprintf('\b\b\b\b!\n')

    
end

%%
%% The closed loop
% Resetting the DM command
dm.coefs = 0;
%%
% Propagation throught the atmosphere to the telescope
ngs=ngs.*tel;
%%
% Saving the turbulence aberrated phase
turbPhase = ngs.meanRmPhase;
%%
% Propagation to the WFS
ngs=ngs*dm*wfs;
%%
% Display of turbulence and residual phase
figure(11)
rad2mic = 1e6/ngs.waveNumber;
h = imagesc([turbPhase,ngs.meanRmPhase]*rad2mic);
ax = gca;
axis equal tight
ylabel(colorbar,'WFE [\mum]')
%%
% Closed loop integrator gain:
loopGain = 0.5;
%%
% closing the loop
nIteration = 2000;
total  = zeros(1,nIteration);
residue = zeros(1,nIteration);
dm.coefs = 0;
ux = ones(1,length(lgs));
tic
for kIteration=1:nIteration
    % Propagation throught the atmosphere to the telescope, +tel means that
    % all the layers move of one step based on the sampling time and the
    % wind vectors of the layers
    ngs=ngs.*+tel;
    lgs = lgs.*tel*dm*wfs;
    % Saving the turbulence aberrated phase
    turbPhase = ngs.meanRmPhase;
    % Variance of the atmospheric wavefront
    total(kIteration) = var(ngs);
    % Propagation to the WFS
    ngs=ngs*dm;
    % Variance of the residual wavefront
    residue(kIteration) = var(ngs);
    % Computing the DM residual coefficients
    meanRmWfsSlopes = bsxfun(@minus,wfs.slopes,mean(wfs.slopes));
    %     residualDmCoefs = M*meanRmWfsSlopes(:);
    % -.- TT sensing
    tt = tt.*tel*dm*tipTiltWfs;
    residualDmTtCoefs = commandTipTilt*tipTiltWfs.slopes;
    %     % TT sensing -.-
    %     % Integrating the DM coefficients
    %     dm.coefs = dm.coefs - loopGain*(residualDmCoefs + residualDmTtCoefs);
    
    % DM slopes
    S_DM = dmWfsCalib.D*dm.coefs;
    % LGS full turbulence slopes estimate
    S_LGS = meanRmWfsSlopes-S_DM*ux;
    % on-axis full turbulence DM command tomography estimate
    C_NGS = M*S_LGS(:);
    % on-axis residual turbulence DM command estimate
    C_res_NGS = C_NGS + dm.coefs;
    % integrator or (may be) low-pass filter (to check)
    dm.coefs = dm.coefs - loopGain*(C_res_NGS+residualDmTtCoefs);
    
    % Display of turbulence and residual phase
    set(h,'Cdata',[turbPhase,ngs.meanRmPhase]*rad2mic)
    title(ax,sprintf('#%4d/%4d',kIteration,nIteration))
    drawnow
end
tic
%%
% Piston removed phase variance
u = (0:nIteration-1).*tel.samplingTime;
atm.wavelength = ngs.wavelength;
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
