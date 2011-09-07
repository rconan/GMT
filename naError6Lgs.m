close all
clear
addpath('/home/rconan/matlab/GMT/mcode/')

%%
nPixelPerLenslet = 30;
wfsNyquistSampling = 0.5;
tel = giantMagellanTelescope('resolution',50*nPixelPerLenslet,'samplingTime',1/500);
% tel = telescope(25,'resolution',50*20,'samplingTime',1/500);
wfs = shackHartmann(50,tel.resolution/2,0.85);
wfs.lenslets.nyquistSampling = wfsNyquistSampling;
wfs.lenslets.fieldStopSize = 30;
ngs = source.*tel*wfs;
setValidLenslet(wfs)
wfs.referenceSlopes = wfs.slopes;
+wfs
wfs.camera.frameListener.Enabled = false;
figure
slopesDisplay(wfs)
wfs.slopesListener.Enabled = false;

zernP = zernike(tel,1:22);
focus = zernike(tel,4);
% zern = zern\wfs;
%%
bif = influenceFunction('monotonic',0.5);
dm = deformableMirror(wfs.lenslets.nLenslet+1,'modes',bif,'resolution',tel.resolution,'validActuator',wfs.validActuator);
ngs = ngs.*tel;
dmWfsCalib = calibration(dm,wfs,ngs,ngs.wavelength*2,100);
 dmWfsCalib.threshold = 2e5;
% %% Calibration
% calibDmCommands = speye(dm.nValidActuator)*ngs.wavelength;
% if dm.nValidActuator>1000
%     steps           = 25;
% else
%     steps = 1;
% end
% nC              = floor(dm.nValidActuator/steps);
% u               = 0;
% poke               = zeros(wfs.nSlope,dm.nValidActuator);
% ngs = ngs.*tel*dm*wfs;
% fprintf(' . actuators range:          ')
% while u(end)<dm.nValidActuator
%     u = u(end)+1:min(u(end)+nC,dm.nValidActuator);
%     fprintf('\b\b\b\b\b\b\b\b\b%4d:%4d',u(1),u(end))
%     dm.coefs = calibDmCommands(:,u);
%     +ngs;
%     poke(:,u) = wfs.slopes/ngs.wavelength;
% end
% fprintf('\n--------------------\n')
%%
bifLowRes = influenceFunction('monotonic',0.5);
dmLowRes = deformableMirror(wfs.lenslets.nLenslet+1,'modes',bifLowRes,'resolution',wfs.lenslets.nLenslet+1,'validActuator',wfs.validActuator);
F = dmLowRes.modes.modes(dmLowRes.validActuator,:);
% P = D/F;
iP = F*dmWfsCalib.M;
iP = repmat( {iP} , 1 ,6 );
iP = blkdiag(iP{:});
%%
naData = load('naUbc0.mat');
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

%% Calibration
zern = zernike(tel,4);
zern.c = ngs.wavelength/4;
ngs = ngs.*zern*wfs;
focus = focus.\wfs;
D4 = focus.c*4/ngs.wavelength;

% focus.c = focus.c/D;
% ngs = ngs.*zern*-focus*wfs;
% focus = focus.\wfs;
zernP = zernike(tel,1:27);
zernP.c = eye(zernP.nMode)*ngs.wavelength/4;
ngs = ngs.*zernP*wfs;
zernP = zernP\wfs;
Dz = zernP.c*4/ngs.wavelength;
% Dz(1,:) = [];
  
lgs = source('height',90e3,'wavelength',photometry.Na); 
% tel.focalDistance = 90e3;
lgs = lgs.*tel*wfs;
zoomResolution = 50;
deltaHeight = -3e3:zoomResolution:3e3;
n = length(deltaHeight);
D = zeros(wfs.nSlope,n);
for k=1:n
    lgs = source('height',90e3+deltaHeight(k),'wavelength',photometry.Na); 
    lgs        = lgs.*tel*wfs;
    D(:,k)     = wfs.slopes;
end

figure
imagesc(D)
colorbar

% %% Calibration test I
% lgsHeight = 92.65e3;
% lgs = source('height',lgsHeight,'wavelength',photometry.Na); 
% tel.focalDistance = 90e3;
% lgs = lgs.*tel*wfs;
% error = bsxfun(@minus,D,wfs.slopes);
% error = sqrt(sum(error.^2));
% figure
% plot(deltaHeight,error)
% [~,index] = min(error);
% fprintf('New focal: %gm\n',tel.focalDistance+deltaHeight(index))
% lgs = source('height',lgsHeight,'wavelength',photometry.Na); 
% tel.focalDistance = 90e3+deltaHeight(index);
% lgs = lgs.*tel*wfs;
% error = bsxfun(@minus,D,wfs.slopes);
% error = sqrt(sum(error.^2));
% [~,index] = min(error);
% fprintf('New focal: %gm\n',tel.focalDistance+deltaHeight(index))
% 
% %% Calibration test II
% lgs = source('height',linspace(85,95,21)*1e3,'wavelength',photometry.Na); 
% tel.focalDistance = 90e3;
% lgs = lgs.*tel*wfs;
% error = bsxfun(@minus,D,wfs.slopes);
% error = sqrt(sum(error.^2));
% [~,index] = min(error);
% figure
% plot(deltaHeight,error,'.-')
% fprintf('New focal: %gm\n',tel.focalDistance+deltaHeight(index))
% 
% % tel.focalDistance = 90e3;
% % focus.c = 0; 
% % for k=1:20
% %     lgs = lgs.*tel*-focus*wfs;
% %     focus = focus.\wfs;
% %     focus.c = focus.c/D;    
% %     fprintf('Focus #%d: %g\n',k,focus.c*1e6)
% %     pause
% % end
% 
% % focus.c = focus.c/D;
% % lgs = lgs.*tel*-focus*wfs;
% % focus = focus.\wfs;
% % fprintf('Focus: %g\n',focus.c)

%%
wfs.camera.frameListener.Enabled = false;
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

%  index = hBin >= 85 & hBin<=95;
a = 90e3:lgsAltBin:95e3;
b = -fliplr(-90e3:lgsAltBin:-85e3);
lgsHeight0 = [b(1:end-1) a];

nLgs = 6;
clear lgs
for kLgs = 1:6
    lgs{kLgs} = source('zenith',arcsec(35),'azimuth',(kLgs-1)*pi/3,...
        'height',lgsHeight0,'wavelength',photometry.Na,'viewPoint',[xL(kLgs),yL(kLgs)]);
    set(lgs{kLgs},'objectiveFocalLength',90e3)
    lgsWfs{kLgs} = shackHartmann(50,tel.resolution/2,0.85);
    lgsWfs{kLgs}.lenslets.nyquistSampling = wfsNyquistSampling;
    lgsWfs{kLgs}.lenslets.fieldStopSize = 30;
    ngs = source.*tel*lgsWfs{kLgs};
    setValidLenslet(lgsWfs{kLgs})
    lgsWfs{kLgs}.referenceSlopes = lgsWfs{kLgs}.slopes;
    +lgsWfs{kLgs};
%     lgs{kLgs} = lgs{kLgs}.*tel*lgsWfs{kLgs};
%     lgsWfs{kLgs}.referenceSlopes = lgsWfs{kLgs}.slopes + lgsWfs{kLgs}.referenceSlopes;
end
% tel.focalDistance = 90e3;

error = bsxfun(@minus,D,wfs.slopes);
error = sqrt(sum(error.^2));
[errorMin,index] = min(error);
% tel.focalDistance = 90e3+deltaHeight(index);
% fprintf('New focal: %gm\n',tel.focalDistance)

lpfN = 1;
lpfCount = 0;
lpfSlopes = zeros(wfs.nSlope,1);

slopes = zeros(wfs.nSlope,nLgs,nT);
a4 = zeros(1,nT);
hNaMean = zeros(1,nT);
naBinnedSubProfile = 1e3*naBinnedSubProfile./max(naBinnedSubProfile(:));
lgsHeight = lgsHeight0'*ones(1,nT+lpfN);
naHeight  = lgsHeight0'*ones(1,nT+lpfN);
deltaHeightTrig = 0;

figure(101)
naProfile = ones(size(lgsHeight0))*max(naBinnedSubProfile(:))';
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

% h = waitbar(0,'patience!');
for kLgs=1:nLgs
    for kT=1:1
        
%         if lpfCount==lpfN
%             
%             lgsWfs{kLgs}.slopes = lpfSlopes/lpfN;
%             
%             error = bsxfun(@minus,D,lgsWfs{kLgs}.slopes);
%             error = sqrt(sum(error.^2));
%             [~,index] = min(error);
%             deltaHeightTrig = deltaHeight(index);
%             
%             lpfCount = 0;
%             lpfSlopes = zeros(lgsWfs{kLgs}.nSlope,1);
%             
% %             line(kT-1,a4(kT-1),'Marker','o','MarkerEdgeColor','r','Parent',ax2);
%             %         disp('TRIGGERING ZOOM!')
%         end
        
        indInt = kT+lpfN;
        lgsHeight(:,indInt) = lgsHeight(:,kT) - gain*deltaHeightTrig;
        naHeight(:,indInt)  = naHeight(:,kT)  + gain*deltaHeightTrig;
        
        hNaMean(kT) = sum(hBin'.*naBinnedSubProfile(:,kT+kT0))/sum(naBinnedSubProfile(:,kT+kT0));
%         set(h3,'xdata',1:kT,'ydata',hNaMean(1:kT));
%         axTimeWindow = [max(0,kT-19),kT+1];
%         set(ax3,'xlim',axTimeWindow)
        
        % Na profile update
        naProfile = interp1(hBin,naBinnedSubProfile(:,kT+kT0),naHeight(:,kT)*1e-3);
%         lgs = source('height',lgsHeight(:,kT),'wavelength',photometry.Na,'viewPoint',[xL,yL]);
        
%         set(h1,'xdata',lgsHeight(:,kT)*1e-3,'ydata',naProfile)
%         set(h1a,'xdata',ones(1,2)*hNaMean(kT),'ydata',get(ax1,'ylim'))
        hNaMeanE = 1e-3*sum(naHeight(:,kT).*naProfile)/sum(naProfile);
%         set(h1b,'xdata',ones(1,2)*hNaMeanE,'ydata',get(ax1,'ylim'))
        
        for k=1:length(lgs{kLgs})
%             lgs{kLgs}(k).height = lgsHeight(k,kT);
            lgs{kLgs}(k).nPhoton = naProfile(k);
        end
        lgs{kLgs} = lgs{kLgs}.*tel*lgsWfs{kLgs};
        
        slopes(:,kLgs,kT) = lgsWfs{kLgs}.slopes;
        %     focus = focus.\wfs;
        %     a4(kT) = 1e9*focus.c/D4;
%         zernP = zernP.\wfs;
%         zc = Dz\zernP.c(2:end);
%         a4(kT) = sqrt(sum(zc(2:end).^2,1))*1e9;
        
%         set(h2,'xdata',1:kT,'ydata',a4(1:kT),'MarkerFaceColor','b')
%         set(ax2,'xlim',axTimeWindow)
        
%         lpfCount  = lpfCount + 1;
%         lpfSlopes = lpfSlopes + lgsWfs{kLgs}.slopes;
        
%         drawnow
    end
end
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
asm = nan(dm.nActuator^2*nSs,1);
asm(repmat(dm.validActuator(:),nSs,1)) = bsxfun(@minus,buf,mean(buf));
asm = reshape(asm,[51,51*nSs]);
figure('Name','Tomographic LGS Aberrations')
imagesc(asm*1e9)
axis xy equal tight
xlabel(colorbar('location','southOutside'),'nm')
set(gca,'xtickLabel',[],'ytickLabel',[])
title(sprintf('wavefront rms, @%d": %.2fnm - @%d": %.2fnm - @%d": %.2fnm',[(0:30:60)',std(buf)'*1e9]') )

% save('naErrorClosedLoop6LGSCentralLaunch.mat','slopes')