tel = giantMagellanTelescope('resolution',51);
atm = atmosphere(photometry.V,15e-2,30);
bif = influenceFunction('monotonic',0.5);
dm = deformableMirror(51,'modes',bif,...
    'resolution',tel.resolution,...
    'validActuator',tel.pupilLogical);
[x,y] = meshgrid(linspace(-1,1,tel.resolution)*tel.R);
rho = x+1i*y;
B0 = phaseStats.covarianceToeplitzMatrix(atm,rho,'mask',tel.pupilLogical);
F = dm.modes.modes(tel.pupilLogical,:);
% F = bsxfun(@minus,F,mean(F));
A = F'*F;
B = F'*B0*F;
fprintf(' >> DM Generalized Eigen Vectors ....')
[V,D] = eig(full(A),B);
fprintf('\b\b\b\b!!\n')
FV = F*V;
d = diag(D);
[ds,ix] = sort(real(d));
FV = FV(:,ix);
% FV = FV/diag(max(FV));
%%
BV = FV'*B0*FV;
U = FV'*FV;
varKL = diag(BV)./diag(U)/dm.nValidActuator;
varZE = zernikeStats.variance(zernike(telescope(25),1:500),atm);
figure(101)
subplot(1,2,1)
loglog(1:dm.nValidActuator,varKL,'.',...
    1:500,varZE,'r+')
grid
subplot(1,2,2)
loglog(1:dm.nValidActuator,cumsum(varKL),'.',...
    1:500,cumsum(varZE),'r+')
line(get(gca,'xlim'),phaseStats.variance(atm)*ones(1,2),...
    'color','k')
grid
set(gca,'ylim',[1e2,1e3])
%%
map = zeros(tel.resolution^2,dm.nValidActuator);
map(tel.pupilLogical,:) = real(FV);
map = reshape(map,...
    [tel.resolution,tel.resolution,dm.nValidActuator]);
figure(102)
for k=1:dm.nValidActuator
    imagesc(map(:,:,k));
    axis square
    colorbar('location','northoutside')
    title(sprintf('#%4d/%4d',k,dm.nValidActuator))
    drawnow;
end

% %%
% [Mp,D2] = eig(full(A));
% D = sqrt(D2);
% M = Mp/D;
% Cp = D*Mp'*B*Mp*D;
% [Ap,S] = eig(Cp);
% KL = M*Ap;
% FVp = F*KL;
% BV = FVp'*B0*FV;
% varKL = diag(BV);
% 
% %%
% map = zeros(tel.resolution^2,dm.nValidActuator);
% map(tel.pupilLogical,:) = FVp;
% map = reshape(map,...
%     [tel.resolution,tel.resolution,dm.nValidActuator]);
% figure
% for k=1:dm.nValidActuator
%     u = dm.nValidActuator - k + 1;
%     imagesc(map(:,:,u));
%     axis square
%     colorbar('location','northoutside')
%     title(sprintf('#%4d/%4d',u,dm.nValidActuator))
%     pause;
% end
% 
