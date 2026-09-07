clear all
clc
%  ========================================================================
%% description
%  ========================================================================
% Displacement Vector Field Imaging, 3d vector computation 
% author: Shota Hodono
% Date: Oct 2025
% cleaned up: Sep 2026

%  ========================================================================
%% settings
%  ========================================================================
bSEPIA  = 0;
Order   = 1;
bFilt   = 1; % gaussian filter YES NO

gamma = 2*pi*42.59e6;
delta =  0.0065;
G     = 69e-3;

phi2dis =1/(2*gamma*G*delta);

%  ========================================================================
%% load data
%  ========================================================================
sPath     = './_data/';

sPathOUT  = [sPath,'_processed/'];

sData  = 'dvfi_sphe_Gz1p0_Gx1p0_TR6';
a      = dir([sPath,sData,'*.nii']);
mag    = single(niftiread([sPath,a(1).name]));
phs    = single(niftiread([sPath,a(2).name]));
I      = niftiinfo([sPath,a(1).name]);
b      = dir([sPath,sData,'*.bvec']);


% alternatively you can call txet file you made
% dvs          = load('~/Desktop/NL/work/FUS/DVFI/vSets/2trans/128.txt');

dvstmp = load([sPath,b(1).name]);
dvstmp = dvstmp';
dvstmp = squeeze(dvstmp(3:2:end,:));

% if coronoal acquisition, you may need to play around a bit here,
% depending on how you defined gradient axis
dvs(:,1) = dvstmp(:,1);
dvs(:,2) = dvstmp(:,3);
dvs(:,3) = dvstmp(:,2);

dvs=dvs*-1;
nDir         = size(dvs,1);

phs    = phs / 4096*pi;
b0     = 2; % nu,ber of b0 volumes in the beginning 

mag    = squeeze(mag(:,:,:,b0+1:end));
phs    = squeeze(phs(:,:,:,b0+1:end));

data  = mag.*exp(1i.*phs);

[nRe, nPh, nSl, nMe] = size(data);
mask     = single(niftiread([sPath,'mask']));
mask     = squeeze(mask(:,:,1:nSl));


%  ========================================================================
%% calculate differences
%  ========================================================================

C1 = data(:,:,:,1:2:end);
C2 = data(:,:,:,2:2:end);
D  = phi2dis*angle(C1.*conj(C2));

if bFilt
    for ii = 1:size(D,4)
        D(:,:,:,ii) = imgaussfilt3(squeeze(D(:,:,:,ii)),0.5);
    end
end

%  ========================================================================
%% vector calculations
%  ========================================================================
% reshape D and mask into vectors
Dmat    = reshape(D, [], size(D,4));   % (nVoxels x nDir)
maskvec = mask(:) > 0;
Dmat    = Dmat(maskvec,:);             % only masked voxels

% solve using least-squares for many right-hand sides:
% dvs is nDir x 3, Dmat' is nDir x nVoxelsMasked
% Solve for X (3 x nVoxelsMasked): dvs \ Dmat'  => (3 x m)
X = dvs \ Dmat';   % least-squares solution for each voxel

% allocate dv and place back
dv            = zeros(nRe*nPh*nSl,3);
dv(maskvec,:) = X';
dv            = reshape(dv, nRe, nPh, nSl, 3);
dvf_x         = dv(:,:,:,1);
dvf_y         = dv(:,:,:,2); 
dvf_z         = dv(:,:,:,3);
dvf_mag       = sqrt(dvf_x.^2 + dvf_y.^2 + dvf_z.^2);


%  ========================================================================
%% save
%  ========================================================================

I.Datatype        = 'double';
I.ImageSize       = [nRe nPh nSl];
I.PixelDimensions = I.PixelDimensions(1:3);

ss =[sPathOUT,sData,'_dvfi_mag'];
if bFilt
    ss = [ss,'_SM']; % smoothed
end
niftiwrite(dvf_mag,ss,I)


%  ========================================================================
%% show
%  ========================================================================
iSl = 7; 
figure
subplot(121)
imagesc(mean(abs(squeeze(data(:,:,iSl,:))),3)), axis image, colorbar
im = mean(abs(data),4);
title('mag image')
subplot(122)
imagesc(1/phi2dis.*(squeeze(mean(D(:,:,iSl,:),4))),[-.1 .1]), axis image, colorbar
title('mean difference [rad]')




x_comp     = (squeeze(dv(:,:,:,1)));
y_comp     = (squeeze(dv(:,:,:,2)));
z_comp     = (squeeze(dv(:,:,:,3)));
figure
subplot(131)
imagesc(squeeze(dvf_mag(:,:,iSl))), axis image
hold on
[Z,X] = meshgrid(1:nRe,1:nPh);
quiver(Z,X,squeeze(dvf_z(:,:,iSl)),squeeze(dvf_x(:,:,iSl)),3,'r','linewidth',1)
vz = squeeze(dvf_z(:,:,iSl));
vx = squeeze(dvf_x(:,:,iSl));
iSl=39;

subplot(132)% XY vectors
imagesc(squeeze(dvf_mag(:,iSl,:))), axis image
hold on
[X,Y] = meshgrid(1:nSl,1:nPh);
quiver(X,Y,squeeze(dvf_y(:,iSl,:)),squeeze(dvf_x(:,iSl,:)),3,'r','linewidth',1)

iSl=40;
subplot(133)
[Z,Y] = meshgrid(1:nSl,1:nRe);
imagesc(squeeze(dvf_mag(iSl,:,:))), axis image
hold on
quiver(Z,Y,squeeze(y_comp(iSl,:,:)),squeeze(z_comp(iSl,:,:)),3,'r','linewidth',1)
%% show encoding
r_sphere = 0.5;
figure; hold on; axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z'); view(3)

% --- Draw transparent sphere ---
[xs, ys, zs] = sphere(50);
surf(r_sphere*xs, r_sphere*ys, r_sphere*zs, ...
    'FaceAlpha',0.2,'EdgeColor','none','FaceColor',[0.6 0.8 1]);
% --- Draw coordinate axes ---
L = 0.8;
quiver3(0,0,0,L,0,0,'r','LineWidth',2); text(L,0,0,'x')
quiver3(0,0,0,0,L,0,'g','LineWidth',2); text(0,L,0,'y')
quiver3(0,0,0,0,0,L,'b','LineWidth',2); text(0,0,L,'z')
% --- Draw rays ---
for ii = 1:size(dvs,1)
    dir = dvs(ii,:);
    P = r_sphere * dir;
    plot3([0 P(1)], [0 P(2)], [0 P(3)], 'k-', 'LineWidth', 1);
    plot3(P(1), P(2), P(3), 'ro', 'MarkerSize', 5);
end


%% show normalized
figure
iSl = 7;
imshow(squeeze(dvf_mag(:,:,iSl)), [0 1.5e-6]), axis image
hold on
% hImg = imagesc(squeeze(dvf_mag(:,:,iSl)), [0 1.5e-6]); axis image
% hold on
% set(hImg, 'AlphaData', 0.8); % Adjust 0.5 to change transparency (0 = clear, 1 = solid)

 cmap_path = '/usr/local/fsl/pkgs/fsleyes-1.10.4-pyh31c8845_0/site-packages/fsleyes/assets/colourmaps/brain_colours/x_rain_iso.cmap';
% 2. Read the file by forcing MATLAB to treat it as text
fsl_x_rain_iso = readmatrix(cmap_path, 'FileType', 'text');
% Extract and apply the colormap matrix from the toolbox
% colormap(fsl_x_rain_iso);

[Z, X] = meshgrid(1:nRe, 1:nPh);

% Extract the 2D slice components
dz = squeeze(dvf_z(:,:,iSl));
dx = squeeze(dvf_x(:,:,iSl));

% Calculate magnitude and normalize vectors to unit length
mag = sqrt(dz.^2 + dx.^2);
% Avoid division by zero for points with 0 magnitude
mag(mag == 0) = 1; 

dz_norm = dz ./ mag;
dx_norm = dx ./ mag;

% Plot normalized quiver arrows (set auto-scaling off using 0 or adjust arrow scale)
quiver(Z, X, dz_norm, dx_norm, 0.5, 'r', 'linewidth', 1), axis off