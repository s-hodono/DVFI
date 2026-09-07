%% Full Visualization: Cone Sampling Inside Sphere
clear all

%% setting
bSAVE = false; % saving?
sPath = '/savingpath/';

axis_dir      = [1, 0, 0];    % [x,y,z]
num_samples   = 64;           % number of encofing
theta_max_deg = 30;           % half-angle of cone in degrees

sName  = ['3D_cone_',num2str(num_samples),'x2_Dir_x',num2str(axis_dir(1)),'_y',num2str(axis_dir(2)),'_z',num2str(axis_dir(3)),'_theta_',num2str(theta_max_deg),'.dvs'];
sNameT = [sName(1:end-3),'txt'];
%%
axis_dir = axis_dir / norm(axis_dir);  % normalize

% --- Parameters ---
r_sphere = 0.5;          % sphere radius (diameter = 1)
theta_max = deg2rad(theta_max_deg);


figure; hold on; axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z'); view(3)

% --- Draw transparent sphere ---
[xs, ys, zs] = sphere(50);
surf(r_sphere*xs, r_sphere*ys, r_sphere*zs, ...
    'FaceAlpha',0.2,'EdgeColor','none','FaceColor',[0.6 0.8 1]);

% --- Draw coordinate axes ---
L = 0.8;
quiver3(0,0,0,L,0,0,'r','LineWidth',2); text(L,0,0,'X')
quiver3(0,0,0,0,L,0,'g','LineWidth',2); text(0,L,0,'Y')
quiver3(0,0,0,0,0,L,'b','LineWidth',2); text(0,0,L,'Z')

%% --- Rotation matrix to align +Z with axis_dir ---
z_axis = [0;0;1];
v = cross(z_axis, axis_dir(:));
s = norm(v);
c = dot(z_axis, axis_dir);
if s < 1e-8
    R = eye(3); % already aligned
else
    vx = [   0   -v(3)  v(2);
           v(3)    0   -v(1);
          -v(2)  v(1)   0 ];
    R = eye(3) + vx + vx^2 * ((1-c)/(s^2)); % Rodrigues' rotation
end

%% --- Uniform random sampling of rays in cone ---
u = rand(num_samples,1);
v_rand = rand(num_samples,1);
cosTheta = (1 - u) + u*cos(theta_max);
theta = acos(cosTheta);
phi = 2*pi * v_rand;

dx = sin(theta).*cos(phi);
dy = sin(theta).*sin(phi);
dz = cos(theta);
dirs_local = [dx dy dz]';   % 3 x N
dirs_world = R * dirs_local;
negs  = ones(1,num_samples);
negs(1,2:2:end) = -1;
dirs_world = dirs_world.*negs;
% --- Draw rays ---
for i = 1:num_samples
    dir = dirs_world(:,i);
    P = r_sphere * dir;
    plot3([0 P(1)], [0 P(2)], [0 P(3)], 'k-', 'LineWidth', 0.8);
    plot3(P(1), P(2), P(3), 'ro', 'MarkerSize', 3);
end

%% --- Draw spherical cap (cone footprint on sphere) ---
n_cap = 50; % mesh resolution
theta_cap = linspace(0, theta_max, n_cap);
phi_cap = linspace(0, 2*pi, n_cap);
[TH, PH] = meshgrid(theta_cap, phi_cap);

Xc = r_sphere * sin(TH) .* cos(PH);
Yc = r_sphere * sin(TH) .* sin(PH);
Zc = r_sphere * cos(TH);

dirs_cap = [Xc(:)'; Yc(:)'; Zc(:)'];
dirs_cap_rot = R * dirs_cap;
Xc_rot = reshape(dirs_cap_rot(1,:), size(Xc));
Yc_rot = reshape(dirs_cap_rot(2,:), size(Yc));
Zc_rot = reshape(dirs_cap_rot(3,:), size(Zc));

% surf(Xc_rot, Yc_rot, Zc_rot, 'FaceColor',[1 0 0], ...
%     'FaceAlpha',0.3,'EdgeColor','none');

%% --- Draw transparent cone volume ---
r_base = r_sphere * sin(theta_max);
h_cone = r_sphere * cos(theta_max);
n_cone = 50;
theta_circle = linspace(0, 2*pi, n_cone);

Xb = r_base * cos(theta_circle);
Yb = r_base * sin(theta_circle);
Zb = h_cone * ones(1,n_cone);
Xa = zeros(1,n_cone);
Ya = zeros(1,n_cone);
Za = zeros(1,n_cone);

Xc_cone = [Xa; Xb];
Yc_cone = [Ya; Yb];
Zc_cone = [Za; Zb];

dirs_cone = [Xc_cone(:)'; Yc_cone(:)'; Zc_cone(:)'];
dirs_cone_rot = R * dirs_cone;
Xc_rot = reshape(dirs_cone_rot(1,:), size(Xc_cone));
Yc_rot = reshape(dirs_cone_rot(2,:), size(Yc_cone));
Zc_rot = reshape(dirs_cone_rot(3,:), size(Zc_cone));

surf(Xc_rot, Yc_rot, Zc_rot, 'FaceColor',[1 0.5 0], ...
    'FaceAlpha',0.3,'EdgeColor','none');

%% --- Compute metrics ---
% Sampling ratio
ratio = (1 - cos(theta_max)) / 2;
fprintf('Sampling ratio = %.4f (%.2f%% of full sphere)\n', ratio, 100*ratio);

% Spherical cap area
cap_area = 2 * pi * r_sphere^2 * (1 - cos(theta_max));
fprintf('Spherical cap area = %.4f\n', cap_area);

% Cone volume
V_cone = (1/3) * pi * r_base^2 * h_cone;
fprintf('Cone volume = %.4f\n', V_cone);

title(['Cone sampling: axis = [' num2str(axis_dir) ...
       '], \theta_{max} = ' num2str(theta_max_deg) '°']);
%% save
if bSAVE
dirs_world = dirs_world';%%
fid = fopen([sPath,sName], 'w');

fprintf(fid, ['[directions=',num2str(size(dirs_world,1)*2+1),']\n']);
fprintf(fid, 'CoordinateSystem = xyz\n');
fprintf(fid, 'Normalisation = none\n');
% fprintf(fid, 'Vector[0] = ( 0.00, 0.0, 0.0 ) \n');


for ii=1:size(dirs_world,1)
    fprintf(fid, ['Vector[',num2str((ii-1)*2),'] = ( ',num2str(dirs_world(ii,1)),', ',num2str(dirs_world(ii,2)),', ',num2str(dirs_world(ii,3)),' ) \n']);
    fprintf(fid, ['Vector[',num2str(ii*2-1),'] = ( ',num2str(dirs_world(ii,1)),', ',num2str(dirs_world(ii,2)),', ',num2str(dirs_world(ii,3)),' ) \n']);
end
fclose(fid);




fid = fopen([sPath,sNameT], 'w');
for ii=1:size(dirs_world,1)
    fprintf(fid, [num2str(dirs_world(ii,1)),' ',num2str(dirs_world(ii,2)),' ',num2str(dirs_world(ii,3)),' \n']);
end
fclose(fid);

end