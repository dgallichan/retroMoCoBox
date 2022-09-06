% simulate effect of motion on a 2D-image 

clear 
close all
run('addRetroMoCoBoxToPath.m')

rng(0); % resets random number generator to get the same 'random' numbers each time code is run

% load example
image_original = imread(fullfile('images','col27_slice80.png'));
image_original = double(image_original(end:-1:1,:,1)'); % png specifies R,G,B even for grayscale image

[Nx, Ny] = size(image_original);
% force the image to have odd dimensions as it makes life much easier...
if mod(Nx,2) == 1
    Nx = Nx+1;
    image_original(Nx,1) = 0;
end
if mod(Ny,2) == 1
    Ny = Ny+1;
    image_original(1,Ny) = 0;
end

% Add some random slowly varying phase to make it a bit more realistic
imPhase = imresize((2*rand(5)')-1,[Nx,Ny]);

image_original = image_original.*exp(1i*imPhase);

kdata_original = fft2s(image_original); 

figure(1)
set(gcf,'Position',[50 50 900 450])
clf
subplot1(1,3)
subplot1(1)
imab(abs(image_original))
colorbar
title('Original image (magnitude)')
subplot1(2)
imab(angle(image_original))
colorbar
title('Original image (phase)')
subplot1(3)
imab(logabs(kdata_original))
colormap(gray)
title('k-space (log magnitude)')
%% Simulate some motion


% Generate some 3 dof motion parameters, translations in x and y plus rotations

mpars = zeros(3,Ny); % translations (vox) and rotations (degrees)
% mpars(2,:) = 5*linspace(-1,1,Ny);
mpars(3,:) = 5*sin(pi*4*linspace(-1,1,Ny));
% mpars(1,round(.2*Ny):round(.31*Ny)) =10;
% mpars(1,:) = 2*rand(1,Ny);
% mpars(3,round(.3*Ny):end) =10;
% mpars(1,:) = 3*perlinNoise1D(Ny,2);
% mpars(2,:) = 3*perlinNoise1D(Ny,2);
% mpars(3,:) = 3*perlinNoise1D(Ny,8);
% mpars(3,:) = 15*perlinNoise1D(Ny,3.^[0:8]);
% mpars(1,:) = ones(1,Ny);
% mpars(3,:) = sin(2*pi*linspace(-1,1,Ny));
% mpars(3,:) = 3*sin(pi*5*linspace(-1,1,Ny));

% mpars(3,round(.2*Ny):round(.4*Ny)) =3;
% mpars(3,round(.4*Ny):round(.6*Ny)) =-2;
% mpars(3,round(.6*Ny):round(.8*Ny)) =1;

% apply threshold to motion, if desired, to make motion more abrupt:
mpar_threshold = 3;
mpars(mpars>mpar_threshold) = 15;
mpars(mpars<mpar_threshold) = 0;


% subtract the mean from each
mpars = mpars - mean(mpars,2);

% Perform the simulated motion corruption
[kdata_simMotion, st] = apply_2DsimulatedMotion(kdata_original,mpars);
% 'st' is the object used by the NUFFT and contains the full k-space
% sampling pattern
image_simMotion = ifft2s(kdata_simMotion);


figure(2)
set(gcf,'Position',[  0 0 950 950])
clf
s1 = subplot(231);
imab(image_original)
title('Original image')
s2 = subplot(232);
imab(image_simMotion)
title('Simulated motion-corruption')
s3 = subplot(233);
imab((real(image_simMotion)-real(image_original))*100/max(image_original(:)))
cbar = colorbar;
s1.Position = [.111 .5 .26 .39];
s2.Position = [.385 .5 .26 .39];
s3.Position = [.662 .5 .26 .39];
title('Difference (%)')
colormap(gray)
subplot(212)
plot(mpars.','linewidth',2)
grid on; grid minor
axis tight
title("Simulated motion parameters")
legend("x displacement","y displacement","rotation")
xlabel("k-space lines in PE direction")
ylabel('Displacement (voxels) or Rotation (degrees)')

% View effective k-space sampling pattern after motion
figure(3)
set(gcf,'Position',[0 0 950 950])
clf
scatter(st.om(:,1),st.om(:,2),15,'.')
axis equal tight
grid on; grid minor
title('Effective k-space sampling due to rotations')

%% Now test motion correction of the simulated corrupted data


image_simCorrected = apply_2DMotionCorrection(kdata_simMotion,mpars,1,2);
image_simCorrected_10iters = apply_2DMotionCorrection(kdata_simMotion,mpars,10,2);

figure(4)
set(gcf,'Position',[50 50 1600 520])
subplot1(1,4)
subplot1(1)
imab(image_original)
title('Original Image')
subplot1(2)
imab(image_simMotion)
title('Simulated motion-corruption')
subplot1(3)
imab(image_simCorrected)
title('Simulated correction after motion-corruption')
subplot1(4)
imab(image_simCorrected_10iters)
title({'Simulated correction after motion-corruption', '10 CG iterations'})

colormap(gray)


