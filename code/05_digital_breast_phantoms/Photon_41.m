function [] = Photon_41 ()
% Repository group: 05_digital_breast_phantoms
% Version role: primary simulation or study entry point
%
% Research note: This file is retained as a documented historical version.
% Review configuration, external data dependencies, and output paths before use.
% Impact of Tumor/Signal on DPF/PPF Series Simulations
% this version: reading phantoms prepared by Diego

clc
close all
set(0,'DefaultFigureWindowStyle','docked')

% Cup B / resolution: 0.2 mm / cancer diameter = 3 or 8 mm
% Differernt percent of
% Air = 0.0000
% Fat (Adipose Tissue) =  0.0010
% Fibroglandular Tissue = 0.9990
Nz = 1020; Ny = 257; Nx = 323; dl = 0.2/10; r_lowres = 0.35;

for Pd = [25,50,75]
    for ix = [0,1,2]
        % Define the input file path
        the_path = ['/home/kaiser/Phantoms/Diego/Phantoms/CupB_Pd',num2str(Pd),'_',num2str(ix),'_y_-z_x.obj'];
        % Open the file for reading
        fid = fopen(the_path, 'r');
        % Read the binary data as single precision (float32)
        data = fread(fid,'float32');
        % Close the file
        fclose(fid);
        % normalization of codes
        data(        data<=0.0) = 0.0;
        data(0.0<data&data<0.9) = 0.1;
        data(0.9<=data        ) = 0.2;
        data = round(data*10);
        % showing stat
        disp(['PD = ',num2str(Pd),', v = ',num2str(ix)])
        disp(['unique codes = ',num2str(unique(data).')])
        N_air = sum(data==0);
        N_fat = sum(data==1);
        N_gln = sum(data==2);
        disp(['N_air + N_fiber + N_Fat = ',num2str(N_air+N_gln+N_fat)])
        disp(['% of Fibroglandular = ',num2str(N_gln/(N_gln+N_fat))])
        disp(['% of Fat = ',num2str(N_fat/(N_gln+N_fat))])
        disp('-----------------------------------------------------------')
        clearvars fid N_air N_fat N_gln the_path

        % Reshape the data into a 3D array
        TheImage = reshape(data, [Nz, Ny, Nx]);
        TheImage_low = round(imresize3(TheImage, 0.35));
        disp(['unique codes = ',num2str(unique(TheImage_low(:)).')])
        N_air = nnz(TheImage_low(:)==0);
        N_fat = nnz(TheImage_low(:)==1);
        N_gln = nnz(TheImage_low(:)==2);
        disp(['N_air + N_fiber + N_Fat = ',num2str(N_air+N_gln+N_fat)])
        disp(['% of Fibroglandular = ',num2str(N_gln/(N_gln+N_fat))])
        disp(['% of Fat = ',num2str(N_fat/(N_gln+N_fat))])
        disp(['load = ',num2str((size(TheImage_low,1)*size(TheImage_low,2)*size(TheImage_low,3))/(191*191*100))])
        disp('-----------------------------------------------------------')
        disp('-----------------------------------------------------------')
        clearvars fid N_air N_fat N_gln the_path
        clearvars data

        for ik = 1:1:min([size(TheImage_low,1),size(TheImage_low,2),size(TheImage_low,3)])
            s_x = round(ik/min([size(TheImage_low,1),size(TheImage_low,2),size(TheImage_low,3)])*size(TheImage_low,1));
            s_y = round(ik/min([size(TheImage_low,1),size(TheImage_low,2),size(TheImage_low,3)])*size(TheImage_low,2));
            s_z = round(ik/min([size(TheImage_low,1),size(TheImage_low,2),size(TheImage_low,3)])*size(TheImage_low,3));

            % extract a XZ-slab
            slab = squeeze(TheImage_low(:,s_y,:));
            subplot(1,3,1), pcolor((1:size(TheImage_low,3)).*dl/r_lowres,(1:size(TheImage_low,1)).*dl/r_lowres,slab), axis equal, axis tight; axis xy; colorbar; shading interp
            hold on, plot([s_z s_z].*dl/r_lowres,[1 size(TheImage_low,1)].*dl/r_lowres,'r-.'), hold off
            hold on, plot([1 size(TheImage_low,3)].*dl/r_lowres,[s_x s_x].*dl/r_lowres,'r-.'), hold off
            xlabel('X (cm)')
            ylabel('Z (cm)')
            title(['Pd = ',num2str(Pd),'%, X = ',num2str(ix)])
            set(gca,'FontSize',18)
            clearvars slab

            % extract a YZ-slab
            slab = squeeze(TheImage_low(:,:,s_z));
            subplot(1,3,2), pcolor((1:size(TheImage_low,2)).*dl/r_lowres,(1:size(TheImage_low,1)).*dl/r_lowres,slab), axis equal, axis tight; axis xy; colorbar; shading interp
            hold on, plot([s_y s_y].*dl/r_lowres,[1 size(TheImage_low,1)].*dl/r_lowres,'r-.'), hold off
            hold on, plot([1 size(TheImage_low,2)].*dl/r_lowres,[s_x s_x].*dl/r_lowres,'r-.'), hold off
            xlabel('Y (cm)')
            ylabel('Z (cm)')
            title(['Pd = ',num2str(Pd),'%, X = ',num2str(ix)])
            set(gca,'FontSize',18)
            clearvars slab

            % extract a XY-slab
            slab = squeeze(TheImage_low(s_x,:,:));
            subplot(1,3,3), pcolor((1:size(TheImage_low,3)).*dl/r_lowres,(1:size(TheImage_low,2)).*dl/r_lowres,slab), axis equal, axis tight; axis xy; colorbar; shading interp
            hold on, plot([s_z s_z].*dl/r_lowres,[1 size(TheImage_low,2)].*dl/r_lowres,'r-.'), hold off
            hold on, plot([1 size(TheImage_low,3)].*dl/r_lowres,[s_y s_y].*dl/r_lowres,'r-.'), hold off
            xlabel('X (cm)')
            ylabel('Y (cm)')
            title(['Pd = ',num2str(Pd),'%, X = ',num2str(ix)])
            set(gca,'FontSize',18)
            clearvars slab

            pause(0.01)
            colormap parula
            clearvars s_x s_y s_z
        end
        clearvars iy
    end
end
end
