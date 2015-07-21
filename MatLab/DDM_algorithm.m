%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential Dynamics code 1
% authors: David Germain, Mathieu Leocmach and Thomas Gibaud
% 
% Image analysis of the DDM experiment
%
% OUTPUT
% D(dt,q)=abs(TFF(I(t+dt,q)-I(t,q)))^2 is generated
% from 2 stacks of images at Frequency1 and Frequency2
% D is radial averaged on q.
% dt is log spaced, 15 points/decade.
% There is a time average over MaxNCouples images at max.
% Data is saved in FolderSave
% 
% You have to put your own parameters in the INPUT section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
projectpath = genpath(pwd);
addpath(projectpath);

%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Pass to Folders containing stacks at Frequency1 and Frequency2.
% - Frequency1>Frequency2.
% - Timewise, the end of the first stack must overlap with the begining of the second
% stack so that the stack can be merged.
% - Both stacks must contain the same number of images, NbImage, with the same
% dimensions, ImageSize, and the same pixel size, PixelSize.
% - The stacks are decomposed in NbImage images in tif format in the folders
% ImageName= the name that preceed the numbers 00001.tif
% Folder1 and Folder2:

Folder1='C:\thomas\research\david\GraphColloides\Bacteria\1103_50µL_1h30_512x512_4000Im_400_6h30att\'; % Folders containing stacks at Frequency1
Folder2='C:\thomas\research\david\GraphColloides\Bacteria\1103_50µL_1h30_512x512_4000Im_4_6h30att\'; % Folders containing stacks at Frequency2
ImageName1='11031_'; %image base name in folder1
ImageName2='11032_'; %image base name in folder1
NbImage = 4000; %total number of images you want to analyse in the stacks
Frequency1 = 400; %(Hz)
Frequency2 = 4; %(Hz)
PixelSize = 6450/10; %Conversion 1 pixel -> nm
MaxNCouples=300; %Analysis parameters: average is performed over MaxNCouples images at most for each dt, MaxNCouples=<NbImage-1

FolderSave='C:\thomas\research\david\GraphColloides\results2\'; %Folder where you want to save the OUTPUT
mkdir(FolderSave); %create folder




%% Acquisition parameters obtained from INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd Folder1;  
Folder = {Folder1,Folder2};
ImageName = {ImageName1, ImageName2};
Im0Info = imfinfo([Folder{1},ImageName{1},sprintf('%05d',1),'.tif']); % first image to get info
Frequency = [Frequency1,Frequency2]; %acquisition frequencies of the 2 stack (Hz)
ImageSize = Im0Info.Width; %pixel
qs = ((2*pi)/(ImageSize*PixelSize))*(1:ImageSize/2); %generate wavevector q in physical units (nm-^1)

% Generate log scale for the times idt, 15 points per decade
Base = 10^(1/12); %power scale to get dt
Ndt = floor(log(NbImage)/log(Base)); %max power for NbImage
idts=unique(floor(Base.^(1:Ndt)));

% Initialize DDM matrix and dt
DDM = zeros(2,length(idts), ImageSize/2);
dt = zeros(2,length(idts));

%%DDM algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop on stack 1 and 2
for i=1:2; 
    dt(i,:) = idts/Frequency(i); %dt in real time
    
    % Loop on idts
    for p = 1:length(idts); 
        idt = idts(p); %real time interval between 2 images
        disp([i idt]);
        AvgFFT  = zeros(ImageSize,ImageSize);
        
        % Spread initial times over the available range
        Increment = max(floor((NbImage - idt)/MaxNCouples), 1);
        InitialTimes = 1:Increment:NbImage-idt;
        
        % loop on t
        for t=InitialTimes; 
            % Load 2 images to make the difference
            Im2 = im2double(imread([Folder{i},ImageName{i},sprintf('%05d',t+idt),'.tif'])); 
            Im1 = im2double(imread([Folder{i},ImageName{i},sprintf('%05d',t),'.tif'])); 
            AvgFFT = AvgFFT + abs(fftshift(fft2(Im2-Im1))).^2;% FFT of differences squared, and add up
        end;
       
        AvgFFT = AvgFFT./length(InitialTimes); % Divide by the number of image differences to get the average
        
        % Replace values at q = 0 by the mean of its neighbour (A modifier)
        AvgFFT(:,(ImageSize/2)+1) = 1/2*(AvgFFT(:,ImageSize/2)+AvgFFT(:,(ImageSize/2)+2)).*ones(ImageSize,1);
                
        [DDM(i, p, :), VecQ] =  radialavg(AvgFFT,ImageSize/2); % Radial average of the FFT (A changer VecQ)
    end;
end

%% Merging of the 2 sets of acquisitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the closest time at Frequency1 to the smallest time at Frequency2
[Trash, Boundary] = min(abs(dt(1,:)-dt(2,1)));

% Rescale the value of radial average at 4 Hz according to the value at Boundary for 400Hz
for i=1:ImageSize/2
    DDM(2,:,i) = DDM(2,:,i)*(DDM(1, Boundary,i)/DDM(2,1,i));
end

% Generate a smooth transition between the 2 curves on the first third of their overlap
Overlap = floor(length(DDM(1,Boundary:end,1))/3); 

Spline = zeros(2, ImageSize/2); %initiatize matrix

for j=2:ImageSize/2  
    xx400 = min(dt(1,Boundary:Boundary+Overlap)):.1:max(dt(1,Boundary:Boundary+Overlap));
    yy400 = spline(dt(1,Boundary:Boundary+Overlap),DDM(1, Boundary:Boundary+Overlap,j),xx400);
    xx4 = min(dt(1,Boundary:Boundary+Overlap)):.1:max(dt(1,Boundary:Boundary+Overlap));
    yy4 = spline(dt(2,1:Overlap),DDM(2, 1:Overlap,j),xx4);

    for i =1:length(xx4)
        Weight = (length(xx4)-i)/(length(xx4)-1);
        Spline(i, j)= Weight*yy400(i)+(1-Weight)*yy4(i);
    end
end



% Merge 400Hz, transition and 4Hz
DDMMerge = [squeeze(DDM(1,dt(1,:)<xx4(1),:)) ; Spline; squeeze(DDM(2,dt(2,:)>xx4(end), :))];
dtMerge = [squeeze(dt(1,dt(1,:)<xx4(1))),xx4,squeeze(dt(2,dt(2,:)>xx4(end)))];

% Spline to decimate to 10 points per decades
dtMergeLog = log10(dtMerge);
dtSpline = min(dtMergeLog):.1:max(dtMergeLog);

DDMSpline= zeros(ImageSize/2, length(dtSpline)); %initiatize matrix
DDMMerge(:,1) = zeros(length(DDMMerge(:,1)),1); %La premiere colonne est remplie de NaN et de 0, les NaN posent problèmes alors je force que des 0
for i = 1:length(DDMMerge(2,:))
    DDMSpline(i,:) = spline(dtMergeLog,DDMMerge(:,i),dtSpline);
end
dtMerge = 10.^(dtSpline);
DDMMerge = DDMSpline';

% Remove the first q because it is irelevant
DDMMerge(:,1) = [];
qs(1)=[];


%% Save in the file DDMVariable.mat in the folder FolderSave: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of images per stacks, image size, the 2 acquisition
% frequencies, the pixel size, the wavector qs, the merge time dtMerge, 
% the merged DDM matrix and FolderSave.

cd(FolderSave);
save([FolderSave,'DDMVariables.mat'], 'NbImage', 'ImageSize', 'Frequency', 'PixelSize', 'qs', 'dtMerge', 'DDMMerge','FolderSave');



