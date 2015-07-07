projectpath = genpath(pwd);
addpath(projectpath);



%% Acquisition parameters 
NbImage = 4000
ImageSize = 512
PixelSize = 6450/10 %Conversion 1 pixel = x nano_metres

%% Analysis parameters
MaxNCouples=300;

%% Generate log scale for Delta t, 10 points per decade
Base = 10^(1/12) %Définit l'échelle de puissance utilisée pour obtenir tau, ici on a exactement 15 points par décades
Ndt = floor(log(NbImage)/log(Base)) %Puissance max pour l'intervalle de temps
idts=unique(floor(Base.^([1:Ndt])));
%% Conversion from indices to physical values of q
qs = ((2*pi)/(ImageSize*PixelSize))*[1:ImageSize/2];

%% Folders containing stacks at 400 Hz and 4 Hz
Folder = {'D:\David\Acquisition\25_06\60mgL\400Hz\Chey_0_',...
          'D:\David\Acquisition\25_06\60mgL\4Hz\Chey_0_'};
Frequency = [400,4];

%% Initialize
DDM = zeros(2,length(idts), ImageSize/2);
dt = zeros(2,length(idts));
for i=1:2;
    %%Real time
    dt(i,:) = idts/Frequency(i);
    
    %% Loop on idts
    for p = 1:length(idts); %On itère sur les puissances
        idt = idts(p) %Vrai intervalle entre 2 images
        AvgFFT  = zeros(ImageSize,ImageSize);
        
        %% Spread initial times over the available range
        Increment = max(floor((NbImage - idt)/MaxNCouples), 1);
        InitialTimes = [1:Increment:NbImage-idt];
        
        %% loop on t
        for t=InitialTimes; %On itère sur les images en ne prenant, au maximum, que les MaxNCouples premières
            %% Load images
            Im2 = im2double(imread([Folder{i},sprintf('%05d',t+idt),'.tif'])); %On importe les 2 images que l'on va soustraire
            Im1 = im2double(imread([Folder{i},sprintf('%05d',t),'.tif']));
            %% FFT of differences squared, and add up
            AvgFFT = AvgFFT + abs(fftshift(fft2(Im2-Im1))).^2;
        end;
        %% Divide by the number of couples to get the average
        AvgFFT = AvgFFT./length(InitialTimes);
        
        %% Replace values at q = 0 by the mean of its neighbour (A modifier)
        AvgFFT(:,(ImageSize/2)+1) = 1/2*(AvgFFT(:,ImageSize/2)+AvgFFT(:,(ImageSize/2)+2)).*ones(ImageSize,1);
        
        %% Radial average of the FFT (A changer VecQ)
        [DDM(i, p, :), VecQ] =  radialavg(AvgFFT,ImageSize/2); %On fait la moyenne radiale, VeqQ contient la valeur des q pour chaque bin
    end;
end

%% Merging of the 2 sets of acquisitions

%% Find the closest time at 400Hz to the smallest time at 4Hz
[Trash, Boundary] = min(abs(dt(1,:)-dt(2,1)));

%% Rescale the value of radial average at 4 Hz according to the value at Boundary for 400Hz
for i=1:ImageSize/2
    DDM(2,:,i) = DDM(2,:,i)*(DDM(1, Boundary,i)/DDM(2,1,i));
end


%% generate a smooth transition between the 2 curves on the first third of their overlap
Overlap = floor(length(DDM(1,Boundary:end,1))/3);

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

%A la fin de ce paragraphe de code, j'obtient la zone de raccordement des 2
%courbes


%% Merge 400Hz, transition and 4Hz
DDMMerge = [squeeze(DDM(1,dt(1,:)<xx4(1),:)) ; Spline; squeeze(DDM(2,dt(2,:)>xx4(end), :))];
dtMerge = [squeeze(dt(1,dt(1,:)<xx4(1))),xx4,squeeze(dt(2,dt(2,:)>xx4(end)))];

%% Spline to decimate to 10 points per decades
dtMergeLog = log10(dtMerge);
dtSpline = min(dtMergeLog):.1:max(dtMergeLog);
DDMMerge(:,1) = zeros(length(DDMMerge(:,1)),1); %La premiere colonne est remplie de NaN et de 0, les NaN posent problèmes alors je force que des 0
for i = 1:length(DDMMerge(2,:))
    DDMSpline(i,:) = spline(dtMergeLog,DDMMerge(:,i),dtSpline);
end
dtMerge = 10.^(dtSpline);
DDMMerge = DDMSpline';

%% Remove the first q because it is irelevant
DDMMerge(:,1) = [];%Je supprime la première colonne, remplie de résultats abérrants
qs(1)=[];


%% Save the merged data
FileID = fopen('EmplacementData.txt');
EmplacementData = textscan(FileID,'%s');
EmplacementData = EmplacementData{1};
EmplacementData = [EmplacementData{1}];
fclose(FileID);


save([EmplacementData,'Variables.mat'], 'NbImage', 'ImageSize', 'Frequency', 'PixelSize', 'qs', 'dtMerge', 'DDMMerge')