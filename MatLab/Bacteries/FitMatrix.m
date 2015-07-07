%% Import data
projectpath = genpath(pwd);
addpath(projectpath);

FileID = fopen('EmplacementData.txt');
EmplacementData = textscan(FileID,'%s');
EmplacementData = EmplacementData{1};
EmplacementData = [EmplacementData{1}];
fclose(FileID);

load([EmplacementData,'Variables.mat'])

%% Maximum time index for fitting
TimeMax = 50 %50 pour les bactéries de l'article, 56 max

DDMCropped = DDMMerge(1:TimeMax,:);
dtCropped = dtMerge(1:TimeMax);

MatrixFit=zeros(size(DDMCropped'));

%%Fit function for D
Pv = @(a,xdata) ( (a(6)+1)*a(5)./(a(6)*xdata) )...
    .* sin( a(6)*atan(xdata /(a(5)*(a(6)+1))))...
    ./( (1+ (xdata/(a(5)*(a(6)+1))).^2).^(a(6)/2));
ISF = @(a,xdata) exp(-xdata/a(3)).*((1-a(4))+a(4)*Pv(a, xdata));
LogDDM = @(a,xdata) log(a(1)*(1-ISF(a, xdata))+a(2));

Params = zeros(ImageSize/2-1, 7);

for Qinter =1:ImageSize/2-1

    %% Initial guess for the fit
    x0 = [max(DDMCropped(:,Qinter))-min(DDMCropped(:,Qinter)), min(DDMCropped(:,Qinter)), 1/((qs(Qinter).^2)*17e4), 0.5, 1/(10000*qs(Qinter)), 1];
    %% lower and upper bounds for the fit
    lb=[(max(DDMCropped(:,Qinter))-min(DDMCropped(:,Qinter)))*0.8, min(DDMCropped(:,Qinter))*0.8, 0, 0, 0, -500];
    ub=[(max(DDMCropped(:,Qinter))-min(DDMCropped(:,Qinter)))*1.2, min(DDMCropped(:,Qinter))*1.2, 1000, 1, 1000, 500];
    
    %% Fit DDMCropped by D
    Params(Qinter, 2:7) = lsqcurvefit(FitBacterie,x0,dtCropped,log(DDMCropped(:,Qinter)'),lb,ub);    
    Params(Qinter, 1) = qs(Qinter);

    MatrixFit(Qinter,:) = exp(FitBacterie(Params(Qinter,2:end),dtCropped));
end;

save([EmplacementData,'Variables.mat'], 'NbImage', 'ImageSize', 'Frequency', 'PixelSize', 'qs', 'dtMerge', 'dtCropped', 'DDMMerge', 'DDMCropped', 'MatrixFit','Params')
clear all