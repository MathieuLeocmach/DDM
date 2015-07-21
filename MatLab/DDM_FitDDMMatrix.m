%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential Dynamics code 2
% authors: David Germain, Mathieu Leocmach and Thomas Gibaud
% 
% Fit the DDM matrix with the appropriate model
%
% OUTPUT: Fit parameters for Colloids or Bacteria
% 
% You have to put your own parameters in the INPUT section, in the
% variable Params0 which initializes the fit and lb and ub which difine the
% upper and lower boundary of the fit parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
projectpath = genpath(pwd);
addpath(projectpath);

%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load DDMVariable.mat obtained from the first DDM codes located in the folder SaveFolder
FolderSave='C:\thomas\research\david\GraphColloides\results2\';
load([FolderSave,'DDMVariables.mat']);

%Choose the correct function to fit the ISF
FitChoice = 1; %1=Bacteria, 2=Colloids
dtLimit = 52; %fit the DDM matrix up to dt=dtLimit (max value = length(dtMerge))



%% Fit function for the DDMMerge matrix
MatrixFit=zeros(ImageSize/2-1,dtLimit); %initialization
ISF_Fit=MatrixFit;%initialization
Noise = mean(DDMMerge(:,ImageSize/2-1)); % Noise floor
%fit function for the Bacteria (ref: )
Pv = @(a,xdata) ( (a(6)+1)*a(5)./(a(6)*xdata) )...
    .* sin( a(6)*atan(xdata /(a(5)*(a(6)+1))))...
    ./( (1+ (xdata/(a(5)*(a(6)+1))).^2).^(a(6)/2));
ISF_Bact = @(a,xdata) exp(-xdata/a(3)).*((1-a(4))+a(4)*Pv(a, xdata));
FitBacteria = @(a,xdata) log(a(1)*(1-ISF_Bact(a, xdata))+a(2));
%fit function for the colloids
FitColloid = @(a,xdata) log(a(1)*(1-exp(-xdata/a(3)))+a(2));    


%% Fit Bacteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FitChoice == 1  
    Params=zeros(ImageSize/2-1,6); %initialization
   
    for Qinter =1:ImageSize/2-1
        %fit parameters initialization, lower and upper boundaries
            Params0 = [max(DDMMerge(1:dtLimit,Qinter))-min(DDMMerge(1:dtLimit,Qinter)),...
                        Noise, 1./(0.2*1e6*qs(Qinter).^2),...
                        0.5,...
                        100*qs(Qinter),...
                        1];
            lb = [(max(DDMMerge(1:dtLimit,Qinter))-Noise)*0.8, Noise*0.8, 0, 0, 0, 0];
            ub = [(max(DDMMerge(1:dtLimit,Qinter))-Noise)*1.2, Noise*1.2, 1000, 1, 10000*qs(Qinter)*10, 5];

            Params(Qinter, :) = lsqcurvefit(FitBacteria,Params0,dtMerge(1:dtLimit),log(DDMMerge(1:dtLimit,Qinter)'),lb,ub); %fit
            MatrixFit(Qinter,:) = exp(FitBacteria(Params(Qinter, :),dtMerge(1:dtLimit))); % DDM matrix fit
            ISF_Fit(Qinter,:) =1+(Params(Qinter, 2)-MatrixFit(Qinter,:))/Params(Qinter, 1); % ISF fit
    end
end;



%% Fit colloids    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FitChoice == 2 
    Params=zeros(ImageSize/2-1,3); % initialization
      
    for Qinter =1:ImageSize/2-1
        %fit parameters initialization, lower and upper boundaries
        Params0 = [max(DDMMerge(1:dtLimit,Qinter))*2,...
                    Noise,...
                    1];
        lb = [(max(DDMMerge(1:dtLimit,Qinter))-Noise)*0.8, Noise*0.8,0];
        ub = [(max(DDMMerge(1:dtLimit,Qinter))-Noise)*1.2, Noise*1.2,1000];   
        Params(Qinter, :) = lsqcurvefit(FitColloid,Params0,dtMerge(1:dtLimit),log(DDMMerge(1:dtLimit,Qinter)'),lb,ub); %fit
        MatrixFit(Qinter,:) = exp(FitColloid(Params(Qinter,:),dtMerge(1:dtLimit))); % DDM matrix fit
        ISF_Fit(Qinter,:)= exp(-(dtMerge(1:dtLimit)./Params(Qinter,3))); % ISF fit
    end
end;



%% save fit parameters in DDMFit as a function of the wavevector qs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dlmwrite([FolderSave,'DataMergeTabqFit.txt'],MatrixFit)
cd(FolderSave);
save([FolderSave,'DDMFitResults.mat'], 'Params','FitColloid','FitBacteria','ISF_Bact','Pv','MatrixFit','ISF_Fit','dtLimit');
