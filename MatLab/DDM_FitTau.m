%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential Dynamics code 3
% authors: David Germain, Mathieu Leocmach and Thomas Gibaud
%
% Fit the decorellation characteristic time(s) with the appropriate model
%
% OUTPUT: diffusion coefficient (colloids and bacteria), 
% mean velocity (bacteria)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
projectpath = genpath(pwd);
addpath(projectpath);

%% INPUT
% load results obtained from the DDM codes 1 and 2 located in the folder SaveFolder
FolderSave='C:\thomas\research\david\GraphColloides\results2\';
load([FolderSave,'DDMVariables.mat']);
load([FolderSave,'DDMFitResults.mat']);

%fit choice
FitChoice = 1; %1=Bacteria, 2=Colloids

%Define the q-boundary to fit tau
qMin = .15; %um^-1
qMax = 2.5; %um^-1

nMin= find(1000*qs < qMin, 1, 'last' );
nMax= find(1000*qs < qMax, 1, 'last' );
Velocity=0;
DiffusionCoeff=0;

%% Fit colloids
if FitChoice == 2
%fit in log scale
x0Diff = 1;
FitDiffusion = @(a,xdata) -2*xdata + a;
xDif = lsqcurvefit( FitDiffusion, x0Diff,log10(qs(nMin:nMax)),log10(Params(nMin:nMax,3))' );
DiffusionCoeff = 10^(-xDif)*1e-6; %um2/s
end

%% Fit bacteria
if FitChoice == 1
%fit in log scale diffusion
x0Diff = 1;
FitDiffusion = @(a,xdata) -2*xdata + a;
xDif = lsqcurvefit( FitDiffusion, x0Diff,log10(qs(nMin:nMax)),log10(Params(nMin:nMax,3))' );
DiffusionCoeff = 10^(-xDif)*1e-6; %um2/s
%fit in log scale velocity
x0Vel = 20;
Fitvelocity = @(a,xdata) -1*xdata + a;
xVel = lsqcurvefit( Fitvelocity, x0Vel,log10(qs(nMin:nMax)),log10(Params(nMin:nMax,5))' );
Velocity= 10^(-xVel)*1e-3; %um/s


end

save([FolderSave,'DDMFitTau.mat'], 'DiffusionCoeff', 'Velocity', 'nMin', 'nMax')
