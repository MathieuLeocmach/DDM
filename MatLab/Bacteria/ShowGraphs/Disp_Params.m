projectpath = genpath(pwd);
addpath(projectpath);

FileID = fopen('EmplacementData.txt');
EmplacementData = textscan(FileID,'%s');
EmplacementData = EmplacementData{1};
EmplacementData = [EmplacementData{1}];
fclose(FileID);

load([EmplacementData,'Variables.mat']);


%Choix des bornes pour le fit
TabqPixelLog = log(qs);
ParamsLog = log(Params(:,6));
loglog(qs, Params(:,6));

%Vérification des bornes
BorneInf = 0.114; %Pour passer aux vraies valeurs de q, prendre exp(.)*10^3
BorneSup = 1.636;


figure(1);
clf;
set(gcf,'paperpositionmode','auto','position',[700 100 600 900]);
Fig1 = subplot(5,1,2)
PositionFig1 = get(Fig1,'position') 
subplot('position',[0.15 PositionFig1(2) 0.8 0.333]) %dimension in the window

loglog(qs*10^3, Params(:,6), 'o');
ylimite = get(gca,'ylim');


IndiceBorneInf = find(abs(qs*1e3-BorneInf*ones(1,ImageSize/2-1)) == min(abs(qs*1e3-BorneInf*ones(1,ImageSize/2-1))))
IndiceBorneSup = find(abs(qs*1e3-BorneSup*ones(1,ImageSize/2-1)) == min(abs(qs*1e3-BorneSup*ones(1,ImageSize/2-1))))

TabqPixelLog = TabqPixelLog(IndiceBorneInf:IndiceBorneSup);
ParamsLog = ParamsLog(IndiceBorneInf:IndiceBorneSup);

x0Diff = 1;
FitVitesse = @(a,xdata) -xdata+a;
xVitesse = lsqcurvefit(FitVitesse,x0Diff,TabqPixelLog',ParamsLog);

v = exp(-xVitesse) % vitesse moyenne en nm/s;

hold on
Theor = loglog(1e3*qs(IndiceBorneInf:IndiceBorneSup), exp(FitVitesse(xVitesse, TabqPixelLog)), 'k','LineWidth',2);

xlimite = get(gca,'xlim');
ylimite = get(gca,'ylim');

PourcentagePosition = 0.08;

set(gca,'xticklabel',[]);
set(gca,'XTick',[0.2 0.4 0.8 1.5 3])
set(gca,'FontSize',18);
set(gca,'fontname','times');
save([EmplacementData,'VitesseMoyenne.mat'], 'v');

run('FitDiffusion2.m');
clear all
run('Alpha_et_Z_2.m');
clear all

print -depsc2 VitesseDiffusionParams.eps;