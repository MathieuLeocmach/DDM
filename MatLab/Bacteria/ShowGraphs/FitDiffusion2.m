clear all

projectpath = genpath(pwd);
addpath(projectpath);

FileID = fopen('EmplacementData.txt');
EmplacementData = textscan(FileID,'%s');
EmplacementData = EmplacementData{1};
EmplacementData = [EmplacementData{1}];
fclose(FileID);

load([EmplacementData,'Variables.mat'])

hold on

%Choix des bornes pour le fit
TabqPixelLog = log(qs);
ParamsLog = log(Params(:,4));
% plot(TabqPixelLog, ParamsLog)

%Vérification des bornes
BorneInf = 0.114; %Pour passer aux vraies valeurs de q, prendre exp(.)*10^3
BorneSup = 1.636;

loglog(qs*10^3, Params(:,4), 'ro')
ylim([7e-3 2500])
ylimite = get(gca,'ylim');
% xlimite = [min(qs*10^3) max(qs*10^3)]
xlimite = [BorneInf BorneSup]



IndiceBorneInf = find(abs(qs*1e3-BorneInf*ones(1,ImageSize/2-1)) == min(abs(qs*1e3-BorneInf*ones(1,ImageSize/2-1))));
IndiceBorneSup = find(abs(qs*1e3-BorneSup*ones(1,ImageSize/2-1)) == min(abs(qs*1e3-BorneSup*ones(1,ImageSize/2-1))));

TabqPixelLog = TabqPixelLog(IndiceBorneInf:IndiceBorneSup);
ParamsLog = ParamsLog(IndiceBorneInf:IndiceBorneSup);

x0Diff = 1;
FitDiffusion = @(a,xdata) -2*xdata+a;
xDiffusion = lsqcurvefit(FitDiffusion,x0Diff,TabqPixelLog',ParamsLog);
CoeffDiffusion = exp(-xDiffusion)*10^(-6)

%fit libre
% x0Diff = [1, 1];
% FitDiffusion = @(a,xdata) a(1)*xdata+a(2);
% x = lsqcurvefit(FitDiffusion,x0Diff,TabqPixelLog',ParamsLog)

hold on
Theor = loglog(1e3*qs(IndiceBorneInf:IndiceBorneSup), exp(FitDiffusion(xDiffusion, TabqPixelLog)), 'k','LineWidth',2);

% xlabel('$q$ {($\mu$m$^{-1}$)}','interpreter', 'latex','fontsize',18);
ylabel('$\tau_d$ {(s)}, $\tau_r$ {(s)}','interpreter', 'latex','fontsize',18);
xlim(xlimite);


xlimite = get(gca,'xlim');

PourcentagePositionX = 0.08;
PourcentagePositionY = 0.95
PositionXFig = 10^(log10(xlimite(1)) + PourcentagePositionX*(4/6)*(log10(xlimite(2)/xlimite(1))));
PositionYFig = 10^(log10(ylimite(1)) + PourcentagePositionY*(log10(ylimite(2)/ylimite(1))));
text(PositionXFig,PositionYFig,'(a)','fontname','times','fontsize',15, 'HorizontalAlignment', 'center');

% set(gca,'xticklabel',[]);
set(gca,'YTick',10.^[-1:3]);
set(gca,'FontSize',18);
set(gca,'fontname','times');
set(gca,'YTick',10.^[-2:3])

save([EmplacementData,'CoeffDiffusion.mat'], 'CoeffDiffusion')