projectpath = genpath(pwd);
addpath(projectpath);

FileID = fopen('EmplacementData.txt');
EmplacementData = textscan(FileID,'%s');
EmplacementData = EmplacementData{1};
EmplacementData = [EmplacementData{1}];
fclose(FileID);

load([EmplacementData,'Variables.mat'])
IntervalA = Params(:,2);

figure(1)
clf
set(gcf,'paperpositionmode','auto','position',[700 0 600 400])
set(0,'DefaultAxesColorOrder',jet(length([1:4:56])));
axes('position',[0.15 0.17 0.74 0.75]) %dimension in the window

loglog(qs*10^3, IntervalA,'Linestyle' , '.', 'color', 'black', 'Marker', '*')
hold on;
loglog(qs*10^3, Params(:,3), 'Linestyle' , '.', 'color', 'blue', 'Marker', '+')
xlabel('$q$ {($\mu$m$^{-1}$)}','interpreter', 'latex','fontsize',18)
ylabel('$A$, $B$ {(a.u.)}','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18,'fontname','times')
set(gca,'YTick',10.^[0:5])

xlim([min(qs*10^3) max(qs*10^3)])
ylim([0.8 5e2])
xlimite = get(gca,'xlim')
ylimite = get(gca,'ylim')

BorneMoyennageMin = 0.1142
BorneMoyennageMax = 1.636
hold on;
semilogx([1 1]*BorneMoyennageMin, ([ylimite(1) ylimite(2)]),'k:')
hold on;
semilogx([1 1]*BorneMoyennageMax, ([ylimite(1) ylimite(2)]),'k:')

PositionX1 = 10^((log10(xlimite(1))+log10(BorneMoyennageMin))/2)
PositionX2 = 10^((log10(BorneMoyennageMin)+log10(BorneMoyennageMax))/2)
PositionX3 = 10^((log10(xlimite(2))+log10(BorneMoyennageMax))/2)
PositionY = 10^(0.95*log10(ylimite(2))+0.05*log10(ylimite(1)))

text(PositionX1,PositionY,'(1)','fontname','times','fontsize',15, 'HorizontalAlignment', 'center')
text(PositionX2,PositionY,'(2)','fontname','times','fontsize',15, 'HorizontalAlignment', 'center')
text(PositionX3,PositionY,'(3)','fontname','times','fontsize',15, 'HorizontalAlignment', 'center')


%Toute la suite du code sert à la modélisation
TabqPixelLog = log(qs);
% plot(TabqPixelLog, IntervalA)
% semilogx(qs*10^3, IntervalA)

BorneInf = -8.15;
BorneSup = -6.42;

IndiceBorneInf = find(abs(TabqPixelLog-BorneInf*ones(1,ImageSize/2-1)) == min(abs(TabqPixelLog-BorneInf*ones(1,ImageSize/2-1))));
IndiceBorneSup = find(abs(TabqPixelLog-BorneSup*ones(1,ImageSize/2-1)) == min(abs(TabqPixelLog-BorneSup*ones(1,ImageSize/2-1))));

TabqPixelPlot = exp(TabqPixelLog(IndiceBorneInf:IndiceBorneSup));
IntervalA = IntervalA(IndiceBorneInf:IndiceBorneSup);


x0A = [1,0.034,0.425,0.01,51000] %Données du papier/Marche pour bactéries
% x0A = [1,0.5,0.425,0.01,51000] %Marche pour Colloides
ParamsA = lsqcurvefit(@FitA,x0A,TabqPixelPlot,IntervalA')
hold on;
CourbeA = loglog(qs*10^3, FitA(ParamsA, qs), '--r', 'linewidth', 2);

PourcentagePosition = 0.08
PositionXFig = 10^(log10(xlimite(1)) + PourcentagePosition*(4/6)*(log10(xlimite(2)/xlimite(1))))
PositionYFig = 10^(log10(ylimite(1)) + PourcentagePosition*(log10(ylimite(2)/ylimite(1))))
% text(PositionXFig,PositionYFig,'(d)','fontname','times','fontsize',15, 'HorizontalAlignment', 'center')
clear all

print -depsc2 A+B+Fit_bact