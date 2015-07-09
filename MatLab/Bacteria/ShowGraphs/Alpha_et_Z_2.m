projectpath = genpath(pwd);
addpath(projectpath);

FileID = fopen('EmplacementData.txt');
EmplacementData = textscan(FileID,'%s');
EmplacementData = EmplacementData{1};
EmplacementData = [EmplacementData{1}];
fclose(FileID);

load([EmplacementData,'Variables.mat'])
load([EmplacementData,'VitesseMoyenne.mat']);

%Pourcentage de bactéries mobiles
Taux = 100*mean(Params(1:100,5))

%%Fig 1
Fig1 = subplot(5,1,3)
PositionFig1 = get(Fig1,'position') 
subplot('position',[0.15 PositionFig1(2)+0.01 0.8 0.75/5]) %dimension in the window

BorneMoyennageMin = 0.1142
BorneMoyennageMax = 1.636

A = semilogx(qs*10^3, Params(:,5), 'ko');
ylim([-0.08 1.08])
% ylim([-500 500])
ylimite = get(gca,'ylim')
xlim([BorneMoyennageMin BorneMoyennageMax])


ylabel('$\alpha$','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18)

IndiceBorneInf = find(abs(qs*10^3-BorneMoyennageMin*ones(1,ImageSize/2-1)) == min(abs(qs*10^3-BorneMoyennageMin*ones(1,ImageSize/2-1))));
IndiceBorneSup = find(abs(qs*10^3-BorneMoyennageMax*ones(1,ImageSize/2-1)) == min(abs(qs*10^3-BorneMoyennageMax*ones(1,ImageSize/2-1))));

Alpha = mean(Params(IndiceBorneInf:IndiceBorneSup,5))
hold on;
% semilogx([qs(IndiceBorneInf)*10^3 qs(IndiceBorneSup)*10^3], [Alpha Alpha], 'r', 'linewidth', 2)

xlimite = get(gca,'xlim')


PourcentagePositionX = 0.08
PourcentagePositionY = 0.885
PositionXFig = 10^(log10(xlimite(1)) + PourcentagePositionX*(4/6)*(log10(xlimite(2)/xlimite(1))))
PositionYFig = ylimite(1) + PourcentagePositionY*(ylimite(2)-ylimite(1))

text(PositionXFig,PositionYFig,'(b)','fontname','times','fontsize',15, 'HorizontalAlignment', 'center')


set(gca,'xticklabel',[])
set(gca,'XTick',[0.2 0.4 0.8 1.5 3])
set(gca,'FontSize',18)
set(gca,'fontname','times')

%%Figure 2

Fig2 = subplot(5,1,4)
PositionFig2 = get(Fig2,'position') 
subplot('position',[0.15 PositionFig2(2)+0.02 0.8 0.75/5]) %dimension in the window

ParamsZ = Params(:,7);
semilogx(qs*10^3, ParamsZ, 'ko')
Zmoy = mean(Params(IndiceBorneInf:IndiceBorneSup,7))
hold on
% semilogx([qs(IndiceBorneInf)*10^3 qs(IndiceBorneSup)*10^3], [Zmoy Zmoy], 'r', 'linewidth', 2)
ylabel('$Z$','interpreter', 'latex','fontsize',18)


xlim([BorneMoyennageMin BorneMoyennageMax])
xlimite = get(gca,'xlim')
ylim([-2 22])
ylimite = get(gca,'ylim')

PourcentagePositionX = 0.08
PositionXFig = 10^(log10(xlimite(1)) + PourcentagePositionX*(4/6)*(log10(xlimite(2)/xlimite(1))))
PositionYFig = ylimite(1) + PourcentagePositionY*(ylimite(2)-ylimite(1))

text(PositionXFig,PositionYFig,'(c)','fontname','times','fontsize',15, 'HorizontalAlignment', 'center')


set(gca,'xticklabel',[])
set(gca,'XTick',[0.2 0.4 0.8 1.5 3])
set(gca,'FontSize',18)
set(gca,'fontname','times')

%%Figure 3

Fig3 = subplot(5,1,5)
PositionFig3 = get(Fig3,'position') 
subplot('position',[0.15 PositionFig3(2)+0.03 0.8 0.75/5]) %dimension in the window

sigma = (v*10^(-3))./(sqrt(Params(:,7)+1))
sigmamoy = mean(sigma(IndiceBorneInf:IndiceBorneSup))
semilogx(qs*10^3, sigma, 'ko')
hold on;
% semilogx([qs(IndiceBorneInf)*10^3 qs(IndiceBorneSup)*10^3], [sigmamoy sigmamoy], 'r', 'linewidth', 2)
ylabel('$\sigma$ {($\mu$m/s)}','interpreter', 'latex','fontsize',18)
xlabel('$q$ {($\mu$m$^{-1}$)}','interpreter', 'latex','fontsize',18)

xlim([BorneMoyennageMin BorneMoyennageMax])
xlimite = get(gca,'xlim')
ylim([-5 45])
ylimite = get(gca,'ylim')

PourcentagePositionX = 0.08
PositionXFig = 10^(log10(xlimite(1)) + PourcentagePositionX*(4/6)*(log10(xlimite(2)/xlimite(1))))
PositionYFig = ylimite(1) + PourcentagePositionY*(ylimite(2)-ylimite(1))

text(PositionXFig,PositionYFig,'(d)','fontname','times','fontsize',15, 'HorizontalAlignment', 'center')

% set(gca,'xticklabel',[])
set(gca,'XTick',[0.2 0.4 0.8 1.5 3])
set(gca,'FontSize',18)
set(gca,'fontname','times')
set(gca,'XMinorTick','on')
hold on

save([EmplacementData,'Alpha.mat'], 'Alpha')
save([EmplacementData,'Z.mat'], 'Zmoy')

run('DistributionVitesse.m')