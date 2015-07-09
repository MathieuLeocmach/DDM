projectpath = genpath(pwd);
addpath(projectpath);

FileID = fopen('EmplacementData.txt');
EmplacementData = textscan(FileID,'%s');
EmplacementData = EmplacementData{1};
EmplacementData = [EmplacementData{1}];
fclose(FileID);

load([EmplacementData,'Variables.mat'])



IndiceDebut = 32;
IndiceFin = 128;
IndicePas = 16;


Tabf = 1-((DDMCropped'-Params(:,3)*ones(1,length(dtCropped)))./(Params(:,2)*ones(1,length(dtCropped))));
TabfFit = 1-((MatrixFit-Params(:,3)*ones(1,length(dtCropped)))./(Params(:,2)*ones(1,length(dtCropped))));

figure(1);
clf
set(gcf,'paperpositionmode','auto','position',[700 100 600 400])
couleur = jet(length([IndiceDebut:IndicePas:IndiceFin]))
c=0 %Iteration sur la couleur
set(0,'DefaultAxesColorOrder',couleur);

axes('position',[0.15 0.17 0.8 0.75]) %dimension in the window
for i=[IndiceDebut:IndicePas:IndiceFin]
    c=c+1
%     semilogx(dtCropped, Tabf(i,:), 'o','Color',couleur(c,:))
%     semilogx(dtCropped*((10^3)*qs(i)), Tabf(i,:), 'o','Color',couleur(c,:))
    semilogx(dtCropped.*((10^3)*qs(i)).^2, Tabf(i,:), 'o','Color',couleur(c,:))
    hold on
%     semilogx(dtCropped, TabfFit(i,:),'Color',couleur(c,:))
%     semilogx(dtCropped*((10^3)*qs(i)), TabfFit(i,:),'Color',couleur(c,:))
    semilogx(dtCropped*((10^3)*qs(i)).^2, TabfFit(i,:),'Color',couleur(c,:))
end

% xlabel('$\Delta t$ {(s)}','interpreter', 'latex','fontsize',18)
% xlabel('$\Delta t q$ {(s/$\mu$m}{)}','interpreter', 'latex','fontsize',18)
xlabel('$\Delta t q^2$ {(s/$\mu$m}$^{2}${)}','interpreter', 'latex','fontsize',18)
ylabel('$f(q, \Delta t)$','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18)
set(gca,'XTick',10.^[-4:4]) %

h = colorbar('location','West');
caxis([qs(min([IndiceDebut:IndicePas:IndiceFin]))*10^3 qs(max([IndiceDebut:IndicePas:IndiceFin]))*10^3])
set(h,'Position',[0.85 0.67 .02 .2])

locate = get(h,'title')

set(locate,'pos',[-2.9 1.51 1],'string','$q$ ($\mu$m$^{-1}$)','interpreter', 'latex',  'FontSize',15,'fontname','times','Rotation',90);

% ylim([-0.1 1.1]) %Bonne valeur poure tracer en fonction de t
% xlim([1e-3 7.5e2])

% ylim([-0.1 1.1]) %Bonne valeur poure tracer en fonction de tq
% xlim([0.5e-3 2e3])

ylim([-0.1 1.1]) %Bonne valeur poure tracer en fonction de tq²
xlim([3.5e-4 5e3])

set(gca,'XTick',10.^[-3:3]) %Bon pour Delta t

xlimite = get(gca,'xlim')
ylimite = get(gca,'ylim')

PourcentagePosition = 0.08
PositionXFig = 10^(log10(xlimite(1)) + PourcentagePosition*(4/6)*(log10(xlimite(2)/xlimite(1))))
PositionYFig = (ylimite(1) + PourcentagePosition*(ylimite(2)-ylimite(1)))
text(PositionXFig,PositionYFig,'(c)','fontname','times','fontsize',15, 'HorizontalAlignment', 'center')


set(gca,'fontname','times')
set(h,'FontSize',12,'fontname','times')
clear all

print -depsc2 ISFBacteries3.eps