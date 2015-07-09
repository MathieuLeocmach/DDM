projectpath = genpath(pwd);
addpath(projectpath);

FileID = fopen('EmplacementData.txt');
EmplacementData = textscan(FileID,'%s');
EmplacementData = EmplacementData{1};
EmplacementData = [EmplacementData{1}];
fclose(FileID);

load([EmplacementData,'Variables.mat'])
load([EmplacementData,'VitesseMoyenne.mat']);
load([EmplacementData,'Z.mat']);
Z=Zmoy
% Z = 4


axes('position',[0.30 0.375 0.15 0.05]) %dimension 
DistVitesse = [0:80000];

DistribShulz = Shulz(DistVitesse, Z, v) ./sum(Shulz(DistVitesse, Z, v));
Sigma = sqrt(sum(DistVitesse .* DistVitesse .* DistribShulz) - v^2)
fwhmax = fwhm(DistVitesse, DistribShulz) %Largeur totale à mi hauteur en unité de DistVitesse

Proba = plot(DistVitesse*10^-3, DistribShulz);
hold on;
VitesseMoyenne = plot([v,v]*10^-3, [0,max(DistribShulz)], '--r');

xlabel('$v$ {($\mu$m/s)}','interpreter', 'latex','fontsize',10)
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [0 0.000001 0])

ylabel('P','interpreter', 'latex','fontsize',10) %Pour tau_run
set(gca,'FontSize',10)
% xlim([4e-2 4.5])
ylim([0 4.5e-5])
set(gca,'fontname','times')