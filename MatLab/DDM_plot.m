%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential Dynamics code 4
% authors: David Germain, Mathieu Leocmach and Thomas Gibaud
%
% Graphs of the DDM analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
projectpath = genpath(pwd);
addpath(projectpath);

%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% results are in the folder SaveFolder
FolderSave='C:\thomas\research\david\GraphColloides\results2\';

%Indicate figX=1 if you want to generate the plot

%model independant
fig1=1; %DDM matrix 
fig2=0; %DDM vs q and vs dt 
%Colloids
fig3=0; %Brownian diffusion, tau and A, B vs q 
fig4=0; %ISF Brownian diffusion
%Bacteria
fig5=1; %ISF bacteria
fig6=1; % bact, tau and A, B vs q 
fig7=1;

%% Figure 1: DDM matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig1==1
    load([FolderSave,'DDMVariables.mat']);

    figure(1);
    clf
set(gcf,'paperpositionmode','auto','position',[1300 100 600 600])

h=surf(qs*1e3,dtMerge,DDMMerge,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
view([0 90])
set(gca, 'Layer','top')
grid off
set(get(h,'Parent'),'ZScale','log');
set(get(h,'Parent'),'XScale','log');
set(get(h,'Parent'),'YScale','log');


colormap(jet)
c = colorbar('location','eastoutside','YScale','log');
set(c,'YTick',[1,10,100,1000,10000,100000])
locate = get(c,'title');
ylabel(c, '$\mathcal{D}$ (a.u.)','interpreter', 'latex',  'FontSize',15,'fontname','times')

xlim([min(qs) max(qs)]*1e3)
ylim([min(dtMerge) max(dtMerge)])
zlim=([min(min(DDMMerge)) max(max(DDMMerge)) ]);

xlabel('$q$ ($\mu$m$^{-1}$)','interpreter', 'latex','fontsize',18)
ylabel('$\Delta t$ {(s)}','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18,'fontname','times')
end

%% Figure 2: DDM vs q and vs dt 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig2==1
    load([FolderSave,'DDMVariables.mat']);

figure(2);
clf
set(gcf,'paperpositionmode','auto','position',[1300 100 600 800])

% Figure a
Fig1 = subplot(2,1,1);
PositionFig1 = get(Fig1,'position') ;
skip=5;
set(0,'DefaultAxesColorOrder',jet(length(1:skip:length(dtMerge))));
subplot('position',[0.15 PositionFig1(2) 0.8 0.75/2]) %dimension in the window

loglog(qs*1000,DDMMerge(1:skip:end,:),'LineWidth',2); hold on;
xlabel('$q$ {($\mu$m$^{-1}$)}','interpreter', 'latex','fontsize',18)
ylabel('$\mathcal{D}$ {(a.u.)}','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18)

set(gca,'YTick',10.^(0:round(1+log10(max(max(DDMMerge))))))
xlim(1000*[min(qs)*.8 max(qs)*1.2])
ylim([min(min(DDMMerge))*.5 max(max(DDMMerge))*2])

set(gca,'fontname','times')
h = colorbar('location','West', 'Yscale', 'log');
caxis([dtMerge(1) dtMerge(end)])
set(h,'Position',[0.86 0.83 .02 .1], 'YTick',10.^(-2:2:2))
locate = get(h,'title');
set(locate,'pos',[-2 1 0],'string','$\Delta t$ (s)','interpreter', 'latex',  'FontSize',15,'fontname','times','Rotation',90);
set(h,'FontSize',12,'fontname','times')

% Figure b
Fig2 = subplot(2,1,2);
PositionFig2 = get(Fig2,'position'); 
skip=20;
set(0,'DefaultAxesColorOrder',jet(length(1:skip:length(qs))));
subplot('position',[0.15 PositionFig2(2) 0.8 0.75/2]) %dimension in the window

loglog(dtMerge, DDMMerge(:,1:skip:end),'LineWidth',2);
xlabel('$\Delta t$ {(s)}','interpreter', 'latex','fontsize',18)
ylabel('$\mathcal{D}$ {(a.u.)}','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18)
set(gca,'YTick',10.^(0:round(1+log10(max(max(DDMMerge))))), 'XTick', 10.^(round(log10(min(dtMerge))):round(1+log10(max(dtMerge)))));

h = colorbar('location','West');
caxis([qs(1) qs(end)]*10^3)
set(h,'Position',[0.25 0.35 .02 .1], 'YTick',10.^((log10(1000*min(qs))):(log10(1000*max(qs)))))
locate = get(h,'title');
PositionTitre = get(locate, 'position');
set(locate,'pos',[-2.9 2.55 5],'string','$q$ ($\mu$m$^{-1}$)','interpreter', 'latex',  'FontSize',15,'fontname','times','Rotation',90);

xlim([min(dtMerge)*.8 max(dtMerge)*1.2])
ylim([min(min(DDMMerge))*.5 max(max(DDMMerge))*2])
set(gca,'fontname','times')
set(h,'FontSize',12,'fontname','times')

end

%% Figure 3: Brownian diffusion, tau and A, B vs q 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig3==1
    load([FolderSave,'DDMVariables.mat']);
    load([FolderSave,  'DDMFitTau.mat']);
    load([FolderSave,'DDMFitResults.mat']);

figure(3);
clf
set(gcf,'paperpositionmode','auto','position',[1300 100 600 800])
% Figure a
Fig1 = subplot(2,1,1);
PositionFig1 = get(Fig1,'position'); 
subplot('position',[0.15 PositionFig1(2) 0.8 0.75/2])

loglog(qs*10^3, Params(:,3), 'ko');
hold on;
plot([1 1]*qs(nMin)*1e3,[1e-3 1e5], '--','LineWidth',2)
plot([1 1]*qs(nMax)*1e3,[1e-3 1e5], '--','LineWidth',2)
loglog(qs*1e3,1/DiffusionCoeff./(qs*1e3).^2 , 'r','LineWidth',2);

ylabel('$\tau_d$ {(s)}','interpreter', 'latex','fontsize',18)
xlim(1000*[min(qs)*.8 max(qs)*1.2])
ylim([min(Params(:,3))*.5 max(Params(:,3))*2])

set(gca,'xticklabel',[])
set(gca,'FontSize',18)
set(gca,'fontname','times')

% Figure b
Fig2 = subplot(2,1,2);
PositionFig2 = get(Fig2,'position'); 
subplot('position',[0.15 PositionFig2(2)+0.09 0.8 0.75/2]) %dimension in the window

loglog(qs*10^3, Params(:,1), 'ks');hold on;
loglog(qs*10^3, Params(:,2), 'kd');
plot([1 1]*qs(nMin)*1e3,[1e-3 1e5], '--','LineWidth',2)
plot([1 1]*qs(nMax)*1e3,[1e-3 1e5], '--','LineWidth',2)

xlim(1000*[min(qs)*.8 max(qs)*1.2])
ylim([min(Params(:,1))*.5 max(Params(:,1))*2])

xlabel('$q$ {($\mu$m$^{-1}$)}','interpreter', 'latex','fontsize',18)
ylabel('$A, B$ {(a.u.)}','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18,'fontname','times')
end


%% Figure 4: isf brownian motion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig4==1
    load([FolderSave,'DDMVariables.mat']);
    load([FolderSave,'DDMFitResults.mat']);
    load([FolderSave,'DDMFitTau.mat']);
    
ISF = 1-((DDMMerge'-Params(:,2)*ones(1,length(dtMerge)))./(Params(:,1)*ones(1,length(dtMerge))))';

ISF_fit=0*ISF;
for i=1:length(Params)
ISF_fit(:,i) = exp(-dtMerge/Params(i,3));
end


figure(4);
clf
set(gcf,'paperpositionmode','auto','position',[1300 100 600 600])
skip=10;
couleur = jet(length(nMin:skip:nMax));

xlabel('$\Delta t$ {(s)}','interpreter', 'latex','fontsize',18)
ylabel('$f(q, \Delta t)$','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18)

semilogx(dtMerge,ISF(:,nMin:skip:nMax),'o'); hold on;
semilogx(dtMerge,ISF_fit(:,nMin:skip:nMax))

h = colorbar('location','West');
caxis([qs(nMin)  qs(nMax)]*10^3)
set(h,'Position',[0.85 0.81 .02 .1], 'YTick',10.^(round(log10(1000*min(qs))):round(log10(1000*max(qs)))))
locate = get(h,'title');
PositionTitre = get(locate, 'position');
set(locate,'pos',[-2.9 1.51 1],'string','$q$ ($\mu$m$^{-1}$)','interpreter', 'latex',  'FontSize',15,'fontname','times','Rotation',90);

xlim([min(dtMerge)*.5 max(dtMerge)*2])
ylim([-0.1 1.1])
set(gca,'XTick',10.^(-3:5)) 

set(gca,'fontname','times')
set(h,'FontSize',12,'fontname','times')
xlabel('$\Delta t $ {(s)}','interpreter', 'latex','fontsize',18)
ylabel('$f(q, \Delta t)$','interpreter', 'latex','fontsize',18)

end


%% Figure 5: isf bacteria motion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig5==1
    load([FolderSave,'DDMVariables.mat']);
    load([FolderSave,'DDMFitResults.mat']);
    load([FolderSave,'DDMFitTau.mat']);
    
ISF = 1-((DDMMerge'-Params(:,2)*ones(1,length(dtMerge)))./(Params(:,1)*ones(1,length(dtMerge))))';



figure(5);
clf
set(gcf,'paperpositionmode','auto','position',[1300 100 600 600])
skip=10;
couleur = jet(length(nMin:skip:nMax));

xlabel('$\Delta t$ {(s)}','interpreter', 'latex','fontsize',18)
ylabel('$f(q, \Delta t)$','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18)

semilogx(dtMerge,ISF(:,nMin:skip:nMax),'o'); hold on;
 semilogx(dtMerge(1:dtLimit),ISF_Fit(nMin:skip:nMax,:)')

h = colorbar('location','West');
caxis([qs(nMin)  qs(nMax)]*10^3)
set(h,'Position',[0.85 0.81 .02 .1], 'YTick',10.^(round(log10(1000*min(qs))):round(log10(1000*max(qs)))))
locate = get(h,'title');
PositionTitre = get(locate, 'position');
set(locate,'pos',[-2.9 1.51 1],'string','$q$ ($\mu$m$^{-1}$)','interpreter', 'latex',  'FontSize',15,'fontname','times','Rotation',90);

xlim([min(dtMerge)*.5 max(dtMerge)*2])
ylim([-0.1 1.1])
set(gca,'XTick',10.^(-3:5)) 

set(gca,'fontname','times')
set(h,'FontSize',12,'fontname','times')
xlabel('$\Delta t $ {(s)}','interpreter', 'latex','fontsize',18)
ylabel('$f(q, \Delta t)$','interpreter', 'latex','fontsize',18)

end



%% Figure 6: bact, tau and A, B vs q 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig6==1
    load([FolderSave,'DDMVariables.mat']);
    load([FolderSave,  'DDMFitTau.mat']);
    load([FolderSave,'DDMFitResults.mat']);

figure(6);
clf
set(gcf,'paperpositionmode','auto','position',[1300 100 600 800])
% Figure a
Fig1 = subplot(2,1,1);
PositionFig1 = get(Fig1,'position'); 
subplot('position',[0.15 PositionFig1(2) 0.8 0.75/2])

loglog(qs*10^3, Params(:,3), 'ko');hold on;

loglog(qs*10^3, Params(:,5), 'ks');

plot([1 1]*qs(nMin)*1e3,[1e-3 1e5], '--','LineWidth',2)
plot([1 1]*qs(nMax)*1e3,[1e-3 1e5], '--','LineWidth',2)
loglog(qs*1e3,1/DiffusionCoeff./(qs*1e3).^2 , 'r','LineWidth',2);
loglog(qs*1e3,1/Velocity./(qs*1e3) , 'r','LineWidth',2);

ylabel('$\tau_d$ {(s)}','interpreter', 'latex','fontsize',18)
xlim(1000*[min(qs)*.8 max(qs)*1.2])
ylim([min(Params(:,5))*.5 max(Params(:,3))*5])

set(gca,'xticklabel',[])
set(gca,'FontSize',18)
set(gca,'fontname','times')

% Figure b
Fig2 = subplot(2,1,2);
PositionFig2 = get(Fig2,'position'); 
subplot('position',[0.15 PositionFig2(2)+0.09 0.8 0.75/2]) %dimension in the window

loglog(qs*10^3, Params(:,1), 'ks');hold on;
loglog(qs*10^3, Params(:,2), 'kd');
plot([1 1]*qs(nMin)*1e3,[1e-3 1e5], '--','LineWidth',2)
plot([1 1]*qs(nMax)*1e3,[1e-3 1e5], '--','LineWidth',2)

xlim(1000*[min(qs)*.8 max(qs)*1.2])
ylim([min(Params(:,1))*.5 max(Params(:,1))*2])

xlabel('$q$ {($\mu$m$^{-1}$)}','interpreter', 'latex','fontsize',18)
ylabel('$A, B$ {(a.u.)}','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18,'fontname','times')
end


%% Figure 7: bact, alpha, Z 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fig7==1
    load([FolderSave,'DDMVariables.mat']);
    load([FolderSave,  'DDMFitTau.mat']);
    load([FolderSave,'DDMFitResults.mat']);

figure(7);
clf
set(gcf,'paperpositionmode','auto','position',[1300 100 600 800])
% Figure a
Fig1 = subplot(2,1,1);
PositionFig1 = get(Fig1,'position'); 
subplot('position',[0.15 PositionFig1(2) 0.8 0.75/2])

semilogx(qs*10^3, Params(:,4), 'ko');hold on;


plot([1 1]*qs(nMin)*1e3,[1e-3 1e5], '--','LineWidth',2)
plot([1 1]*qs(nMax)*1e3,[1e-3 1e5], '--','LineWidth',2)

ylabel('$\alpha$ ','interpreter', 'latex','fontsize',18)
xlim(1000*[min(qs)*.8 max(qs)*1.2])
ylim([0 1.1])

set(gca,'xticklabel',[])
set(gca,'FontSize',18)
set(gca,'fontname','times')

% Figure b
Fig2 = subplot(2,1,2);
PositionFig2 = get(Fig2,'position'); 
subplot('position',[0.15 PositionFig2(2)+0.09 0.8 0.75/2]) %dimension in the window

loglog(qs*10^3, Params(:,6), 'ks');hold on;
plot([1 1]*qs(nMin)*1e3,[1e-3 1e5], '--','LineWidth',2)
plot([1 1]*qs(nMax)*1e3,[1e-3 1e5], '--','LineWidth',2)

xlim(1000*[min(qs)*.8 max(qs)*1.2])
ylim([.1 max(Params(:,6))*2])

xlabel('$q$ {($\mu$m$^{-1}$)}','interpreter', 'latex','fontsize',18)
ylabel('$Z$ ','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18,'fontname','times')
end

