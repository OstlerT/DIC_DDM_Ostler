function [] = DDMPlots(DDM,qs,dt,nMin,nMax,DiffusionCoeff,Params)
% close all
figure(1);
set(gcf,'paperpositionmode','auto','position',[.1 .2 .5 .5])
% Figure a
Fig1 = subplot(2,1,1);
% PositionFig1 = get(Fig1,'position'); 
% subplot('position',[0.15 PositionFig1(2) 0.8 0.75/2])

loglog(qs, Params(:,3), 'ko');
hold on;
plot([1 1]*qs(nMin),[1e-3 1e5], '--','LineWidth',2)
plot([1 1]*qs(nMax),[1e-3 1e5], '--','LineWidth',2)
loglog(qs,1/DiffusionCoeff./(qs).^2 , 'r','LineWidth',2);
loglog(qs,(1/.5)./(qs).^2 , 'k','LineWidth',2);

ylabel('$\tau_d$ {(s)}','interpreter', 'latex','fontsize',18)
xlim([min(qs)*.8 max(qs)*1.2])
ylim([min(Params(:,3))*.5 max(Params(:,3))*2])

set(gca,'xticklabel',[])
set(gca,'FontSize',18)
set(gca,'fontname','times')

% Figure b
Fig2 = subplot(2,1,2);
% PositionFig2 = get(Fig2,'position'); 
% subplot('position',[0.15 PositionFig2(2)+0.09 0.8 0.75/2]) %dimension in the window

loglog(qs, Params(:,1), 'ks');hold on;
loglog(qs, Params(:,2), 'kd');
plot([1 1]*qs(nMin),[1e-3 1e5], '--','LineWidth',2)
plot([1 1]*qs(nMax),[1e-3 1e5], '--','LineWidth',2)

xlim([min(qs)*.8 max(qs)*1.2])
ylim([min(Params(:,1))*.5 max(Params(:,1))*2])

xlabel('$q$ {($\mu$m$^{-1}$)}','interpreter', 'latex','fontsize',18)
ylabel('$A, B$ {(a.u.)}','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18,'fontname','times')

ISF = 1-((DDM'-Params(:,2)*ones(1,length(dt)))./(Params(:,1)*ones(1,length(dt))))';

ISF_fit=0*ISF;
for i=1:length(Params)
ISF_fit(:,i) = exp(-dt/Params(i,3));
end


figure(2);
% clf
set(gcf,'paperpositionmode','auto','position',[.1 .2 .5 .5])
skip=1;
couleur = jet(length(nMin:skip:nMax));

xlabel('$\Delta t$ {(s)}','interpreter', 'latex','fontsize',18)
ylabel('$f(q, \Delta t)$','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18)
semilogx(dt,ISF(:,nMin:skip:nMax),'o'); hold on;
semilogx(dt,ISF_fit(:,nMin:skip:nMax))
% semilogx(dt,ISF(:,nMin:skip:nMax),'o'); hold on;
% semilogx(dt,ISF_fit(:,nMin:skip:nMax))

h = colorbar('location','West');
caxis([qs(nMin)  qs(nMax)]*10^3)
set(h,'Position',[0.85 0.81 .02 .1], 'YTick',10.^(round(log10(1000*min(qs))):round(log10(1000*max(qs)))))
locate = get(h,'title');
PositionTitre = get(locate, 'position');
set(locate,'pos',[-2.9 1.51 1],'string','$q$ ($\mu$m$^{-1}$)','interpreter', 'latex',  'FontSize',15,'fontname','times','Rotation',90);

xlim([min(dt)*.5 max(dt)*2])
ylim([-0.1 1.1])
set(gca,'XTick',10.^(-3:5)) 

set(gca,'fontname','times')
set(h,'FontSize',12,'fontname','times')
xlabel('$\Delta t $ {(s)}','interpreter', 'latex','fontsize',18)
ylabel('$f(q, \Delta t)$','interpreter', 'latex','fontsize',18)
end
