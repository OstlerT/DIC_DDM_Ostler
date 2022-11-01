

clear all
PC_ID = GetPC_ID();
rng('default')
s = rng
%%%%Meta-Parameters%%%
WindowSize = 490;
Shear_Length = 1;
Shear_Direction = 0; %radians

Num_Particles = 150;
Time = 4000;
Frequency = 10;
Pixel_Size = 1; %Conversion to um
deltat = 1/Frequency;
D=.5;
sigma = 3; %Gaussian Beam variance
alpha=50;%Gaussian Beam brightness
%%%%Image generation params%%%%
highres=1;
%
% Generating Particle locations
Domain_Length = 600;

Particle_Locations = repmat(Domain_Length*rand(Num_Particles,2),1,1,Time); %Randomly place particles
Particle_Locations(:,:,2:end) = Particle_Locations(:,:,2:end)+Diffusion(D,Time,deltat,Num_Particles,Pixel_Size); %Adds diffusive motion for all timepoints after the first
x_0 = Shear_Length * [sin(Shear_Direction)*ones(Num_Particles,1,Time) cos(Shear_Direction)*ones(Num_Particles,1,Time)];
Image0 = Trajectory2Image_Grad(mod(Particle_Locations,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Image0_shift = Trajectory2Image_Grad(mod(Particle_Locations+x_0,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Particle_Locations_shift = Particle_Locations+x_0;
Image1=Image0_shift-Image0;
MaxNCouples=2000;
%
[DDM0, NonRadAvg0, dt, qs]= DDMAlgorithm(Image0/255,Pixel_Size,MaxNCouples,Frequency);


 [DDM1, NonRadAvg1, ~,~]= DDMAlgorithm(Image1/255,Pixel_Size,MaxNCouples,Frequency);
%
 qMin = .3;%.15; %um^-1
qMax = 1.1; %um^-1
%
nMin= find(qs < qMin, 1, 'last' );
nMax= find(qs < qMax, 1, 'last' );
qsmax=150;
%
[DiffCoeff0, Params0, ~, ~] = DDM_fit(DDM0,dt,length(dt)-2, qs,qMin,qMax);

[DiffCoeff1, Params1, nMin, nMax] = DDM_fit(DDM1,dt,length(dt)-2, qs,qMin,qMax);

%%
% Visualise fitting paramaters for DDM matrix on X
clf
DDMPlots(DDM0,qs,dt,nMin,nMax,DiffCoeff0,Params0)
%%
% Visualise fitting paramaters for DDM matrix on I

DDMPlots(DDM1,qs,dt,nMin,nMax,DiffCoeff1,Params1)

%%
%Visualising one frame of the simulations for X and I
%Figure 4
t=tiledlayout(1,2);
nexttile
imshow(Image0(:,:,1),[])
title('$X$')
colorbar
nexttile
imshow(Image1(:,:,1),[ ])
colorbar
title('$I$')
if highres==1
export_fig(['C:\Users\C1602597\Dropbox\DIC\','Images','\FirstFrame.png'], '-r300')
else
export_fig(['C:\Users\C1602597\Dropbox\DIC\','Images_lowres','\FirstFrame.png'])
end

%%
%Generating 'Kidney bean'diagrams for each matrix
clf
t=1;
max1 = (max(NonRadAvg1(:,:,1),[],'all'));
max0= (max(NonRadAvg0(:,:,1),[],'all'));
titleset={'$\mathcal{D}_{X}$','$\mathcal{D}_{I}$'}
fig=figure(1);

for i=0:1
    subplot(1,2,i+1)
    if i==0
contourf([-flip(qs(1:100)) -0.000001 0.000001 qs(1:100)], [-flip(qs(1:100)) -0.000001 0.000001 qs(1:100)], squeeze(NonRadAvg0(145:end-144,145:end-144,t)))
    colorbar
caxis([0 max0])

    else
        contourf([-flip(qs(1:100)) -0.000001 0.000001 qs(1:100)], [-flip(qs(1:100)) -0.000001 0.000001 qs(1:100)], squeeze(NonRadAvg1(145:end-144,145:end-144,t)))
        colorbar
caxis([0 max1])
    end

title(titleset{i+1},'interpreter','latex')
xlabel('$q_x$','interpreter','latex')
ylabel('$q_y$','interpreter','latex')
end
h = axes(fig,'visible','off'); 

axis square

set(gca,'FontSize',25)



%%
% Figure 5
clf
close all
figure('position',[0.1 0.1 0.6 0.6])
end_time=14;
for t=1:2:end_time
    if t==1
%     plot(1000*qs, DDM0(t,:),'-.m')
    semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[0 0 0])
    hold on
    semilogy(qs(1:2:qsmax),2*(1-besselj(0,Shear_Length*qs(1:2:qsmax))).*DDM0(t,1:2:qsmax),'x','Color',[0 0 0])
% xline(qMin, '--','LineWidth',3) 
% xline(qMax, '--','LineWidth',3) 
       legend('$\mathcal{D}_I$','$2(1-J_0(q x_0)) \mathcal{D}_X$','interpreter','latex','AutoUpdate','off')

       semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[.5 1-dt(t)/dt(end_time) 1-dt(t)/dt(end_time)])
    hold on
    semilogy(qs(1:2:qsmax),2*(1-besselj(0,Shear_Length*qs(1:2:qsmax))).*DDM0(t,1:2:qsmax),'x','MarkerSize',10,'Color',[.5 1-dt(t)/dt(end_time) 1-dt(t)/dt(end_time)])
    else
        semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[.5 1-dt(t)/dt(end_time) 1-dt(t)/dt(end_time)])
    hold on
    semilogy(qs(1:2:qsmax),2*(1-besselj(0,Shear_Length*qs(1:2:qsmax))).*DDM0(t,1:2:qsmax),'x','MarkerSize',10,'Color',[.5 1-dt(t)/dt(end_time) 1-dt(t)/dt(end_time)])
    end
        
% plot(1000*qs,Cross_Corr_RadAvg(t,:),'-b')
% plot(1000*qs,.5*(Params0(:,1)).*(x1y1' - x1y2(:,t)) ,'--r');

end



c=colorbar;
mapmesh=50;
map = [.5*ones(mapmesh,1) 1-linspace(0,1,mapmesh)' 1-linspace(0,1,mapmesh)'];
colormap(map)
caxis([dt(1) dt(14)])
c.Label.String = '$\Delta t$, s';
c.Label.Interpreter = 'latex';


% xline(qs(nMin))
hold on
% xline(qs(nMax))
xlabel('$q, \mu \textrm{m}^{-1}$','interpreter','latex')
% ylabel('a.u','interpreter','latex')
set(gca,'FontSize',25)

if highres==1
export_fig(['C:\Users\C1602597\Dropbox\DIC\','Images','\DDMrelation.png'], '-r300')
else
export_fig(['C:\Users\C1602597\Dropbox\DIC\','Images_lowres','\DDMrelation.png'])
end



%%
%Generate data to compare fitting accuracy
clear all
rng('default');
PC_ID = GetPC_ID();
Dfit=zeros(50,2);
%%%%Meta-Parameters%%%
WindowSize = 490;
Shear_Length = 1;
Shear_Direction = 0; %radians

Num_Particles = 150;
Time = 4000;
Frequency = 10;
Pixel_Size = 1; %Conversion to um
deltat = 1/Frequency;
D=.5;
sigma = 3; %Gaussian Beam variance
alpha=50;%Gaussian Beam brightness
%%%%Image generation params%%%%
highres=1;
%
% Generating Particle locations
Domain_Length = 600;
%
for i=1:50
Particle_Locations = repmat(Domain_Length*rand(Num_Particles,2),1,1,Time); %Randomly place particles
Particle_Locations(:,:,2:end) = Particle_Locations(:,:,2:end)+Diffusion(D,Time,deltat,Num_Particles,Pixel_Size); %Adds diffusive motion for all timepoints after the first
x_0 = Shear_Length * [sin(Shear_Direction)*ones(Num_Particles,1,Time) cos(Shear_Direction)*ones(Num_Particles,1,Time)];
Image0 = Trajectory2Image_Grad(mod(Particle_Locations,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Image0_shift = Trajectory2Image_Grad(mod(Particle_Locations+x_0,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Particle_Locations_shift = Particle_Locations+x_0;
Image1=Image0_shift-Image0;
MaxNCouples=2000;
%
% [DDM0, NonRadAvgDDM0,dt, qs] = DDMAlgorithm_framegen(Particle_Locations,Pixel_Size,MaxNCouples,Frequency,alpha,sigma,WindowSize,Domain_Length,Num_Particles,Time)
[DDM0, NonRadAvg0, dt, qs]= DDMAlgorithm(Image0/255,Pixel_Size,MaxNCouples,Frequency);
 [DDM1, NonRadAvg1, ~,~]= DDMAlgorithm(Image1/255,Pixel_Size,MaxNCouples,Frequency);

 qMin = .5;%.15; %um^-1
qMax = 1.2; %um^-1
%
nMin= find(qs < qMin, 1, 'last' );
nMax= find(qs < qMax, 1, 'last' );
% figure
qsmax=150;
%
[DiffCoeff0, Params0, ~, ~] = DDM_fit(DDM0,dt,length(dt), qs,qMin,qMax);

 [DiffCoeff1, Params1, nMin, nMax] = DDM_fit(DDM1,dt,length(dt), qs,qMin,qMax);
Dfit(i,:)=[DiffCoeff0 DiffCoeff1];
end
save('Boxplot_data.mat','Dfit')
%%
load('Boxplot_data.mat')
%%
% Figure 6
clf
boxplot(Dfit,'Labels',{'$\mathcal{D}_X$','$\mathcal{D}_I$'}); hold on;
set(gca,'TickLabelInterpreter','latex');
scatter(ones(50,1),Dfit(:,1)); scatter(2*ones(50,1),Dfit(:,2 ))
max(abs(Dfit(:,1)-Dfit(:,2)))
[mean(Dfit(:,1)) mean(Dfit(:,2)) ]
[std(Dfit(:,1)) std(Dfit(:,2))]
ylabel('$D$, $\mu \textrm{m}^2 \textrm{s}^{-1}$','interpreter','latex')

export_fig(['C:\Users\C1602597\Dropbox\DIC\','Images','\DFitBoxPLots.png'], '-r300')
%%
clf
histogram(Dfit(:,1)-Dfit(:,2))
xlabel('$D_X-D_I$, $\mu \textrm{m}^2 \textrm{s}^{-1}$','interpreter','latex')
ylabel('Count')
export_fig(['C:\Users\C1602597\Dropbox\DIC\','Images','\DFitBoxDiff.png'], '-r300')

%%
%Advection-Diffusion simulations
clf

clear all
PC_ID = GetPC_ID();
%%%%Meta-Parameters%%%
WindowSize = 490;
Shear_Length = 1;
Shear_Direction = 0; %radians

Num_Particles = 150;
Time = 2000;
Frequency = 4;
Pixel_Size = 1; %Conversion to um
deltat = 1/Frequency;
D=.5;
v=.5;
phi_v=pi/2;
sigma = 3; %Gaussian Beam variance
alpha=50;%Gaussian Beam brightness
%%%%Image generation params%%%%
highres=1;
%
% Generating Particle locations
Domain_Length = 600;

Particle_Locations = repmat(Domain_Length*rand(Num_Particles,2),1,1,Time); %Randomly place particles
Particle_Locations(:,:,2:end) = Particle_Locations(:,:,2:end)+Diffusion(D,Time,deltat,Num_Particles,Pixel_Size); %Adds diffusive motion for all timepoints after the first
Particle_Locations(:,:,2:end) = Particle_Locations(:,:,2:end)+ DirectedAdvection(v,phi_v,Time,deltat,Num_Particles)

x_0 = Shear_Length * [sin(Shear_Direction)*ones(Num_Particles,1,Time) cos(Shear_Direction)*ones(Num_Particles,1,Time)];
Image0 = Trajectory2Image_Grad(mod(Particle_Locations,Domain_Length),WindowSize,Domain_Length,Time,Num_Particles,sigma,alpha);
Image0_shift = Trajectory2Image_Grad(mod(Particle_Locations+x_0,Domain_Length),WindowSize,Domain_Length,Time,Num_Particles,sigma,alpha);
Particle_Locations_shift = Particle_Locations+x_0;
Image1=Image0_shift-Image0;
MaxNCouples=800;
%
[DDM0, NonRadAvg0, dt, qs]= DDMAlgorithm(Image0/255,Pixel_Size,MaxNCouples,Frequency);
 [DDM1, NonRadAvg1, ~,~]= DDMAlgorithm(Image1/255,Pixel_Size,MaxNCouples,Frequency);
%%
%Figure 7
 qMin = .3;%.15; %um^-1
qMax = 1.1; %um^-1
%
nMin= find(qs < qMin, 1, 'last' );
nMax= find(qs < qMax, 1, 'last' );
qsmax=150;
clf
close all
figure('position',[0.1 0.1 0.6 0.6])

for t=1:2:14
    if t==1
%     plot(1000*qs, DDM0(t,:),'-.m')
    semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[0 0 0])
    hold on
    semilogy(qs(1:2:qsmax),2*(1-besselj(0,Shear_Length*qs(1:2:qsmax))).*DDM0(t,1:2:qsmax),'x','Color',[0 0 0])
% xline(qMin, '--','LineWidth',3) 
% xline(qMax, '--','LineWidth',3) 
       legend('$\mathcal{D}_I$','$2(1-J_0(q x_0)) \mathcal{D}_X$','interpreter','latex','AutoUpdate','off')

    semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[.5 1-dt(t)/dt(14) 1-dt(t)/dt(14)])
    hold on
    semilogy(qs(1:2:qsmax),2*(1-besselj(0,Shear_Length*qs(1:2:qsmax))).*DDM0(t,1:2:qsmax),'x','MarkerSize',10,'Color',[.5 1-dt(t)/dt(14) 1-dt(t)/dt(14)])
    else
    semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[.5 1-dt(t)/dt(14) 1-dt(t)/dt(14)])
    hold on
    semilogy(qs(1:2:qsmax),2*(1-besselj(0,Shear_Length*qs(1:2:qsmax))).*DDM0(t,1:2:qsmax),'x','MarkerSize',10,'Color',[.5 1-dt(t)/dt(14) 1-dt(t)/dt(14)])
    end
        
% plot(1000*qs,Cross_Corr_RadAvg(t,:),'-b')
% plot(1000*qs,.5*(Params0(:,1)).*(x1y1' - x1y2(:,t)) ,'--r');

end



c=colorbar;
mapmesh=50;
map = [.5*ones(mapmesh,1) 1-linspace(0,1,mapmesh)' 1-linspace(0,1,mapmesh)'];
colormap(map)
caxis([dt(1) dt(14)])
c.Label.String = '$\Delta t$, s';
c.Label.Interpreter = 'latex';


% xline(qs(nMin))
hold on
% xline(qs(nMax))
xlabel('$q, \mu \textrm{m}^{-1}$','interpreter','latex')
% ylabel('a.u','interpreter','latex')
set(gca,'FontSize',25)
if highres==1
export_fig(['C:\Users\C1602597\Dropbox\DIC\','Images','\DDMrelation_velocity.png'], '-r300')
else
export_fig(['C:\Users\C1602597\Dropbox\DIC\','Images_lowres','\DDMrelation_velocity.png'])
end
%%
function [DiffCoeff, Params, nMin, nMax] = DDM_fit(DDM,dt,dtLimit, qs,qMin,qMax)

opts = optimset('Display','off');
MatrixFit=zeros(length(qs),dtLimit); %initialization
ISF_Fit=MatrixFit;%initialization
Noise = mean(DDM(:,length(qs))); % Noise floor

FitColloid = @(a,xdata) (a(1)*(1-exp(-xdata/a(3)))+a(2));    


    Params=zeros(length(qs),3); % initialization
      
    for Qinter =1:length(qs)
        %fit parameters initialization, lower and upper boundaries
        Params0 = [max(DDM(1:dtLimit,Qinter))*2,...
                    Noise,...
                    1];
        lb = [0,0,0];
        ub = [Inf,Inf,Inf];  
        Params(Qinter, :) = lsqcurvefit(FitColloid,Params0,dt(1:dtLimit),(DDM(1:dtLimit,Qinter)'),lb,ub,opts); %fit
        MatrixFit(Qinter,:) = exp(FitColloid(Params(Qinter,:),dt(1:dtLimit))); % DDM matrix fit
        ISF_Fit(Qinter,:)= exp(-(dt(1:dtLimit)./Params(Qinter,3))); % ISF fit
    end

nMin= find(qs < qMin, 1, 'last' );
nMax= find(qs < qMax, 1, 'last' );
Velocity=0;
DiffusionCoeff=0;

% Fit colloids
% fit in log scale
x0Diff = 1;
FitDiffusion = @(a,xdata) -2*xdata + a;
xDif = lsqcurvefit( FitDiffusion, x0Diff,log10(qs(nMin:nMax)),log10(Params(nMin:nMax,3))',[],[],opts );
DiffCoeff = 10^(-xDif); %um2/s

end

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
