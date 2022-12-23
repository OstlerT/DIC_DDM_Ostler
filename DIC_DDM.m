

clear all
%Set random seed to default for reproducibility
rng('default')
seed = rng
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
Domain_Length = 600;

%
% Generating Particle locations
Particle_Locations = repmat(Domain_Length*rand(Num_Particles,2),1,1,Time); %Randomly place particles
Particle_Locations(:,:,2:end) = Particle_Locations(:,:,2:end)+Diffusion(D,Time,deltat,Num_Particles,Pixel_Size); %Adds diffusive motion for all timepoints after the first
s = Shear_Length * [sin(Shear_Direction)*ones(Num_Particles,1,Time) cos(Shear_Direction)*ones(Num_Particles,1,Time)];
X = Trajectory2Image_Grad(mod(Particle_Locations,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Y = Trajectory2Image_Grad(mod(Particle_Locations+s,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Particle_Locations_shift = Particle_Locations+s;
I=Y-X;
MaxNCouples=2000;
%Generate DDM matrices
[DDM0, NonRadAvg0, dt, qs]= DDMAlgorithm(X/255,Pixel_Size,MaxNCouples,Frequency);
 [DDM1, NonRadAvg1, ~,~]= DDMAlgorithm(I/255,Pixel_Size,MaxNCouples,Frequency);
%Define fitting interval
qMin = .5;%.15; %um^-1
qMax = 1.2; %um^-1
%Find fitting interval indices from list of frequencies
nMin= find(qs < qMin, 1, 'last' );
nMax= find(qs < qMax, 1, 'last' );
qsmax=150;
%Perform fitting
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
imshow(X(:,:,1),[])
title('$X$')
colorbar
nexttile
imshow(I(:,:,1),[ ])
colorbar
title('$I$')




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
    semilogy(qs(1:4:qsmax),2*(1-besselj(0,Shear_Length*qs(1:4:qsmax))).*DDM0(t,1:4:qsmax),'*','MarkerSize',10,'Color',[0 0 0])
% xline(qMin, '--','LineWidth',3) 
% xline(qMax, '--','LineWidth',3) 
       legend('$\mathcal{D}_I$','$2(1-J_0(q s)) \mathcal{D}_X$','interpreter','latex','AutoUpdate','off')

       semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[.5 1-dt(t)/dt(end_time) 1-dt(t)/dt(end_time)])
    hold on
    semilogy(qs(1:4:qsmax),2*(1-besselj(0,Shear_Length*qs(1:4:qsmax))).*DDM0(t,1:4:qsmax),'*','MarkerSize',10,'Color',[.5 1-dt(t)/dt(end_time) 1-dt(t)/dt(end_time)])
    else
        semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[.5 1-dt(t)/dt(end_time) 1-dt(t)/dt(end_time)])
    hold on
    semilogy(qs(1:4:qsmax),2*(1-besselj(0,Shear_Length*qs(1:4:qsmax))).*DDM0(t,1:4:qsmax),'*','MarkerSize',10,'Color',[.5 1-dt(t)/dt(end_time) 1-dt(t)/dt(end_time)])
    end
        
% plot(1000*qs,Cross_Corr_RadAvg(t,:),'-b')
% plot(1000*qs,.5*(Params0(:,1)).*(x1y1' - x1y2(:,t)) ,'--r');

end



c=colorbar;
mapmesh=50;
map = [.5*ones(mapmesh,1) 1-linspace(0,1,mapmesh)' 1-linspace(0,1,mapmesh)'];
colormap(map)
caxis([dt(1) dt(14)])
c.Label.String = '$\Delta t$ (s)';
c.Label.Interpreter = 'latex';


hold on
xlabel('$q$ $\left(\mu \textrm{m}^{-1}\right)$','interpreter','latex')
set(gca,'FontSize',25)



%%
Datagen = 0; %1 to generate new data, 0 to load
if Datagen==1
%Generate data to compare fitting accuracy
clear all
s = rng('shuffle');
PC_ID = GetPC_ID();
NumSims=150;
Dfit=zeros(NumSims,2);
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

%
% Generating Particle locations
Domain_Length = 600;
%

for i=1:NumSims
Particle_Locations = repmat(Domain_Length*rand(Num_Particles,2),1,1,Time); %Randomly place particles
Particle_Locations(:,:,2:end) = Particle_Locations(:,:,2:end)+Diffusion(D,Time,deltat,Num_Particles,Pixel_Size); %Adds diffusive motion for all timepoints after the first
s = Shear_Length * [sin(Shear_Direction)*ones(Num_Particles,1,Time) cos(Shear_Direction)*ones(Num_Particles,1,Time)];
X = Trajectory2Image_Grad(mod(Particle_Locations,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Y = Trajectory2Image_Grad(mod(Particle_Locations+s,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Particle_Locations_shift = Particle_Locations+s;
I=Y-X;
MaxNCouples=2000;
%
[DDM0, NonRadAvg0, dt, qs]= DDMAlgorithm(X/255,Pixel_Size,MaxNCouples,Frequency);
 [DDM1, NonRadAvg1, ~,~]= DDMAlgorithm(I/255,Pixel_Size,MaxNCouples,Frequency);

 qMin = .5;%.15; %um^-1
qMax = 1.2; %um^-1
%
nMin= find(qs < qMin, 1, 'last' );
nMax= find(qs < qMax, 1, 'last' );
% figure
% qsmax=150;
%
[DiffCoeff0, Params0, ~, ~] = DDM_fit(DDM0,dt,length(dt), qs,qMin,qMax);

 [DiffCoeff1, Params1, nMin, nMax] = DDM_fit(DDM1,dt,length(dt), qs,qMin,qMax);
Dfit(i,:)=[DiffCoeff0 DiffCoeff1];
end

save('Boxplot_data.mat','Dfit','s')
else
%
load('Boxplot_data.mat')
end
%%
% Figure 6a
clf
boxplot(Dfit,'Labels',{'$\mathcal{D}_X$','$\mathcal{D}_I$'}); hold on;
set(gca,'TickLabelInterpreter','latex');
scatter(ones(150,1),Dfit(:,1)); scatter(2*ones(150,1),Dfit(:,2 ))
max(abs(Dfit(:,1)-Dfit(:,2)))
[mean(Dfit(:,1)) mean(Dfit(:,2)) ]
[std(Dfit(:,1)) std(Dfit(:,2))]
ylabel('$D$ $\left(\mu \textrm{m}^2 \textrm{s}^{-1}\right)$','interpreter','latex')

%6b
clf
histogram(Dfit(:,1)-Dfit(:,2),'NumBins',14)
xlabel('$D_X-D_I$ $\left(\mu \textrm{m}^2 \textrm{s}^{-1}\right)$','interpreter','latex')
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
Time = 4000;
Frequency = 8;
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

s = Shear_Length * [sin(Shear_Direction)*ones(Num_Particles,1,Time) cos(Shear_Direction)*ones(Num_Particles,1,Time)];
X = Trajectory2Image_Grad(mod(Particle_Locations,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Y = Trajectory2Image_Grad(mod(Particle_Locations+s,Domain_Length),WindowSize,Domain_Length,Time,sigma,alpha);
Particle_Locations_shift = Particle_Locations+s;
I=Y-X;
MaxNCouples=800;
%
[DDM0, NonRadAvg0, dt, qs]= DDMAlgorithm(X/255,Pixel_Size,MaxNCouples,Frequency);
 [DDM1, NonRadAvg1, ~,~]= DDMAlgorithm(I/255,Pixel_Size,MaxNCouples,Frequency);
%%
%Figure 7
 qMin = .5;%.15; %um^-1
qMax = 1.2; %um^-1
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
    semilogy(qs(1:4:qsmax),2*(1-besselj(0,Shear_Length*qs(1:4:qsmax))).*DDM0(t,1:4:qsmax),'*','MarkerSize',10,'Color',[0 0 0])
% xline(qMin, '--','LineWidth',3) 
% xline(qMax, '--','LineWidth',3) 
       legend('$\mathcal{D}_I$','$2(1-J_0(q s)) \mathcal{D}_X$','interpreter','latex','AutoUpdate','off')

    semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[.5 1-dt(t)/dt(14) 1-dt(t)/dt(14)])
    hold on
    semilogy(qs(1:4:qsmax),2*(1-besselj(0,Shear_Length*qs(1:4:qsmax))).*DDM0(t,1:4:qsmax),'*','MarkerSize',10,'Color',[.5 1-dt(t)/dt(14) 1-dt(t)/dt(14)])
    else
    semilogy(qs(1:qsmax),DDM1(t,1:qsmax),'-','Color',[.5 1-dt(t)/dt(14) 1-dt(t)/dt(14)])
    hold on
    semilogy(qs(1:4:qsmax),2*(1-besselj(0,Shear_Length*qs(1:4:qsmax))).*DDM0(t,1:4:qsmax),'*','MarkerSize',10,'Color',[.5 1-dt(t)/dt(14) 1-dt(t)/dt(14)])
    end
        
% plot(1000*qs,Cross_Corr_RadAvg(t,:),'-b')
% plot(1000*qs,.5*(Params0(:,1)).*(x1y1' - x1y2(:,t)) ,'--r');

end



c=colorbar;
mapmesh=50;
map = [.5*ones(mapmesh,1) 1-linspace(0,1,mapmesh)' 1-linspace(0,1,mapmesh)'];
colormap(map)
caxis([dt(1) dt(14)])
c.Label.String = '$\Delta t$ (s)';
c.Label.Interpreter = 'latex';
hold on
xlabel('$q$ $\left(\mu \textrm{m}^{-1}\right)$','interpreter','latex')
set(gca,'FontSize',25)


% %

