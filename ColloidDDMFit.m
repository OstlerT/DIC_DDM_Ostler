clear all
close all
%To generate data fresh, datagen=1
%or to load existing data (recommended), datagen=0
datagen=0
if datagen==1
    clear all


    projectpath = genpath(pwd);
    addpath(projectpath);

    G_num = 3;
    %%%file name is the location AND name of the image stack you want to
    %%%process
    fileName = ['E:\Emily_Data\G',num2str(G_num),'_Final_Negative_Removed\G',num2str(G_num),'_Final_Negative_Removed.tif'];
    %%Foldersave is the name of the folder you want the output to be
    %%stored. If it doesnt exist, it gets made.
    FolderSave=['DIC_DDM_Ostler\ColloidData\']; %Folder where you want to save the OUTPUT
    mkdir(FolderSave); %create folder to store output data
    %%%%%%%%%%%%%%%Determine loading parameters%%%%%%%%
    fp = fopen(fileName , 'rb');
    info = imfinfo(fileName);
    numFramesStr = regexp(info.ImageDescription, 'images=(\d*)', 'tokens');


    NbImage = 3110;
    framenum=3;
    imData=cell(1,framenum);
    ftell_store = zeros(1,3);
    for cnt = 1:framenum
        ftell_store(cnt) = ftell(fp);
        imData{cnt} = fread(fp, [info.Width info.Height], 'uint8', 0, 'ieee-be')';
    end
    x=diff(ftell_store);
    if x(1) ~= x(2)
        error('The image stack is made up of images of different size')
    end
    fclose(fp);

    %%Meta Parameters
    PixelSize = 322.5; %Conversion 1 pixel -> nm
    Frequency = 7680/900; %Number of images/Timeframe

    %% Acquisition parameters obtained from INPUT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    MaxNCouples=700; %Analysis parameters: average is performed over MaxNCouples images at most for each dt, MaxNCouples=<NbImage-1
    ImageSize = info.Height; %The length of one dimension of the square image, Preferably is also the radius of the embryo.
    qs = ((2*pi)/(ImageSize*PixelSize))*(1:ImageSize/2); %generate wavevector q in physical units (nm-^1)


    % Generate log scale for the times idt, 15 points per decade
    Base = 10^(1/12); %power scale to get dt
    Ndt = floor(log(NbImage)/log(Base)); %max power for NbImage
    idts=unique(floor(Base.^(1:Ndt)));

    % Initialize DDM matrix and dt
    DDM = zeros(length(idts), ImageSize/2); %Image size must be an even number in order for the code to run
    dt = zeros(1,length(idts));


    %DDM algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Loops over each individual trajectory
    %loop on stack 1 and 2

    dt(1,:) = idts/Frequency; %dt in real time
    w=blackman(ImageSize)*blackman(ImageSize)';



    parfor p = 1:length(idts)
        idt = idts(p); %real time interval between 2 images
        disp([num2str(100*(length(idts)-p)/(length(idts))),'% remaining']);%Doesnt necessarily increment in order, but gives a guide on how much work is left
        AvgFFT  = zeros(ImageSize,ImageSize);
        % Spread initial times over the available range
        Increment = max(floor((NbImage - idt)/MaxNCouples), 1);
        if (NbImage - idt)< MaxNCouples
            InitialTimes=1:(NbImage-idt);
        else
            Pairlist=randperm(NbImage - idt);
            InitialTimes=Pairlist(1:MaxNCouples);
        end


        % loop on t
        for t_temp=1:length(InitialTimes)
            t=InitialTimes(t_temp);
            %             disp(t)
            fp = fopen(fileName , 'rb');

            fseek(fp, info.StripOffsets+x(2)*(t+idt-1), 'bof');
            Im2 = fread(fp, [info.Width info.Height], 'uint8', 0, 'ieee-be')';
            Im2 = Im2(:,200:1223)/255;
            fseek(fp, info.StripOffsets+x(2)*(t-1), 'bof');
            Im1 = fread(fp, [info.Width info.Height], 'uint8', 0, 'ieee-be')';
            Im1 = Im1(:,200:1223)/255;
            fclose(fp);
            AvgFFT = AvgFFT + abs(fftshift(fft2(w.*(Im2-Im1)))).^2;% FFT of differences squared, and add up

        end


        AvgFFT = AvgFFT./length(InitialTimes); % Divide by the number of image differences to get the average

        % Replace values at q = 0 by the mean of its neighbour
        AvgFFT(:,(ImageSize/2)+1) = 1/2*(AvgFFT(:,ImageSize/2)+AvgFFT(:,(ImageSize/2)+2)).*ones(ImageSize,1);

        [DDM(p, :), VecQ] =  radialavg(AvgFFT,ImageSize/2); % Radial average of the FFT


    end


    DDM(:,1) = [];
    qs(1)=[];


    save([FolderSave,'\DDM.mat'],'NbImage', 'ImageSize', 'Frequency', 'PixelSize', 'qs', 'dt', 'DDM');

    dtLimit = 30; %(max value = length(dt))

    %% Fit function for the DDM matrix
    MatrixFit=zeros(ImageSize/2-1,dtLimit); %initialization
    ISF_Fit=MatrixFit;%initialization
    Noise = mean(DDM(:,ImageSize/2-1)); % Noise floor
    %fit function for the Bacteria (ref: )
    Pv = @(a,xdata) ( (a(6)+1)*a(5)./(a(6)*xdata) )...
        .* sin( a(6)*atan(xdata /(a(5)*(a(6)+1))))...
        ./( (1+ (xdata/(a(5)*(a(6)+1))).^2).^(a(6)/2));
    ISF_Bact = @(a,xdata) exp(-xdata/a(3)).*((1-a(4))+a(4)*Pv(a, xdata));
    FitBacteria = @(a,xdata) log(a(1)*(1-ISF_Bact(a, xdata))+a(2));
    %fit function for the colloids
    FitColloid = @(a,xdata) log(a(1)*(1-exp(-xdata/a(3)))+a(2));


    Params=zeros(ImageSize/2-1,3); % initialization

    for Qinter =1:ImageSize/2-1
        %fit parameters initialization, lower and upper boundaries
        Params0 = [max(DDM(1:dtLimit,Qinter))*2,...
            Noise,...
            1];
        lb = [(max(DDM(1:dtLimit,Qinter))-Noise)*0.8, Noise*0.8,0];
        ub = [(max(DDM(1:dtLimit,Qinter))-Noise)*1.2, Noise*1.2,1000];
        Params(Qinter, :) = lsqcurvefit(FitColloid,Params0,dt(1:dtLimit),log(DDM(1:dtLimit,Qinter)'),lb,ub); %fit
        MatrixFit(Qinter,:) = exp(FitColloid(Params(Qinter,:),dt(1:dtLimit))); % DDM matrix fit
        ISF_Fit(Qinter,:)= exp(-(dt(1:dtLimit)./Params(Qinter,3))); % ISF fit
    end



    qMin = 1.1; %um^-1
    qMax = 3.3; %um^-1

    nMin= find(1000*qs < qMin, 1, 'last' );
    nMax= find(1000*qs < qMax, 1, 'last' );

    Velocity=0;
    DiffusionCoeff=0;

    %fit in log scale
    x0Diff = 1;
    FitDiffusion = @(a,xdata) -2*xdata + a;
    xDif = lsqcurvefit( FitDiffusion, x0Diff,log10(qs(nMin:nMax)),log10(Params(nMin:nMax,3))' );
    DiffusionCoeff = 10^(-xDif)*1e-6; %um2/s
    save([FolderSave,'\DDMFit.mat'], 'Params','FitColloid','FitBacteria','ISF_Bact','Pv','MatrixFit','ISF_Fit','dtLimit');
    save([FolderSave,'\DDMFitTau.mat'], 'Params','FitColloid','FitBacteria','ISF_Bact','Pv','MatrixFit','ISF_Fit','dtLimit');

else
    %Replace below with location of saved colloid dataset
    FolderSave=['F:\DIC_Paper_Github\DIC_DDM_Ostler\ColloidData\']; %Folder where you want to save the OUTPUT

    cd(FolderSave)
    load([FolderSave,'\DDM.mat']);
    load([FolderSave,'\DDMFit.mat']);
    load([FolderSave,'\DDMFitTau.mat']);

end

%%
% Figure 7


figure(3);
clf
set(gcf,'paperpositionmode','auto','position',[0.1 0.1 0.6 0.6])
% Figure a
Fig1 = subplot(2,1,1);
PositionFig1 = get(Fig1,'position');
subplot('position',[0.15 PositionFig1(2) 0.8 0.75/2])

loglog(qs*10^3, Params(:,3), 'ko');
hold on;
loglog(qs*1e3,1/DiffusionCoeff./(qs*1e3).^2 , 'r','LineWidth',2);
%loglog(qs*1e3,1/(.667)./(qs*1e3).^2,'-m')
%     loglog(qs*1e3,1/0.830./(qs*1e3).^2 , '--k','LineWidth',2);
plot([1 1]*qs(nMin)*1e3,[1e-3 1e5], '--','LineWidth',2,'Color',[0 0 .8])
plot([1 1]*qs(nMax)*1e3,[1e-3 1e5], '--','LineWidth',2,'Color',[0 0 .8])
legend('Fitted Values','Linear regression','interpreter','latex','autoupdate','off')

ylabel('$\tau_{\mathcal{D}_I}$ {(s)}','interpreter', 'latex','fontsize',18)
xlim(1000*[min(qs)*.8 max(qs)*1.2])
ylim([1e-3 1e5])

set(gca,'xticklabel',[])
set(gca,'FontSize',18)
set(gca,'fontname','times')

% Figure b
Fig2 = subplot(2,1,2);
PositionFig2 = get(Fig2,'position');
subplot('position',[0.15 PositionFig2(2)+0.09 0.8 0.75/2]) %dimension in the window

loglog(qs*10^3, Params(:,1), 'ks');hold on;
loglog(qs*10^3, Params(:,2), 'kd');
legend('$A_I(q)$','$B_I(q)$','interpreter','latex','autoupdate','off','location','southwest')
plot([1 1]*qs(nMin)*1e3,[1e-6 1e3], '--','LineWidth',2,'Color',[0 0 .8])
plot([1 1]*qs(nMax)*1e3,[1e-6 1e3], '--','LineWidth',2,'Color',[0 0 .8])

xlim(1000*[min(qs)*.8 max(qs)*1.2])
ylim([1e-6 1e3])

xlabel('$q$ $\left(\mu \textrm{m}^{-1}\right)$','interpreter', 'latex','fontsize',18)
ylabel('$A, B$ {(a.u.)}','interpreter', 'latex','fontsize',18)
set(gca,'FontSize',18,'fontname','times')
