function [DDM, NonRadAvgDDM,dt, qs] = DDMAlgorithm(Image,PixelSize,MaxNCouples,Frequency)


NbImage = size(Image,3);%total number of images you want to analyse in the stacks
Base = 10^(1/12); %power scale to get dt
Ndt = floor(log(NbImage)/log(Base)); %max power for NbImage
idts=unique(floor(Base.^(1:Ndt)));
ImageSize = size(Image,1); %The length of one dimension of the square image, Preferably is also the radius of the embryo.
qs = ((2*pi)/(ImageSize*PixelSize))*(1:ImageSize/2); %generate wavevector q in physical units (nm-^1)
w=blackman(size(Image,1))*blackman(size(Image,1))'; %Requires signal processing toolbox
Image=w.*Image;%Apply blackman filter to image stack
DDM = zeros(length(idts), ImageSize/2); %Image size must be an even number in order for the code to run
dt = zeros(1,length(idts));
NonRadAvgDDM=zeros(size(Image));
%DDM algorithm


dt(1,:) = idts/Frequency; %dt in real time
for p = 1:length(idts)
    idt = idts(p); 
    
    AvgFFT  = zeros(ImageSize,ImageSize);
    % Spread initial times over the available range
    if (NbImage - idt)< MaxNCouples
        InitialTimes=1:(NbImage-idt);
    else
    Pairlist=randperm(NbImage - idt);
    InitialTimes=Pairlist(1:MaxNCouples);
    end
    for t=InitialTimes
        AvgFFT = AvgFFT + abs(fftshift(fft2((Image(:,:,t+idt)-Image(:,:,t))))).^2;% FFT of differences squared, and add up
        
    end
    
    AvgFFT = AvgFFT./length(InitialTimes); % Divide by the number of image differences to get the average
%     disp(length(InitialTimes));
    % Replace values at q = 0 by the mean of its neighbour (A modifier)
    AvgFFT(:,(ImageSize/2)+1) = 1/2*(AvgFFT(:,ImageSize/2)+AvgFFT(:,(ImageSize/2)+2)).*ones(ImageSize,1);
    NonRadAvgDDM(:,:,p) = AvgFFT;
    [DDM(p, :), VecQ] =  radialavg(AvgFFT,ImageSize/2);
end

% Remove the first q 
DDM(:,1) = [];
qs(1)=[];