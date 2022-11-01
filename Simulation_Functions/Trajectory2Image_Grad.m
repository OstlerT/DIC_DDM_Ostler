
function [ImageStack] = Trajectory2Image_Grad(Particle_Locations,WindowSize,domainlength,time,sigma,alpha)
NParticles=size(Particle_Locations,1);
ImageStack = zeros(WindowSize,WindowSize,time);
% ImageStackfft = zeros(WindowSize,WindowSize,time);
% w=blackman(size(ImageStack,1))*blackman(size(ImageStack,1))';
% pwing = floor((particlesize-1)/2);
% Particle_Locations(Particle_Locations<= pwing)=Particle_Locations(Particle_Locations<=pwing)+pwing;
R=repmat(1:domainlength,NParticles,1);
if domainlength==WindowSize
    startind=1;endind=domainlength;
else
    
startind = floor((domainlength-WindowSize)/2);  endind = startind+WindowSize-1;
end
R=R(:,startind:endind,:);
for t=1:time
%      x0 = squeeze(Particle_Locations(p,1,:));
%      y0 = squeeze(Particle_Locations(p,2,:));
    ImageStack(:,:,t) = (alpha*(exp(-((R-repmat(squeeze(Particle_Locations(:,1,t)),1,WindowSize)).^2)/(2*sigma^2))' * (exp(-((R-repmat(squeeze(Particle_Locations(:,2,t)),1,WindowSize)).^2)/(2*sigma^2)))));
%     ImageStackfft(:,:,t) = fftshift(fft2(w.*ImageStack(:,:,t)));

end

