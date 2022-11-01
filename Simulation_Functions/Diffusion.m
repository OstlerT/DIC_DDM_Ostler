
function Steps_Diffusion = Diffusion(D_cartesian,time,deltat,NParticles,Pixelsize)
D_pixel=D_cartesian/(Pixelsize^2); %convert units to pixels from um
% Diff_Lengths = abs(normrnd(0,sqrt(2*D*deltat) , [NParticles 1 time-1]));
% Diff_angles = 2*pi*rand([NParticles 1 time-1]);
% Steps_Diffusion_X= Diff_Lengths.*cos(Diff_angles);
% Steps_Diffusion_Y= Diff_Lengths.*sin(Diff_angles);
Steps_Diffusion_X=normrnd(0,sqrt(2*D_pixel*deltat) , [NParticles 1 time-1]);
Steps_Diffusion_Y=normrnd(0,sqrt(2*D_pixel*deltat) , [NParticles 1 time-1]);

Steps_Diffusion = cumsum(cat(2,Steps_Diffusion_X,Steps_Diffusion_Y),3);
