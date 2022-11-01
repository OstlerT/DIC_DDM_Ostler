
function Steps_Velocity = DirectedAdvection(v,phi_v,time,deltat,NParticles)

Steps_Velocity_X = deltat*v*cos(phi_v)*ones(NParticles, 1, time-1);
Steps_Velocity_Y = deltat*v*sin(phi_v)*ones(NParticles, 1, time-1);
% Steps_Velocity = cat(2,Steps_Velocity_X,Steps_Velocity_Y);
Steps_Velocity = cumsum(cat(2,Steps_Velocity_X,Steps_Velocity_Y),3);
