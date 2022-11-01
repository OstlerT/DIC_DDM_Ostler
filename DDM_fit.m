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
