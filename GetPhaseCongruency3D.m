function [pcSum,Vesselness,Neuritenees] = GetPhaseCongruency3D(im,sigmaOnf,k,noiseMethod)
%% Phase
nscale          = 4;%3     % Number of wavelet scales.    
norient         = 20;%10    % Number of filter orientations.
minWaveLength   = 3; %3    % Wavelength of smallest scale filter.    
mult            = 2.1;   % Scaling factor between successive filters.
g               = 10;%10    % Controls the sharpness of the transition in
                         % the sigmoid function used to weight phase
                         % congruency for frequency spread.                
cutOff          = 0.5;%0.5   % The fractional measure of frequency spread
                         % below which phase congruency values get penalized.      
                         
if nargin == 1 % 如果输入参数只有一个，就是只有图片，就采用默认参数
    sigmaOnf        = 0.01;%0.45  % Ratio of the standard deviation of the
                             % Gaussian describing the log Gabor filter's
                             % transfer function in the frequency domain
                             % to the filter center frequency.        
    k               = 5;%5     % No of standard deviations of the noise
                             % energy beyond the mean at which we set the
                             % noise threshold point. 
    noiseMethod     = 0.2;%-1    % Choice of noise compensation method. 
end
    
%% Phase Congruency 3D
[PC,EO,T,pcSum,orients] = PhaseCongruency3D(im,nscale,norient,minWaveLength,mult,sigmaOnf,k,cutOff,g,noiseMethod);
%%
pcSum = Normalize(pcSum);
imq = zeros(size(im,1),size(im,2),size(im,3),size(PC,1));
for i=1:size(PC,1)
    imq(:,:,:,i) = PC{i,1};
end
imq = Normalize(imq);
%% Phase Vesselness 3D
alpha = 0.5; beta = 0.5; c = 15;
Vesselness = PhaseVesselness3D(imq,orients,alpha,beta,c);

%% Phase Neuriteness 3D
sigma = [2 2 2] ; 
Neuritenees = NeuriteneesFilter3D(imq,sigma);

end