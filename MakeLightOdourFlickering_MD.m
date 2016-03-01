% makes control paradigms so that the effective flow through the odor in the
% main air stream is uniformly distrubted, with mean $dil_mean and minimum
% $dil_mean - $dil_rande and maximum $dil_mean + $dil_range
% 
% this is the script used to make the control paradigms in the first mean
% shifted experiment. it has been slighlty modified: to a slower correlaton
% time (100ms), and now there is no step on to a flicker -- it always
% flickers. also, because the configutation has changed, some parts have
% been modified to talk to the MFCs correctly. 
%
% this specific script makes one odour stimulus, and combines it with
% differnet levels of light activation. 
%
% ControlParadigm Order (should match setup, and config. in Kontroller)
%
% 1. (AO)   to LED Dimmer (0-10V)
% 2. (AO)   to MFC500 (0-5V)
% 3. (DO)   to switch controlling main air @ 2L/min

clear ControlParadigm
ControlParadigm(1).Name = 'CleanAir';
ControlParadigm(1).Outputs = zeros(3,1e4);
ControlParadigm(1).Outputs(3,:) = 1;


dt  =1e-4;
T = 60;

%switching time
tc = .1; % 50ms is too fast for the MFCs to follow

nsteps= T/tc;
% noise = rand(nsteps,1);
% load frozen noise
load('noise.mat')

% define light levels
% light_mean = [0 logspace(log10(1.8),1,4)];
light_mean = [0];

% for MFC500
dil_min  = 0/100;
dil_max  = 20/100;

dil_mean = (dil_min +dil_max)/2;
V = 2000; %mL/min
MFC_Scale = 100; % 1V=100mL/min
min_flow = 5/200; % guesstimate from turn down ratio

for i = 1:length(light_mean)
    
    n = strcat('Odour+',oval(light_mean(i),2),'V Light');
    ControlParadigm2(i).Name = n;
    
    ControlParadigm2(i).Outputs= zeros(3,T/dt);
    
    % set to a random number every tc
 
    a = 1;
    z = a + tc/dt;
    for j = 1:nsteps
        this_dil = noise(j)*(dil_max - dil_min) + dil_min;
        
        ControlParadigm2(i).Outputs(2,a:z) = (this_dil*V/(1 - this_dil))/MFC_Scale;
        
        
        % increment
        a = z;
        z = min([a + tc/dt T/dt]);
    end
    
    m = mean(ControlParadigm2(i).Outputs(2,:));
    ControlParadigm2(i).Outputs(2,end) = m; % don't leave the ORN with some high or low value
    ControlParadigm2(i).Outputs(2,ControlParadigm2(i).Outputs(2,:) > 5) = 5; % clip for sanity
    ControlParadigm2(i).Outputs(2,ControlParadigm2(i).Outputs(2,:) < min_flow) = min_flow; % clip for sanity

    % add the main air
    ControlParadigm2(i).Outputs(3,:) = 1;
    
    % add the light
    ControlParadigm2(i).Outputs(1,:) = light_mean(i);
    
    
end

ControlParadigm = [ControlParadigm ControlParadigm2];

l = length(ControlParadigm)+1;

ControlParadigm(l).Name = 'end';
ControlParadigm(l).Outputs = zeros(3,1e4);

%% save
n = ('C:\Users\mahmut\Documents\MATLAB\LightOdourFlicker_Kontroller_paradigm.mat');
save(n,'ControlParadigm')

