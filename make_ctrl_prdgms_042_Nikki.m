function ControlParadigm = make_ctrl_prdgms_042_Nikki(odor_vol,mix_vol,inittime,width,save_name,saveit)
%  make_ctrl_prdgms_042_Nikki(odor_vol,mix_vol,inittime,width,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_042_Nikki(odor_vol,mix_vol,inittime,width,save_name,saveit)
% This functions constructs a control paradigm for a single pulse. initial
% pulse is applied at the inittime (sec).
%
% Sampling rate for this paradigm is 10000 Hz. 
%
%%   Outputs
%
%   LED
%   MFC2 (1000 ml/min)
%   MFC3 (500ml/min)
%   MFC500
%   puff valve
%   main air
%
%%   Inputs
%   odor_vol:   (vector) volume (ml/min) of the air flow rate passing through the odor
%   bottle 
%   mix_vol:    (vector) volume of the clean air flow rate mixed with odorant 
%   inittime:   the time prior and after the puffs in seconds
%   width:      (vector) with of the puffs in sconds
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk



sr = 10000;  % digitization sampling rate Hz 



    


    % now initiate the output matrix
    % 1: LED, 2: Clean, 3: Odor, 4: Flicker, 5: odor_puff, 6: main_air
    voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
    voltages(6,:) = 1;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;

    % start generating paradigms
    parad_num = 1;

    for wn = 1:length(width)
        for cn = 1:length(odor_vol)
            if width(wn)<=5
                nop = (2*inittime+width(wn))*sr;
            else
                nop = (2*inittime+width(wn)+10)*sr;
            end
            voltages = zeros(6,nop); % allocate space for the output matrix and set to zero
            odor_volt = odor_vol(cn)/500*5;   % voltage for odor MFC (Alicat 500ml)
            mix_volt = mix_vol(cn)/1000*5;   % voltage for air mix MFC (Alicat 1000ml)
            voltages(6,:) = 1;  % write voltage for the main air MFC
            voltages(3,1:end-1) = odor_volt;  % write voltage for the odor MFC
            voltages(3,end) =200/500*5 ;  % write voltage for the odor MFC
            voltages(2,1:end-1) = mix_volt;  % write voltage for the mix air MFC
            voltages(2,end) = 200/1000*5;  % write voltage for the mix air MFC, leave mix on for clening
            voltages(5,inittime*sr+1:(inittime+width(wn))*sr+1) = 1;  % write voltage for the valve turn on

            parad_num = parad_num + 1; 
            ControlParadigm(parad_num).Name = ['w_', num2str(width(wn)), '_o/c_' num2str(odor_vol(cn)) '/' num2str(mix_vol(cn))];
            ControlParadigm(parad_num).Outputs = voltages;

       end
    end
    % load gaussian noise for flickering
%         ga_co_par = load('LightOdourFlicker_Kontroller_paradigm.mat');
        ogfkp = load('C:\srinivas\Optimised_Gaussian_Flicker_Kontroller_Paradigm.mat');
        for cpe = 1:numel(ogfkp.ControlParadigm)-1
            voltages = zeros(6,length(ogfkp.ControlParadigm(cpe).Outputs(1,:))); % allocate space for the output matrix and set to zero
            ControlParadigm(parad_num + 1).Name = ogfkp.ControlParadigm(cpe).Name;
            voltages(4,:) = ogfkp.ControlParadigm(cpe).Outputs(2,:);
            voltages(6,:) = 1;  % write voltage for the main air MFC
            ControlParadigm(parad_num + 1).Outputs = voltages;
            parad_num = parad_num + 1;
        end
    
        voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
        ControlParadigm(parad_num + 1).Name = 'end';
        ControlParadigm(parad_num + 1).Outputs = voltages;
        
        
    if saveit
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
    end