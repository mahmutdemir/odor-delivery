function ControlParadigm = make_ctrl_prdgms_032(odor_vol,mix_vol,bckg_vol,inittime,width,lag,noflags,t_step,nofpulses,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_032(odor_vol,mix_vol,bckg_vol,inittime,width,lag,noflags,t_sep,nofpulses,save_name,saveit)
% function ControlParadigm =make_ctrl_prdgms_032(odor_vol,mix_vol,bckg_vol,inittime,width,lag,noflags,t_sep,nofpulses,save_name,saveit)
% This functions constructs a control paradigm for two pulses seperated with 
% lag seconds and construct trials with lags upto noflags*lag. initial
% pulse is applied at the inittime (sec). Total number of pulse couples are
% nofpulse. The pulse couples are separated by t_sep seconds. For each las
% time two paradigms are constructed such that; 1. paradigm: tall-short (ts)
% and 2. paradigm is short-tall (st) pulse couples

% Sampling rate for this paradigm is 10000 Hz. 

%%   Inputs
%
%   odor_vol:   volume (ml/min) of the air flow rate passing through the odor
%   bottle
%   mix_vol:    volume of the clean air flow rate mixed with odorant 
%   bckg_vol:    volume of the air flow rate passing from background odor bottle 
%   inittime:   the time prior and after the puffs in seconds
%   width:      with of the puffs in sconds
%   lag:        time interval between pulses (sec)
%   noflags:    number of lags to be generated
%   t_step:     time (sec) between pulse couples
%   nofpulses:  total number of pulse couples in a run
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk



sr = 10000;  % digitization sampling rate Hz 



    


    % now initiate the output matrix
    % 1: main, 2: Odor, 3: Clean, 4: Background, 5: odor_puff, 6: bckg_puff
    voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
    voltages(1,:) = 2;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;

    % start generating paradigms

    
    for i = 1:noflags
        
        nop = (2*inittime+(2*width+lag*i+t_step)*nofpulses-t_step)*sr;
        voltages = zeros(6,nop); % allocate space for the output matrix and set to zero
        odor_volt = odor_vol/500*5;   % voltage for odor MFC (Alicat 500ml)
        mix_volt = mix_vol/500*5;   % voltage for air mix MFC (Alicat 500ml)
        bckg_volt = bckg_vol/200*5;   % voltage for air mix MFC (Alicat 200ml)
        voltages(1,:) = 2;  % write voltage for the main air MFC
        voltages(2,1:end-1) = odor_volt;  % write voltage for the odor MFC
        voltages(3,1:end-1) = mix_volt;  % write voltage for the mix air MFC
        voltages(4,1:end-1) = bckg_volt;  % write voltage for the mix air MFC
        
       
        for j = 1: nofpulses
            voltages(5,round((inittime+(j-1)*(t_step+lag*i))*sr)+1:round((inittime+width+(j-1)*(t_step+lag*i))*sr)+1) = 1;  % write voltage for the puff valve
            voltages(6,round((inittime+(j-1)*(t_step+lag*i)+width+lag*i)*sr)+1:round((inittime+width+(j-1)*(t_step+lag*i)+width+lag*i)*sr)+1) = 1;  % write voltage for the puff valve
        end
      
        ControlParadigm(2*i).Name = ['ts_' num2str((i)*lag*1000) '_ms_lag'];
        ControlParadigm(2*i).Outputs = voltages;
        
        voltages(6,:) = 0;  % write voltage for the puff valve
        voltages(5,:) = 0;  % write voltage for the puff valve
        for j = 1: nofpulses
           voltages(6,round((inittime+(j-1)*(t_step+lag*i))*sr)+1:round((inittime+width+(j-1)*(t_step+lag*i))*sr)+1) = 1;  % write voltage for the puff valve
           voltages(5,round((inittime+(j-1)*(t_step+lag*i)+width+lag*i)*sr)+1:round((inittime+width+(j-1)*(t_step+lag*i)+width+lag*i)*sr)+1) = 1;  % write voltage for the puff valve
        end
      
        ControlParadigm(2*i+1).Name = ['st_' num2str((i)*lag*1000) '_ms_lag'];
        ControlParadigm(2*i+1).Outputs = voltages;
        
    end
         
    
        voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
        ControlParadigm(2*i+2).Name = 'end';
        ControlParadigm(2*i+2).Outputs = voltages;
        
    if saveit
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
    end