function ControlParadigm = make_ctrl_prdgms_025(odor_rat,mix_rat,inittime,width,lag,noflags,duration,gate,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_025(odor_rat,mix_rat,inittime,width,lag,noflags,duration,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_025(odor_rat,mix_rat,inittime,width,lag,noflags,duration,save_name,saveit)
% This functions constructs a control paradigm for two identical pulses seperated with 
% lag milliseconds and construct trials with lags upto noflags*lag. initial
% pulse is applied at the inittime (sec). Total length of the record is
% duration (sec) and with is in seconds. Sampling rate for this paradigm is 10000 Hz. 

%%   Inputs
%
%   odor_rat:     ratio of the air flow rate passing through the odor
%   bottle
%   mix_rat:    ratio of the clean air flow rate mixed with odorant 
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in sconds
%   lag:        time interval between pulses
%   noflags:    number of lags to be generated
%   duration:   total duration of the recording
%   gate:       0: closed, low dose line, 1: open, high dose line
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk



sr = 10000;  % digitization sampling rate Hz 


    nop = (duration)*sr;
    voltages = zeros(5,10000); % allocate space for the output matrix and set to zero
    voltages(5,:) = 1;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;

    % start generating paradigms
        voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
        odor_volt = odor_rat/20;   % voltage for odor MFC. This will be written to AO0
        mix_volt = mix_rat/20;   % voltage for air mix MFC. This will be written to AO1
        voltages(1,1:end-1) = odor_volt;  % write voltage for the odor MFC
        voltages(2,1:end-1) = mix_volt;  % write voltage for the mix air MFC
        voltages(1,1:.5*sr) = 2;  % write voltage for the odor MFC
        voltages(2,1:.5*sr) = 2;  % write voltage for the mix air MFC
        voltages(3,.5*sr+1:(duration-.5)*sr) = gate;  % write voltage for the odor gate. it is closed for low dose
        voltages(5,:) = 1;  % write voltage for the main air MFC
    
    for i = 1:noflags

        voltages(4,:) = 0; % allocate space for the output matrix and set to zero
        voltages(4,inittime*sr+1:(inittime+width)*sr) = 1;  % write voltage for the puff valve
        voltages(4,(inittime+width)*sr+(i-1)*lag*10:(inittime+2*width)*sr-1+(i-1)*lag*10) = 1;  % write voltage for the puff valve
      
        ControlParadigm(i+1).Name = [num2str((i-1)*lag) '_ms_lag'];
        ControlParadigm(i+1).Outputs = voltages;
    end
         
    
        voltages = zeros(5,10000); % allocate space for the output matrix and set to zero
        ControlParadigm(i+2).Name = 'end';
        ControlParadigm(i+2).Outputs = voltages;
        
    if saveit
            argnames = [{'initial time'},{'pulse width'},{'duration'}];
            arg_values = [{inittime},{width},{duration}];
            Parameters = struct('arg_names',argnames,'arg_values',arg_values);
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
            save ([save_name,'_parameters.mat'], 'Parameters')
    end