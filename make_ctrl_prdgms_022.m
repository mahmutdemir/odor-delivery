function ControlParadigm = make_ctrl_prdgms_022(maxfreq,maxtime,corr_length,gate_corr,inittime,width,duration,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_022(maxfreq,maxtime,corr_length,gate_corr,inittime,width,duration,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_022(maxfreq,maxtime,corr_length,gate_corr,inittime,width,duration,save_name,saveit)
% This functions constructs a control paradigm for natural flickering odor stimuli
% pulse width (sec) and total duration of "duration" (sec). The total number of 
% pulses is one.  Sampling rate for this paradigm is 10000 Hz. 

%%   Inputs
%
%   maxfreq:     max oscillation frequency for the MFC's
%   maxtime:    max oscillation time for the a given frequency
%   corr_length: minimum length of the time that valve is open or closed in ms
%   gate_corr: the factor for odor gate valve correlation 
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in sconds
%   duration:   total duration of the recording
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk



sr = 10000;  % digitization sampling rate Hz 
max_rate = 30;  % the max flow rate to be 300 ml/min


    nop = (duration)*sr;
    binary_series = make_pseudo_binary_series(corr_length,sr,width);
    binary_series_odor_gate = make_pseudo_binary_series(gate_corr*corr_length,sr,width);
    odor_rate = max_rate*make_oscillation_series(maxfreq,maxtime,sr,duration);
    mix_rate = -odor_rate;
    
    
    voltages = zeros(5,10000); % allocate space for the output matrix and set to zero
    voltages(5,:) = 1;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;

    % start generating paradigms

        odor_volt = odor_rate/20;   % voltage for odor MFC. This will be written to AO0
        mix_volt = mix_rate/20;   % voltage for air mix MFC. This will be written to AO1
        voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
        voltages(1,:) = odor_volt;  % write voltage for the odor MFC
        voltages(2,:) = mix_volt;  % write voltage for the mix air MFC
        voltages(5,:) = 1;  % write voltage for the main air MFC
        voltages(4,inittime*sr+1:(inittime+width)*sr) = binary_series';  % write voltage for the puff valve
        voltages(3,inittime*sr+1:(inittime+width)*sr) = binary_series_odor_gate';  % write voltage for the odor gate. it is closed for low dose
      
        ControlParadigm(2).Name = 'natural_flick';
        ControlParadigm(2).Outputs = voltages;
        
                       
    
    
    
        voltages = zeros(5,10000); % allocate space for the output matrix and set to zero
        ControlParadigm(3).Name = 'end';
        ControlParadigm(3).Outputs = voltages;
        
    if saveit
            argnames = [{'initial time'},{'pulse width'},{'duration'}];
            arg_values = [{inittime},{width},{duration}];
            Parameters = struct('arg_names',argnames,'arg_values',arg_values);
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
            save ([save_name,'_parameters.mat'], 'Parameters')
    end