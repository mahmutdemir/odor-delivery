function ControlParadigm = make_ctrl_prdgms_020(odor_rat,mix_rat,odor_gate_on,corr_length,inittime,width,duration,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_020(odor_rat,mix_rat,valve_input,inittime,width,duration,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_020(odor_rat,mix_rat,valve_input,inittime,width,duration,save_name,saveit)
% This functions constructs a control paradigm for flickering odor stimuli
% pulse width (sec) and total duration of "duration" (sec). The total number of 
% pulses is one.  Sampling rate for this paradigm is 10000 Hz. 

%%   Inputs
%
%   odor_rat:      odor flow rate ratio in the total flow. Assuming 200
%   mlpmin total
%   mix_rat:      clean air flow rate ratio in the total flow. Assuming 200
%   mlpmin total
%   total_flow:       total flow rate of the mixture
%   corr_length: minimum length of the time that valve is open or closed in ms
%   odor_gate_on: to open gate for high dose line set to 1, otherwise 0 
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in sconds
%   duration:   total duration of the recording
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk



sr = 10000;  % digitization sampling rate Hz 
conv_factor = (odor_rat(1)+mix_rat(1))/20;


    nop = (duration)*sr;
    binary_series = make_pseudo_binary_series(corr_length,sr,width);
    
    voltages = zeros(5,10000); % allocate space for the output matrix and set to zero
    voltages(5,:) = 1;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;

    % start generating paradigms
    for i = 1:length(odor_rat)

        odor_volt = conv_factor*odor_rat(i)/(odor_rat(i)+mix_rat(i));   % voltage for odor MFC. This will be written to AO0
        mix_volt = conv_factor*mix_rat(i)/(odor_rat(i)+mix_rat(i));   % voltage for air mix MFC. This will be written to AO1
        voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
        voltages(1,1:(inittime+width+1)*sr+1) = odor_volt;  % write voltage for the odor MFC
        voltages(2,1:(inittime+width+1)*sr+1) = mix_volt;  % write voltage for the mix air MFC
        voltages(5,:) = 1;  % write voltage for the main air MFC
        voltages(4,inittime*sr+1:(inittime+width)*sr) = binary_series';  % write voltage for the puff valve
        voltages(3,10000:end-10000) = odor_gate_on;  % write voltage for the odor gate. it is closed for low dose
      
        ControlParadigm(i+1).Name = ['flick_' num2str(odor_rat(i)) '_' num2str(mix_rat(i))];
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