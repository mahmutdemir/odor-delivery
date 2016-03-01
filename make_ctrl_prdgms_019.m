function ControlParadigm = make_ctrl_prdgms_019(dim_volt,inittime,width,duration,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_018(dim_volt,inittime,width,duration,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_018(dim_volt,inittime,width,duration,save_name,saveit)
% This functions constructs a control paradigm for dimmebale LED control
% pulse width (sec) and total duration of "duration" (sec). The total number of 
% pulses is one.  Sampling rate for this paradigm is 10000 Hz. 

%%   Inputs
%
%   dim_volt:    vector containing voltage values for dimming
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in sconds
%   duration:   total duration of the recording
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk
%
%   the following arguments are internal and set to zero.
%   boost_width:   the time in seconds to boost odor MFC to 5 volts in the
%                  beginning of the dose in order to get it reach the set value before the
%                  puff.
%   odor_rat:      2X1 Cell [{low} {high}] containing vectors for odor air ratio in the final mixture
%   total_flow:       total flow rate of the mixture
%   keep_air_on: 1 for keeping depletion air on, 0 for turning air off at
%   the end of the pulse


sr = 10000;  % digitization sampling rate Hz 


    nop = (duration)*sr;

    % start generating paradigms
    for i = 1:length(dim_volt)

        voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
        voltages(1,1:end-5000) = dim_volt(i);  % write voltage for the odor MFC
        voltages(4,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the light pulse on
        voltages(5,:) = 1;  % write voltage for the main air MFC
      
        ControlParadigm(i).Name = ['dim_V_' num2str(dim_volt(i))];
        ControlParadigm(i).Outputs = voltages;
        
                
    end
    
    
    
        voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
        ControlParadigm(i+1).Name = 'end';
        ControlParadigm(i+1).Outputs = voltages;
        
    if saveit
            argnames = [{'initial time'},{'pulse width'},{'duration'}];
            arg_values = [{inittime},{width},{duration}];
            Parameters = struct('arg_names',argnames,'arg_values',arg_values);
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
            save ([save_name,'_parameters.mat'], 'Parameters')
    end