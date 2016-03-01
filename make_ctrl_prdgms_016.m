function ControlParadigm = make_ctrl_prdgms_016(boost_volt,boost_width,odor_rat,mix_rat,inittime,width,duration,odor_gate_on,save_name,saveit)
%  function ControlParadigm = make_ctrl_prdgms_015(boost_volt,boost_width,odor_rat,mix_rat,inittime,width,duration,n_of_trials,predilfac,name,save_name,save)
% function ControlParadigm = make_ctrl_prdgms_015(boost_volt,boost_width,odor_rat,mix_rat,inittime,width,duration,n_of_trials,predilfac,name,save_name,save)
% This functions constructs a control paradigm for an odor puff with given
% pulse width (sec) and total duration of "duration" (sec). The total number of 
% pulses is one. The mixing ratio is defined as odor_rat:mix_rat which
% is odor and clean air partial ratios in the mix. odor mixing is initiated
% inittime seconds prior to the recording. Sampling rate for this paradigm is
% 10000 Hz. 
% Air is not turned off at the end of each trial
% odor gate valve is added in order to switch between odor bottles
%%   Inputs
%
%   boost_volt:    vector containing voltage values
%   boost_width:   the time in seconds to boost odor MFC to 5 volts in the
%                  beginning of the dose in order to get it reach the set value before the
%                  puff.
%   odor_rat:      Vector containing odor air ratio in the final mixture
%   mix_rat:       Vector containing clean air ratio in the final mixture
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in sconds
%   duration:   total duration of the recording
%   odor_gate_on: 0= gate is closed 1= gate if open 
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk

sr = 10000;  % digitization sampling rate Hz 
noparad = length(odor_rat)+2;
conv_factor = (odor_rat(1)+mix_rat(1))/20;


    nop = (duration)*sr;
    voltages = zeros(5,10000); % allocate space for the output matrix and set to zero
    voltages(5,:) = 1;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;
    


    
    for i = 2:noparad-1
        
        
        odor_volt = conv_factor*odor_rat(i-1)/(odor_rat(i-1)+mix_rat(i-1));   % voltage for odor MFC. This will be written to AO0
        mix_volt = conv_factor*mix_rat(i-1)/(odor_rat(i-1)+mix_rat(i-1));   % voltage for air mix MFC. This will be written to AO1
        voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
        voltages(1,1:end) = odor_volt;  % write voltage for the odor MFC
        voltages(1,(inittime+width+1)*sr+1:end) = 1;  % write voltage for the odor MFC
        voltages(2,:) = 1;  % write voltage for the mix air MFC
        
%         if odor_volt < 8/20
%             voltages(1,end) = 0;  % write voltage for the odor MFC
%         end
        
        if odor_volt <1 && odor_volt~=0
            voltages(1,1:boost_width*sr+1) = boost_volt(i-1);  % write voltage for the odor MFC
        end
        
        voltages(2,1:(inittime+width+1)*sr+1) = mix_volt;  % write voltage for the mix air MFC

        voltages(5,:) = 1;  % write voltage for the main air MFC
        voltages(4,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the puff valve
        voltages(3,10000:end-10000) = odor_gate_on;  % write voltage for the odor gate
      
        ControlParadigm(i).Name = ['dil_' num2str(odor_rat(i-1)) '_' num2str(mix_rat(i-1))];
        ControlParadigm(i).Outputs = voltages;
        
                
    end
    
    voltages = zeros(5,10000); % allocate space for the output matrix and set to zero

    ControlParadigm(noparad).Name = 'end';
    ControlParadigm(noparad).Outputs = voltages;


    if saveit
            argnames = [{'initial time'},{'pulse width'},{'duration'}];
            arg_values = [{inittime},{width},{duration}];
            Parameters = struct('arg_names',argnames,'arg_values',arg_values);
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
            save ([save_name,'_parameters.mat'], 'Parameters')
    end