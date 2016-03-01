function ControlParadigm = make_ctrl_prdgms_017(boost_volt,boost_width,odor_rat,total_flow,keep_air_on,inittime,width,duration,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_017(boost_volt,boost_width,odor_rat,total_flow,keep_air_on,inittime,width,duration,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_017(boost_volt,boost_width,odor_rat,total_flow,keep_air_on,inittime,width,duration,save_name,saveit)
% This functions constructs a control paradigm for an odor puff with given
% pulse width (sec) and total duration of "duration" (sec). The total number of 
% pulses is one. The mixing ratio is defined as odor_rat:mix_rat which
% is odor and clean air partial ratios in the mix. odor mixing is initiated
% inittime seconds prior to the recording. Sampling rate for this paradigm is
% 10000 Hz. 
% Air is turned off at the end of each trial
% odor gate valve is added in order to switch between odor bottles
% this function makes control paradigms for low and high doses.
%%   Inputs
%
%   boost_volt:    vector containing voltage values
%   boost_width:   the time in seconds to boost odor MFC to 5 volts in the
%                  beginning of the dose in order to get it reach the set value before the
%                  puff.
%   odor_rat:      2X1 Cell [{low} {high}] containing vectors for odor air ratio in the final mixture
%   total_flow:       total flow rate of the mixture
%   keep_air_on: 1 for keeping depletion air on, 0 for turning air off at
%   the end of the pulse
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in sconds
%   duration:   total duration of the recording
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk

odor_rat_low = odor_rat{1};
odor_rat_high = odor_rat{2};
mix_rat_low = total_flow - odor_rat{1};
mix_rat_high = total_flow - odor_rat{2};


sr = 10000;  % digitization sampling rate Hz 
noparad = length(odor_rat_low)+length(odor_rat_high)+2;
conv_factor = total_flow/20;


    nop = (duration)*sr;
    voltages = zeros(5,10000); % allocate space for the output matrix and set to zero
    voltages(5,:) = 1;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;
    


    % start generating low doses
    for i = 1:length(odor_rat_low)
        
        
        odor_volt = conv_factor*odor_rat_low(i)/(odor_rat_low(i)+mix_rat_low(i));   % voltage for odor MFC. This will be written to AO0
        mix_volt = conv_factor*mix_rat_low(i)/(odor_rat_low(i)+mix_rat_low(i));   % voltage for air mix MFC. This will be written to AO1
        voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
        voltages(1,1:(inittime+width+1)*sr+1) = odor_volt;  % write voltage for the odor MFC
        
        if keep_air_on == 1
            if odor_volt ~=0
            voltages(1,(inittime+width+1)*sr+1:end) = 1;  % write voltage for the odor MFC. set to 200ml/min
            end
        end
        
        if odor_volt <1 && odor_volt~=0
            voltages(1,1:boost_width*sr+1) = boost_volt(i);  % write voltage for the odor MFC
        end
        
        voltages(2,1:(inittime+width+1)*sr+1) = mix_volt;  % write voltage for the mix air MFC

        voltages(5,:) = 1;  % write voltage for the main air MFC
        voltages(4,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the puff valve
        voltages(3,10000:end-10000) = 0;  % write voltage for the odor gate. it is closed for low dose
      
        ControlParadigm(i+1).Name = ['dil_low_' num2str(odor_rat_low(i)) '_' num2str(mix_rat_low(i))];
        ControlParadigm(i+1).Outputs = voltages;
        
                
    end
    
    % start generating high doses
    for i = 1:length(odor_rat_high)
        
        
        odor_volt = conv_factor*odor_rat_high(i)/(odor_rat_high(i)+mix_rat_high(i));   % voltage for odor MFC. This will be written to AO0
        mix_volt = conv_factor*mix_rat_high(i)/(odor_rat_high(i)+mix_rat_high(i));   % voltage for air mix MFC. This will be written to AO1
        voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
        voltages(1,1:(inittime+width+1)*sr+1) = odor_volt;  % write voltage for the odor MFC
        
        if keep_air_on == 1
            voltages(1,end-10000:end) = 1;  % write voltage for the odor MFC. set to 200ml/min
        end
        
        if odor_volt <1 && odor_volt~=0
            voltages(1,1:boost_width*sr+1) = boost_volt(i);  % write voltage for the odor MFC
        end
        
        voltages(2,1:(inittime+width+1)*sr+1) = mix_volt;  % write voltage for the mix air MFC

        voltages(5,:) = 1;  % write voltage for the main air MFC
        voltages(4,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the puff valve
        voltages(3,10000:end-10000) = 1;  % write voltage for the odor gate. it is open for low dose
      
        ControlParadigm(i+length(odor_rat_low)+1).Name = ['dil_high_' num2str(odor_rat_high(i)) '_' num2str(mix_rat_high(i))];
        ControlParadigm(i+length(odor_rat_low)+1).Outputs = voltages;
        
                
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