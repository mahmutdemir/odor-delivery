function make_ctrl_prdgms_005(odor_rat,mix_rat,inittime,width,duration,n_of_trials,name,save_name)
%  make_ctrl_prdgms_005(odor_rat,mix_rat,inittime,width,duration,n_of_trials,name,save_name)
% function make_ctrl_prdgms_005(odor_rat,mix_rat,inittime,width,duration,n_of_trials,name,save_name)
% This functions constructs a control paradigm for an odor puff with given
% pulse width (sec) and total duration of "duration" (sec). The total number of 
% pulses is one. The mixing ratio is defined as odor_rat:mix_rat which
% is odor and clean air partial ratios in the mix. odor mixing is initiated
% inittime seconds prior to the recording. Sampling rate for this paradigm is
% 10000 Hz. Number of intented trials is used to contruct the paradigm such
% that at the last trial the mixing is iniated in order to get it ready for
% the next dose.this function is good for ramp up. For ramp down air flow
% mixing adjustements should be corrected

sr = 10000;  % digitization sampling rate Hz 
noparad = length(odor_rat)+2;

    nop = (inittime+duration)*sr;
    odor_volt = 1*odor_rat(1)/(odor_rat(1)+mix_rat(1));   % voltage for odor MFC. This will be written to AO0
    mix_volt = 1*mix_rat(1)/(odor_rat(1)+mix_rat(1));   % voltage for air mix MFC. This will be written to AO1
    voltages = zeros(4,10000); % allocate space for the output matrix and set to zero
    voltages(1,:) = odor_volt;  % write voltage for the odor MFC
    if odor_volt<1
       voltages(1,1:500) = 2;  % write voltage for the odor MFC 
    end
    voltages(2,:) = mix_volt;  % write voltage for the mix air MFC
    if mix_volt<1
        if mix_volt~=0
            voltages(2,1:500) = 2;  % write voltage for the odor MFC
        end
    end
    voltages(4,:) = 1;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;
    
if n_of_trials==1
    
    for i = 2:noparad-1
        
        
        
    odor_volt = 1*odor_rat(i-1)/(odor_rat(i-1)+mix_rat(i-1));   % voltage for odor MFC. This will be written to AO0
    mix_volt = 1*mix_rat(i-1)/(odor_rat(i-1)+mix_rat(i-1));   % voltage for air mix MFC. This will be written to AO1
    voltages = zeros(4,nop); % allocate space for the output matrix and set to zero
    voltages(1,:) = odor_volt;  % write voltage for the odor MFC
    voltages(2,:) = mix_volt;  % write voltage for the mix air MFC
    voltages(4,:) = 1;  % write voltage for the main air MFC
    voltages(3,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the puff valve

    ControlParadigm(i).Name = ['dilution_' name{i-1}];
    ControlParadigm(i).Outputs = voltages;

    end
    
    voltages = zeros(4,10000); % allocate space for the output matrix and set to zero

    ControlParadigm(noparad).Name = 'end';
    ControlParadigm(noparad).Outputs = voltages;
    
    
elseif n_of_trials>1   
    
    for i = 2:noparad-1
        odor_volt = 1*odor_rat(i-1)/(odor_rat(i-1)+mix_rat(i-1));   % voltage for odor MFC. This will be written to AO0
        mix_volt = 1*mix_rat(i-1)/(odor_rat(i-1)+mix_rat(i-1));   % voltage for air mix MFC. This will be written to AO1
        voltages = zeros(4,nop); % allocate space for the output matrix and set to zero
        voltages(1,:) = odor_volt;  % write voltage for the odor MFC
        voltages(2,:) = mix_volt;  % write voltage for the mix air MFC
        voltages(4,:) = 1;  % write voltage for the main air MFC
        voltages(3,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the puff valve

        ControlParadigm(i*2-2).Name = ['dilution_' name{i-1}];
        ControlParadigm(i*2-2).Outputs = voltages;

        if i~=noparad-1
            if odor_volt<1
                voltages(1,(inittime+width+0.25)*sr+1:(inittime+width+1.25)*sr+1) = 1;  % write voltage for the odor MFC 
            end
                odor_volt = 1*odor_rat(i)/(odor_rat(i)+mix_rat(i));   % voltage for odor MFC. This will be written to AO0
            mix_volt = 1*mix_rat(i)/(odor_rat(i)+mix_rat(i));   % voltage for air mix MFC. This will be written to AO1
            voltages(1,(inittime+width+1.25)*sr+1:end) = odor_volt;  % write voltage for the odor MFC
            voltages(2,(inittime+width+0.5)*sr+1:end) = mix_volt;  % write voltage for the mix air MFC
            voltages(3,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the puff valve

            ControlParadigm(i*2-1).Name = ['dilution_' name{i-1} '_final'];
            ControlParadigm(i*2-1).Outputs = voltages;
        end

    end
    voltages = zeros(4,10000); % allocate space for the output matrix and set to zero

    ControlParadigm(2*noparad-3).Name = 'end';
    ControlParadigm(2*noparad-3).Outputs = voltages;

else error('The number of trilas cannot be zero')
end

    


    argnames = [{'initial time'},{'pulse width'},{'duration'}];
    arg_values = [{inittime},{width},{duration}];
    Parameters = struct('arg_names',argnames,'arg_values',arg_values);
    save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
    save ([save_name,'_parameters.mat'], 'Parameters')