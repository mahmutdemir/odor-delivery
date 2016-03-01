function ControlParadigm = make_ctrl_prdgms_027(odor_rat,total,inittime,width,duration,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_027(odor_rat,total_flow,inittime,width,duration,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_027(odor_rat,total_flow,inittime,width,duration,save_name,saveit)
% This functions constructs a control paradigm for an odor puff with given
% pulse width (sec) and total duration of "duration" (sec). The total number of 
% pulses is one. The mixing ratio is defined as odor_rat:mix_rat which
% is odor and clean air partial ratios in the mix. odor mixing is initiated
% inittime seconds prior to the recording. Sampling rate for this paradigm is
% 10000 Hz. keeps the clean air on for cleaning the line
% This function is constructed for clean air mfc : 500 ml alicat and
% odorant air mfs is 200 ml alicat
%%   Inputs
%
%   odor_rat:         vector contining values for odor air ratio in the final mixture
%   total:            total odor percentage value for the mixture 
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in sconds
%   duration:   total duration of the recording
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk


sr = 10000;  % digitization sampling rate Hz 
noparad = length(odor_rat)+2;

odor_air_mfc_cf = 200/5; % ml/volt
clean_air_mfc_cf = 500/5; % ml/volt
mix_volt_clean = 300/clean_air_mfc_cf;   % voltage for odor MFC. This will be written to AO0


    nop = (duration)*sr;
    voltages = zeros(5,10000); % allocate space for the output matrix and set to zero
    voltages(5,:) = 1;  % write voltage for the main air MFC
    voltages(2,:) = mix_volt_clean;  % write voltage for the mix air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;
    


    % start generating low doses
    for i = 1:length(odor_rat)
        
        
        odor_volt = 200/total*odor_rat(i)/odor_air_mfc_cf;   % voltage for odor MFC. This will be written to AO0
        mix_volt = (200-200/total*odor_rat(i))/clean_air_mfc_cf;   % voltage for odor MFC. This will be written to AO0
        voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
        voltages(1,1:(inittime+width+1)*sr+1) = odor_volt;  % write voltage for the odor MFC
        voltages(2,1:(inittime+width+1)*sr+1) = mix_volt;  % write voltage for the mix air MFC
        voltages(2,(inittime+width+1)*sr+1:end) = mix_volt_clean;  % write voltage for the mix air MFC
        voltages(5,:) = 1;  % write voltage for the main air MFC
        voltages(4,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the puff valve
        voltages(3,10000:end-10000) = 0;  % write voltage for the odor gate. This will not be used for this case
      
        ControlParadigm(i+1).Name = ['dil_' num2str(odor_rat(i)) '_' num2str(total-odor_rat(i))];
        ControlParadigm(i+1).Outputs = voltages;
        
                
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