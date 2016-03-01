function ControlParadigm = make_ctrl_prdgms_031(odor_rat,bckg,inittime,width,duration,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_031(odor_rat,total_flow,inittime,width,duration,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_031(odor_rat,total_flow,inittime,width,duration,save_name,saveit)
% This functions constructs a control paradigm for an odor puff with given
% pulse width (sec) and total duration of "duration" (sec). The total number of 
% pulses is one. The mixing ratio is defined as odor_rat:mix_rat which
% is odor and clean air partial ratios in the mix. odor mixing is initiated
% inittime seconds prior to the recording. Sampling rate for this paradigm is
% 10000 Hz. keeps the clean air on for cleaning the line. The odor line is
% kept on at 50 ml/min to deplete the odor to improve dose repeatibility.
% This function is constructed for the DAQ with 4 Anolog outputs.
% MFC's in use are as following:
%   Main Air: Aalborg 5l/min
%   Odor: Alicat 500ml/min
%   Clean: Alicat 500ml/min
%   Background: Alicat 200ml/min
%%   Inputs
%
%   odor_rat:         vector contining values for odor air volume in the
%   final mixture (total is 200 ml/min)
%   bckg:       background flow (ml/min)
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in sconds
%   duration:   total duration of the recording
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk


sr = 10000;  % digitization sampling rate Hz 
noparad = length(odor_rat)+2;

bckg_air_mfc_cf = 200/5; % ml/volt
odor_air_mfc_cf = 500/5; % ml/volt
mix_volt_clean = 300/odor_air_mfc_cf;   % voltage for odor MFC. This will be written to AO0
odor_deplete_volt = 50/odor_air_mfc_cf;   % voltage for odor MFC. This will be written to AO0


    nop = (duration)*sr;
    % now initiate the output matrix
    % 1: main, 2: Odor, 3: Clean, 4: Background, 5: odor_puff, 6: bckg_puff
    voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
    voltages(1,:) = 2;  % write voltage for the main air MFC
    voltages(2,:) = odor_deplete_volt;  % write voltage for the mix air MFC
    voltages(3,:) = mix_volt_clean;  % write voltage for the mix air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;
    


    % start generating doses
    for i = 1:length(odor_rat)
        
        
        odor_volt = odor_rat(i)/odor_air_mfc_cf;   % voltage for odor MFC. This will be written to AO0
        mix_volt = (200-odor_rat(i))/odor_air_mfc_cf;   % voltage for odor MFC. This will be written to AO0
        voltages = zeros(6,nop); % allocate space for the output matrix and set to zero
        voltages(1,:) = 2;  % write voltage for the main air MFC
        voltages(2,1:(inittime+width+1)*sr+1) = odor_volt;  % write voltage for the odor MFC
        voltages(2,(inittime+width+1)*sr+1:end) = odor_deplete_volt;  % write voltage for the mix air MFC
        voltages(3,1:(inittime+width+1)*sr+1) = mix_volt;  % write voltage for the mix air MFC
        voltages(3,(inittime+width+1)*sr+1:end) = mix_volt_clean;  % write voltage for the mix air MFC
        voltages(4,:) = bckg/bckg_air_mfc_cf;   % write voltage for the background air MFC
        voltages(5,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the puff valve
        if bckg>0
            voltages(6,10000:end-10000) = 1;  % write voltage for the background valve
        end
        ControlParadigm(i+1).Name = ['dil_' num2str(odor_rat(i)) '_' num2str(200-odor_rat(i))];
        ControlParadigm(i+1).Outputs = voltages;
        
                
    end
    
    voltages = zeros(6,10000); % allocate space for the output matrix and set to zero

    ControlParadigm(noparad).Name = 'end';
    ControlParadigm(noparad).Outputs = voltages;


    if saveit
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
    end