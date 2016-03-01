function ControlParadigm = make_ctrl_prdgms_040_smokepulses(vol,width,interval,n_of_pulses,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_040_smokepulses(inittime,width,interval,n_of_pulses,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_040_smokepulses(inittime,width,interval,n_of_pulses,save_name,saveit)
% This functions constructs a control paradigm.  Sampling rate for this paradigm is 1000 Hz. 

%%   Inputs
%
%   vol:    ml/min volume of the air
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in seconds
%   interval:   interval (sec) between pulses
%   n_of_pulses:   number of pulses to apply
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk



sr = 1000;  % digitization sampling rate Hz 


    % start generating paradigms
    voltages = zeros(3,1*sr); % allocate space for the output matrix and set to zero
    voltages(1,1:end-1) = vol/1000*5;  % write voltage for the MFC
    ControlParadigm(1).Name = ['wdth_' num2str(width) '_s_intvl_' num2str(interval) '_s'];
    ControlParadigm(1).Outputs = voltages;




    % start generating paradigms
    voltages = zeros(2,n_of_pulses*(interval+width)*sr+1); % allocate space for the output matrix and set to zero
    voltages(1,1:end-1) = vol/1000*5;  % write voltage for the MFC
    for i = 1:n_of_pulses
        voltages(2,(i*interval+(i-1)*width)*sr+1:i*(interval+width)*sr+1) = 1;
    end
    ControlParadigm(1).Name = ['wdth_' num2str(width) '_s_intvl_' num2str(interval) '_s'];
    ControlParadigm(1).Outputs = voltages;
    
       
    
    if saveit
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
    end