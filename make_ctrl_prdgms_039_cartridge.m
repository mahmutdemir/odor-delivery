function ControlParadigm = make_ctrl_prdgms_039_cartridge(inittime,width,duration,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_039_cartridge(inittime,width,duration,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_039_cartridge(inittime,width,duration,save_name,saveit)
% This functions constructs a control paradigm.  Sampling rate for this paradigm is 10000 Hz. 

%%   Inputs
%
%   inittime:   the time prior the puff in seconds
%   width:      with of the puff in seconds
%   duration:   duration of the recording (sec)
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk



sr = 10000;  % digitization sampling rate Hz 



    % start generating paradigms
    voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
    voltages(6,:) = 1;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;
    
       
    nop = duration*sr;
    voltages = zeros(6,nop); % allocate space for the output matrix and set to zero
    voltages(6,:) = 1;  % write voltage for the main air MFC
    %1: led dim, 2: mfc500, 3: mfc:200, 4:
    %vmfc500, 5: vmfc200, 6:main
    voltages(4,inittime*sr+1:(inittime+width)*sr+1) = 1;  % write voltage for the valve turn on
    voltages(2,:) = 200/500*5;  % write voltage to set mfc to 200ml/min


    ControlParadigm(2).Name = 'pulse';
    ControlParadigm(2).Outputs = voltages;
                

    
    
    
        voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
        ControlParadigm(3).Name = 'end';
        ControlParadigm(3).Outputs = voltages;
        
    if saveit
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
    end