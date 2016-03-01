function ControlParadigm = make_ctrl_prdgms_042_WA_pulses(vol1,vol2,inittime,width,interval,n_of_pulses,save_name,saveit,main_air)
%  ControlParadigm = make_ctrl_prdgms_042_WA_pulses(vol1,vol2,inittime,width,interval,n_of_pulses,save_name,saveit,main_air)
% function ControlParadigm = make_ctrl_prdgms_041_WA_pulses(vol1,vol2,inittime,width,interval,n_of_pulses,save_name,saveit,main_air)
% This functions constructs a control paradigm.  Sampling rate for this paradigm is 1000 Hz.
% Devices that kontrol are:
%   Aalborg 1L/min MFC:
%   Main gate: 2 way valve to open turn on/off main air
%   3-Way Solenoid Valve: for odor modulation
%   Output Channles List:
%       - mfc1
%       - mfc2
%       - main_v
%       - seed_v

%%   Inputs
%
%   vol:    ml/min volume of the air
%   inittime:   the time prior the puff in seconds (odor is on for that
%   time) prior to the pulses
%   width:      with of the puff in seconds
%   interval:   interval (sec) between pulses
%   n_of_pulses:   number of pulses to apply
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk
%   main_air: -1: off, 0-1: start at the point main_air equals to the
%   normnalized time. tmax = 1;

odorvol = vol1;
cleanvol = vol2;

sr = 1000;  % digitization sampling rate Hz 


switch nargin
    case 9
    case 8
        main_air = 0;   % start in the beginning
    case 7
        main_air = 0;   % start in the beginning
        saveit = 1;     % save the generated file
    case 6
        error('Not enough input arguments.')
end
    


% start generating paradigms
voltages = zeros(4,1*sr); % allocate space for the output matrix and set to zero
voltages(3,:) = 1;  % turn on the main air
ControlParadigm(1).Name = 'start';
ControlParadigm(1).Outputs = voltages;


if length(main_air) == length(odorvol)
elseif length(main_air)==1
    main_air = ones(1,length(odorvol))*main_air;
else
    error('main_air has to a single value or same length as dilutions')
end
    
parad_num = 2; 

    for wn = 1:length(width)
        
        nop = (n_of_pulses*(interval+width(wn))+2*inittime)*sr;
        % main air only
        voltages = zeros(4,nop); % allocate space for the output matrix and set to zero
        voltages(3,:) = 1;  % turn on the main air
        ControlParadigm(parad_num).Name = ['main_only_w_', num2str(width(wn))];
        ControlParadigm(parad_num).Outputs = voltages;
        % no air
        voltages = zeros(4,nop); % allocate space for the output matrix and set to zero
        ControlParadigm(parad_num+1).Name = ['no_air_w_', num2str(width(wn))];
        ControlParadigm(parad_num+1).Outputs = voltages;
        parad_num = parad_num + 2;
        for cn = 1:length(odorvol)
            voltages = zeros(4,nop); % allocate space for the output matrix and set to zero
            voltages(1,1:end-1) = cleanvol(cn)/1000*5;  % write voltage for the MFC
            voltages(2,1:end-1) = odorvol(cn)/1000*5;  % write voltage for the MFC
            if main_air(cn) == -1
                voltages(3,:) = 0;  % turn off the main air
            elseif (main_air(cn)>=0)&&(main_air(cn)<=1)
                voltages(3,round((n_of_pulses*(interval+width(wn))+2*inittime)*main_air(cn)*sr)+1:end) = 1;
            else
                error('main_air has to be either -1 or between 0 and 1.')
            end
            voltages(4,1:inittime*sr+1) = 1; % turn the odor initially
            for i = 1:n_of_pulses
                voltages(4,((i*interval+(i-1)*width(wn))+inittime)*sr+1:(i*(interval+width(wn))+inittime)*sr+1) = 1;
            end
            ControlParadigm(parad_num).Name = ['o/c_', num2str(odorvol(cn)), '/', num2str(cleanvol(cn)), '_w_', num2str(width(wn)), '_int_' num2str(interval)];
            ControlParadigm(parad_num).Outputs = voltages;
            
            parad_num = parad_num + 1; 

       end
    end
   
    
    
    


    
    
    % end generating paradigms
    voltages = zeros(4,1*sr); % allocate space for the output matrix and set to zero
    ControlParadigm(parad_num).Name = 'end';
    ControlParadigm(parad_num).Outputs = voltages;
    
       
    
    if saveit
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
    end