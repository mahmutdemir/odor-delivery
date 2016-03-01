function ControlParadigm = make_ctrl_prdgms_033(odor_vol,mix_vol,bckg_vol,inittime,width,lag,interval,nofpulses,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_033(odor_vol,mix_vol,bckg_vol,inittime,width,lag,interval,nofpulses,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_033(odor_vol,mix_vol,bckg_vol,inittime,width,lag,interval,nofpulses,save_name,saveit)
% This functions constructs a control paradigm for two pulses seperated with 
% lag seconds and construct trials with lags upto noflags*lag. initial
% pulse is applied at the inittime (sec). Total number of pulse couples are
% nofpulse. The pulse couples are separated by interval seconds. The total
% number of paradigms to be constructed is: length(odor_vol)*length(bckg_vol)*length(width)*length(lag)

% Sampling rate for this paradigm is 10000 Hz. 

%%   Inputs
%
%   odor_vol:   (vector) volume (ml/min) of the air flow rate passing through the odor
%   bottle 
%   mix_vol:    (vector) volume of the clean air flow rate mixed with odorant 
%   bckg_vol:    (vector) volume of the air flow rate passing from background odor bottle 
%   inittime:   the time prior and after the puffs in seconds
%   width:      (vector) with of the puffs in sconds
%   lag:        (vector) time interval between pulses (sec)
%   interval:     time (sec) between pulse couples
%   nofpulses:  total number of pulse couples in a run
%   save_name:   paradigm name string for saving
%   save:       0= do not save 1= save the paradigm to the disk



sr = 10000;  % digitization sampling rate Hz 



    


    % now initiate the output matrix
    % 1: main, 2: Odor, 3: Clean, 4: Background, 5: odor_puff, 6: bckg_puff
    voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
    voltages(1,:) = 2;  % write voltage for the main air MFC

    ControlParadigm(1).Name = 'start';
    ControlParadigm(1).Outputs = voltages;

    % start generating paradigms
    parad_num = 1;
    cflag = zeros(length(odor_vol),length(width));
    pflag = zeros(length(bckg_vol),length(width));

    for wn = 1:length(width)
        for i = 1:length(lag)
            for cn = 1:length(odor_vol)
                for pn = 1:length(bckg_vol)
                    
                    if ~(odor_vol(cn)==0 && bckg_vol(pn)==0)
                        if odor_vol(cn)==0 && pflag(pn,wn) ==0
                            pflag(pn,wn) = 1;

                            nop = (2*inittime+(2*width(wn)+lag(i)+interval)*nofpulses-interval)*sr;
                            voltages = zeros(6,nop); % allocate space for the output matrix and set to zero
                            odor_volt = odor_vol(cn)/500*5;   % voltage for odor MFC (Alicat 500ml)
                            mix_volt = mix_vol(cn)/500*5;   % voltage for air mix MFC (Alicat 500ml)
                            bckg_volt = bckg_vol(pn)/200*5;   % voltage for air mix MFC (Alicat 200ml)
                            voltages(1,:) = 2;  % write voltage for the main air MFC
                            voltages(2,1:end-1) = odor_volt;  % write voltage for the odor MFC
                            voltages(3,1:end-1) = mix_volt;  % write voltage for the mix air MFC
                            voltages(4,1:end-1) = bckg_volt;  % write voltage for the mix air MFC


                            for j = 1: nofpulses
                                if odor_vol(cn)~=0
                                    voltages(5,round((inittime+(j-1)*(interval+lag(i)))*sr)+1:...
                                        round((inittime+width(wn)+(j-1)*(interval+lag(i)))*sr)+1) = 1;  % write voltage for the puff valve
                                end
                                if bckg_vol(pn)~=0
                                    voltages(6,round((inittime+(j-1)*(interval+lag(i))+width(wn)+lag(i))*sr)+1:...
                                        round((inittime+width(wn)+(j-1)*(interval+lag(i))+width(wn)+lag(i))*sr)+1) = 1;  % write voltage for the puff valve
                                end
                            end

                            parad_num = parad_num + 1; 
                            ControlParadigm(parad_num).Name = ['c_' num2str(odor_vol(cn)) '_' num2str(mix_vol(cn)) ...
                                '_p_' num2str(bckg_vol(pn)) '_l_' num2str(lag(i)*1000) '_w_' num2str(width(wn)*1000)];
                            ControlParadigm(parad_num).Outputs = voltages;
                        
                        elseif odor_vol(cn)==0 && pflag(pn,wn) == 1
                        elseif bckg_vol(pn)==0 && cflag(cn,wn) == 0
                            cflag(cn,wn) = 1;
                            nop = (2*inittime+(2*width(wn)+lag(i)+interval)*nofpulses-interval)*sr;
                            voltages = zeros(6,nop); % allocate space for the output matrix and set to zero
                            odor_volt = odor_vol(cn)/500*5;   % voltage for odor MFC (Alicat 500ml)
                            mix_volt = mix_vol(cn)/500*5;   % voltage for air mix MFC (Alicat 500ml)
                            bckg_volt = bckg_vol(pn)/200*5;   % voltage for air mix MFC (Alicat 200ml)
                            voltages(1,:) = 2;  % write voltage for the main air MFC
                            voltages(2,1:end-1) = odor_volt;  % write voltage for the odor MFC
                            voltages(3,1:end-1) = mix_volt;  % write voltage for the mix air MFC
                            voltages(4,1:end-1) = bckg_volt;  % write voltage for the mix air MFC


                            for j = 1: nofpulses
                                if odor_vol(cn)~=0
                                    voltages(5,round((inittime+(j-1)*(interval+lag(i)))*sr)+1:...
                                        round((inittime+width(wn)+(j-1)*(interval+lag(i)))*sr)+1) = 1;  % write voltage for the puff valve
                                end
                                if bckg_vol(pn)~=0
                                    voltages(6,round((inittime+(j-1)*(interval+lag(i))+width(wn)+lag(i))*sr)+1:...
                                        round((inittime+width(wn)+(j-1)*(interval+lag(i))+width(wn)+lag(i))*sr)+1) = 1;  % write voltage for the puff valve
                                end
                            end

                            parad_num = parad_num + 1; 
                            ControlParadigm(parad_num).Name = ['c_' num2str(odor_vol(cn)) '_' num2str(mix_vol(cn)) ...
                                '_p_' num2str(bckg_vol(pn)) '_l_' num2str(lag(i)*1000) '_w_' num2str(width(wn)*1000)];
                            ControlParadigm(parad_num).Outputs = voltages;
                        elseif bckg_vol(pn)==0 && cflag(cn,wn) == 1
                        else
                            nop = (2*inittime+(2*width(wn)+lag(i)+interval)*nofpulses-interval)*sr;
                            voltages = zeros(6,nop); % allocate space for the output matrix and set to zero
                            odor_volt = odor_vol(cn)/500*5;   % voltage for odor MFC (Alicat 500ml)
                            mix_volt = mix_vol(cn)/500*5;   % voltage for air mix MFC (Alicat 500ml)
                            bckg_volt = bckg_vol(pn)/200*5;   % voltage for air mix MFC (Alicat 200ml)
                            voltages(1,:) = 2;  % write voltage for the main air MFC
                            voltages(2,1:end-1) = odor_volt;  % write voltage for the odor MFC
                            voltages(3,1:end-1) = mix_volt;  % write voltage for the mix air MFC
                            voltages(4,1:end-1) = bckg_volt;  % write voltage for the mix air MFC


                            for j = 1: nofpulses
                                if odor_vol(cn)~=0
                                    voltages(5,round((inittime+(j-1)*(interval+lag(i)))*sr)+1:...
                                        round((inittime+width(wn)+(j-1)*(interval+lag(i)))*sr)+1) = 1;  % write voltage for the puff valve
                                end
                                if bckg_vol(pn)~=0
                                    voltages(6,round((inittime+(j-1)*(interval+lag(i))+width(wn)+lag(i))*sr)+1:...
                                        round((inittime+width(wn)+(j-1)*(interval+lag(i))+width(wn)+lag(i))*sr)+1) = 1;  % write voltage for the puff valve
                                end
                            end

                            parad_num = parad_num + 1; 
                            ControlParadigm(parad_num).Name = ['c_' num2str(odor_vol(cn)) '_' num2str(mix_vol(cn)) ...
                                '_p_' num2str(bckg_vol(pn)) '_l_' num2str(lag(i)*1000) '_w_' num2str(width(wn)*1000)];
                            ControlParadigm(parad_num).Outputs = voltages;
                        end
                        
                    end

                    
                end
                
            end

        end
    end
         
    
        voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
        ControlParadigm(parad_num + 1).Name = 'end';
        ControlParadigm(parad_num + 1).Outputs = voltages;
        
    if saveit
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
    end