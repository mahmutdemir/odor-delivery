function ControlParadigm = make_ctrl_prdgms_030(main,odor,clean,bckg)
%  ControlParadigm = make_ctrl_prdgms_030(main,odor,clean,bckg)
% function ControlParadigm = make_ctrl_prdgms_030(main,odor,clean,bckg)
% This functions constructs a control paradigm for MFC's so that they 
% will be set to a constant value . Sampling rate for this paradigm is
% 10000 Hz. keeps the clean air on for cleaning the line
% This function is constructed for the DAQ with 4 Anolog outputs.
% MFC's in use are as following:
%   Main Air: Aalborg 5l/min
%   Odor: Alicat 500ml/min
%   Clean: Alicat 500ml/min
%   Background: 200ml/min



sr = 10000;  % digitization sampling rate Hz 

       
odor_volt = odor*5/500;   % voltage for odor MFC. This will be written to AO0
mix_volt = clean*5/500;   % voltage for odor MFC. This will be written to AO0
voltages = zeros(6,10000); % allocate space for the output matrix and set to zero
voltages(1,:) = main/5000*5;  % write voltage for the main air MFC
voltages(2,:) = odor_volt;  % write voltage for the odor MFC
voltages(3,:) = mix_volt;  % write voltage for the mix air MFC
voltages(4,:) = bckg*5/200;   % write voltage for the background air MFC
ControlParadigm(1).Name = ['dil_' num2str(odor) '_' num2str(clean)];
ControlParadigm(1).Outputs = voltages;
        
                
  