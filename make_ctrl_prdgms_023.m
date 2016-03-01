function ControlParadigm = make_ctrl_prdgms_023(odor_or_mix,max_flow,duration,odor_gate_on,save_name,saveit)
%  ControlParadigm = make_ctrl_prdgms_023(main_or_mix,max_flow,duration,odor_gate_on,save_name,saveit)
% function ControlParadigm = make_ctrl_prdgms_023(main_or_mix,max_flow,duration,odor_gate_on,save_name,saveit)
% This functions constructs a control paradigm where either the odor or
% mixture flow is increased from zero to max flow in given time
% Sampling rate for this paradigm is 10000 Hz. 
% Air is turned off at the end of each trial
% odor gate valve is added in order to switch between odor bottles
% this function makes control paradigms for low and high doses.
%%   Inputs
%
%   odor_or_mix:   1: odor, 2: mix
%   max_flow:      maximum flow rate to be reached
%   duration:      the duration for the increase process 
%   odor_gate_on:      0: closed, 1: open
%   save_name:   paradigm name string for saving
%   saveit:       0= do not save 1= save the paradigm to the disk


if max_flow > 100
    disp('Cannot apply flow rates larger than 1L/min. It will be set to 1L/min')
    max_flow = 100;
end

sr = 10000; % sampling rate
nop = duration*sr;  % total number of points for the paradigm vectors
flow_volt = max_flow/20;    % voltage to be applied to MFC

    voltages = zeros(5,nop); % allocate space for the output matrix and set to zero
    voltages(3,:) = odor_gate_on;  % write voltage for the odor gate. it is closed for low dose

if odor_or_mix ==1
    % increase the odor flow rate linearly

    voltages(1,1:end-1) = (1:nop-1)*flow_volt/nop;  % write voltage for the odor MFC.
    voltages(2,:) = 0;  % write voltage for the mixing MFC.
elseif odor_or_mix ==2
    % increase the odor flow rate linearly

    voltages(2,1:end-1) = (1:nop-1)*flow_volt/nop;  % write voltage for the odor MFC.
    voltages(1,:) = 0;  % write voltage for the mixing MFC.
else
    error('select 1: Odor or 2: mixture flow.')
end

    ControlParadigm(1).Name = 'SetFlow';
    ControlParadigm(1).Outputs = voltages;
    


    if saveit
            save ([save_name,'_Kontroller_Paradigm.mat'], 'ControlParadigm')
    end