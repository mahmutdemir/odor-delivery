function SetLinIncMFCAir(Odor_Mix,MaxFlow,Duration,Gate)
%       SetLinIncMFCAir(Odor_Mix,MaxFlow,Duration,Gate)
% function SetLinIncMFCAir(Odor_Mix,MaxFlow,Duration,Gate)
% Sets the air flows of the MFC's and the odor gate.
% Odor_Mix:      1: odor, 2: mix
% Max_Flow:      maximum flow rate to be reached
% Duration:      the duration for the increase process 
% Gate:          odor gate upstream off Low and high doses. 1: on, 0: off


ControlParadigm = make_ctrl_prdgms_023(Odor_Mix,MaxFlow,Duration,Gate,'',0);
data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10000);
clear all