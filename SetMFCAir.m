function SetMFCAir(Main,Odor,Mix,Gate)
%       SetMFCAir(Main,Odor,Mix,Gate)
% function SetMFCAir(Main,Odor,Mix,Gate)
% Sets the air flows of the MFC's and the odor gate.
% Main: 2L/min. 0: off, 1:on
% Odor: 0-100
% Mix; 0-100
% Gate: odor gate upstream off Low and high doses. 1: on, 0: off


ControlParadigm = make_ctrl_prdgms_018(Main,Odor,Mix,Gate,'',0);
data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10000);
clear all