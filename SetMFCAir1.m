function SetMFCAir1(Main,Odor,Mix,bckg)
%      SetMFCAir1(Main,Odor,Mix,bckg)
% function SetMFCAir1(Main,Odor,Mix,bckg)
% Sets the air flows of the MFC's
% Main: 0-5000;
% Odor: 0-500
% Mix; 0-500
% Mix; 0-200


ControlParadigm = make_ctrl_prdgms_030(Main,Odor,Mix,bckg);
data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10000);
clear all