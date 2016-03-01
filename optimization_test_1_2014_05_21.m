%% measure the matrix of paramters
%
if 0
clear all;
clc;
    save_date = '2014_05_21_boost_test_1';
    odor_rat = [0, .5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 18, 20];
    mix_rat = 20-odor_rat;
%     boost_width = [0, .1, .2, .5,  1, 2] ;
%     boost_volt = [0, .5, 1, 2, 3, 5];
    boost_width = [.1, .2, .3, .5, 1] ;
    boost_volt = [.5, 1, 2, 3, 5];
    total_rec_time= 10;
    pulse_width = .5; %sec
    init_time = 6; %sec
    odor_gate_on = 0; % at the end of the trial
    saveit = 0;
    
    mean_air = zeros(length(boost_width),length(boost_width));
    max_air = zeros(length(boost_width),length(boost_width));
    
    dil = 2;
    conv_factor = (odor_rat(1)+mix_rat(1))/20;
    odor_volt = conv_factor*odor_rat(dil)/(odor_rat(dil)+mix_rat(dil));   % voltage for odor MFC. This will be written to AO0
    set_flow  = odor_volt*200;
    timevec = (1:total_rec_time*10000)/10000;
    
    ControlParadigm = make_ctrl_prdgms_013(boost_volt(1),boost_width(1),odor_rat(1),mix_rat(1),init_time,pulse_width,total_rec_time,odor_gate_on,[save_date ''],saveit);
    data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10000);
    
    for i = 1:length(boost_width)
        for j = 1: length(boost_volt)
            ControlParadigm = make_ctrl_prdgms_013(boost_volt(j),boost_width(i),odor_rat(dil),mix_rat(dil),init_time,pulse_width,total_rec_time,odor_gate_on,[save_date ''],saveit);

            disp(['Doing the measurement of boost width: ' num2str(boost_width(i)) ' boost_volt: ' num2str(boost_volt(j))])
            data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',2,'w',10000);


            ind = find(data(2).Air_Odor(1,:)>odor_volt);
            plot(timevec,data(2).Air_Odor(1,:),'k',timevec(ind),data(2).Air_Odor(1,ind),'*r')
            mean_air(i,j) = (mean(data(2).Air_Odor(1,(init_time-2)*10000+1:(init_time)*10000+1)))*200;
            max_air(i,j) = max(data(2).Air_Odor(1,:)*200);

            disp(['Set Flow is: ' num2str(set_flow) ' where as'])
            disp(['Air Flow was: ' num2str(mean_air(i,j)) ' max flow: ' num2str(max_air(i,j))])
            
            pause(15)
            
        end
    end

    
    data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',3,'w',10000);
    
    boxplot(mean_air',boost_width,'notch','on')
end
%% find the best boost voltage for the given width
%
if 0
%clear all;
clc;
    save_date = '2014_05_22_boost_test_1';
    odor_rat = [0, .5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 18, 20];
    mix_rat = 20-odor_rat;
%     boost_width = [.1, .2, .3, .5, 1] ;
    boost_width = .5;
    total_rec_time= 10;
    pulse_width = .5; %sec
    init_time = 6; %sec
    odor_gate_on = 0; % at the end of the trial
    saveit = 0;
    
%     mean_air = zeros(length(odor_rat),length(boost_width),10);
%     max_air = zeros(length(odor_rat),length(boost_width),10);
%     boost_val = zeros(length(odor_rat),length(boost_width),10);
    
    
    
    
        conv_factor = (odor_rat(1)+mix_rat(1))/20;
        timevec = (1:total_rec_time*10000)/10000;
%         air_flow_mat = zeros(length(odor_rat),length(timevec));
        

%         boost_volt = 0;
%         ControlParadigm = make_ctrl_prdgms_013(boost_volt,boost_width(1),odor_rat(1),mix_rat(1),init_time,pulse_width,total_rec_time,odor_gate_on,[save_date ''],saveit);
%         data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10000);
%         figure; hold on
        
    for dil = 5: length(odor_rat)
        disp(['Dilution: ' num2str(odor_rat(dil)) '_' num2str(mix_rat(dil))])
        odor_volt = conv_factor*odor_rat(dil)/(odor_rat(dil)+mix_rat(dil));   % voltage for odor MFC. This will be written to AO0
        set_flow  = odor_volt*200;


        for i = 1:length(boost_width)
            
                goon = 1;
                iter = 1;
                iter1 = 0;
                
                if dil == 2
                   boost_volt = 3;
                end

                boost_voltmax = 3;
                boost_voltmin = 0;



            while goon

                ControlParadigm = make_ctrl_prdgms_013(boost_volt,boost_width(i),odor_rat(dil),mix_rat(dil),init_time,pulse_width,total_rec_time,odor_gate_on,[save_date ''],saveit);

                disp(['Doing the measurement of boost width: ' num2str(boost_width(i)) ' boost_volt: ' num2str(boost_volt)])
                data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',2,'w',10000);


                ind = find(data(2).Air_Odor(1,:)>odor_volt);
%                 plot(timevec,data(2).Air_Odor(1,:),'k',timevec(ind),data(2).Air_Odor(1,ind),'*r')
                mean_air(dil,i,iter:end) = (mean(data(2).Air_Odor(1,(init_time-2)*10000+1:(init_time)*10000+1)))*200;
                max_air(dil,i,iter:end) = max(data(2).Air_Odor(1,:)*200);
                
                disp(['Set Flow: ' num2str(set_flow)])
                disp(['This iteration is ' num2str(iter)])

                if mean_air(dil,i,iter) > set_flow+.5

                    boost_voltmax = boost_volt;
                    boost_voltdif = (boost_voltmax-boost_voltmin)/2;
                    boost_volt = boost_voltdif + boost_voltmin;
                    disp(['Set Flow < measured Flow: ' num2str(mean_air(dil,i,iter)) ' and also max flow was: ' num2str(max_air(dil,i,iter))])
                    disp('Decrease the boost voltage')
                    disp(['Set max volt to ' num2str(boost_voltmax)])
                    disp(['min volt = ' num2str(boost_voltmin)])
                    disp(['boost volt = ' num2str(boost_volt)])
                    iter1 = 0;


                elseif mean_air(dil,i,iter) < set_flow-.5

                    boost_voltmin = boost_volt;
                    boost_voltdif = (boost_voltmax-boost_voltmin)/2;
                    boost_volt = boost_voltdif + boost_voltmin;
                    disp(['Set Flow > measured Flow: ' num2str(mean_air(dil,i,iter)) ' and also max flow was: ' num2str(max_air(dil,i,iter))])
                    disp('Increase the boost voltage')
                    disp(['Set min volt to ' num2str(boost_voltmin)])
                    disp(['max volt = ' num2str(boost_voltmax)])
                    disp(['boost volt = ' num2str(boost_volt)])
                    iter1 = 0;
                    
                    if boost_voltdif == 0
                        goon = 0;
                    disp('This width is too less')
                end
                else
                    if iter1 == 0
                    disp('Air flow value is within 0.5 error range')
                    disp('Repeat thwo more times and accept if stable')
                    elseif iter1 == 2
                        goon = 0;
                        disp(['Stability Test - Iteration: ' num2str(iter1)])
                        disp('This value is stable, will be accepted')
                    else
                        disp(['Stability Test - Iteration: ' num2str(iter1)])
                        disp('Seems good...')
                    end
                    iter1 = iter1 +1;
                    
                end

                iter = iter + 1;
                if iter == 10
                    goon = 0;
                    disp('Max number of iterations is reached')
                end
                boost_val(dil,i,iter:end) = boost_volt;
                pause(15)

            end
            air_flow_mat(dil,:) = data(2).Air_Odor(1,:); 
        end
        subplot(4,4,dil)
        plot(squeeze(mean_air(dil,1,:)),'*k')
        title(['Dilution: ' num2str(odor_rat(dil)) '_' num2str(mix_rat(dil))])
        ylabel('air flow (ml/min)')
        xlabel('iter #')
    end


        data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',3,'w',10000);

end
%% Using the boost values from optimization now see how the system works 
% record three trials and see the look at the statistics
if 0
clc;
clear all
    save_date = '2014_05_26_boosted_highdose1';
    odor_rat = [0, .5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 18, 20];
    mix_rat = 20-odor_rat;
    boost_width = .5;
    total_rec_time= 15;
    pulse_width = .5; %sec
    init_time = 5; %sec
    odor_gate_on = 1; % at the end of the trial
    saveit = 1;
    load('2014_05_22_Boost_amplitudes')
    n_of_trials = 3;
    
    
     mean_air = zeros(length(odor_rat),n_of_trials);
     max_air = zeros(length(odor_rat),n_of_trials);
     boost_val = zeros(length(odor_rat),n_of_trials);
    
    
    
    
        conv_factor = (odor_rat(1)+mix_rat(1))/20;
        timevec = (1:total_rec_time*10000)/10000;
        air_flow_mat = zeros(length(odor_rat),length(timevec));
        

        ControlParadigm = make_ctrl_prdgms_014(boost_amp(:,3),boost_width,odor_rat,mix_rat,init_time,pulse_width,total_rec_time,odor_gate_on,[save_date '_off'],saveit);
        ControlParadigm = make_ctrl_prdgms_015(boost_amp(:,3),boost_width,odor_rat,mix_rat,init_time,pulse_width,total_rec_time,odor_gate_on,[save_date '_on'],saveit);
        ControlParadigm = make_ctrl_prdgms_016(boost_amp(:,3),boost_width,odor_rat,mix_rat,init_time,pulse_width,total_rec_time,odor_gate_on,[save_date '_onup200'],saveit);
%         data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10000);
%         figure; hold on
        
%     for dil = 1: length(odor_rat)
%         disp(['Dilution: ' num2str(odor_rat(dil)) '_' num2str(mix_rat(dil))])
%         odor_volt = conv_factor*odor_rat(dil)/(odor_rat(dil)+mix_rat(dil));   % voltage for odor MFC. This will be written to AO0
%         set_flow  = odor_volt*200;
% 
% 
%         for i = 1:length(boost_width)
%             
%                 goon = 1;
%                 iter = 1;
%                 iter1 = 0;
%                 
%                 if dil == 2
%                    boost_volt = 3;
%                 end
% 
%                 boost_voltmax = 3;
%                 boost_voltmin = 0;
% 
% 
% 
%             while goon
% 
%                 ControlParadigm = make_ctrl_prdgms_013(boost_volt,boost_width(i),odor_rat(dil),mix_rat(dil),init_time,pulse_width,total_rec_time,odor_gate_on,[save_date ''],saveit);
% 
%                 disp(['Doing the measurement of boost width: ' num2str(boost_width(i)) ' boost_volt: ' num2str(boost_volt)])
%                 data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',2,'w',10000);
% 
% 
%                 ind = find(data(2).Air_Odor(1,:)>odor_volt);
% %                 plot(timevec,data(2).Air_Odor(1,:),'k',timevec(ind),data(2).Air_Odor(1,ind),'*r')
%                 mean_air(dil,i,iter:end) = (mean(data(2).Air_Odor(1,(init_time-2)*10000+1:(init_time)*10000+1)))*200;
%                 max_air(dil,i,iter:end) = max(data(2).Air_Odor(1,:)*200);
%                 
%                 disp(['Set Flow: ' num2str(set_flow)])
%                 disp(['This iteration is ' num2str(iter)])
% 
%                 if mean_air(dil,i,iter) > set_flow+.5
% 
%                     boost_voltmax = boost_volt;
%                     boost_voltdif = (boost_voltmax-boost_voltmin)/2;
%                     boost_volt = boost_voltdif + boost_voltmin;
%                     disp(['Set Flow < measured Flow: ' num2str(mean_air(dil,i,iter)) ' and also max flow was: ' num2str(max_air(dil,i,iter))])
%                     disp('Decrease the boost voltage')
%                     disp(['Set max volt to ' num2str(boost_voltmax)])
%                     disp(['min volt = ' num2str(boost_voltmin)])
%                     disp(['boost volt = ' num2str(boost_volt)])
%                     iter1 = 0;
% 
% 
%                 elseif mean_air(dil,i,iter) < set_flow-.5
% 
%                     boost_voltmin = boost_volt;
%                     boost_voltdif = (boost_voltmax-boost_voltmin)/2;
%                     boost_volt = boost_voltdif + boost_voltmin;
%                     disp(['Set Flow > measured Flow: ' num2str(mean_air(dil,i,iter)) ' and also max flow was: ' num2str(max_air(dil,i,iter))])
%                     disp('Increase the boost voltage')
%                     disp(['Set min volt to ' num2str(boost_voltmin)])
%                     disp(['max volt = ' num2str(boost_voltmax)])
%                     disp(['boost volt = ' num2str(boost_volt)])
%                     iter1 = 0;
%                     
%                     if boost_voltdif == 0
%                         goon = 0;
%                     disp('This width is too less')
%                 end
%                 else
%                     if iter1 == 0
%                     disp('Air flow value is within 0.5 error range')
%                     disp('Repeat thwo more times and accept if stable')
%                     elseif iter1 == 2
%                         goon = 0;
%                         disp(['Stability Test - Iteration: ' num2str(iter1)])
%                         disp('This value is stable, will be accepted')
%                     else
%                         disp(['Stability Test - Iteration: ' num2str(iter1)])
%                         disp('Seems good...')
%                     end
%                     iter1 = iter1 +1;
%                     
%                 end
% 
%                 iter = iter + 1;
%                 if iter == 10
%                     goon = 0;
%                     disp('Max number of iterations is reached')
%                 end
%                 boost_val(dil,i,iter:end) = boost_volt;
%                 pause(15)
% 
%             end
%             air_flow_mat(dil,:) = data(2).Air_Odor(1,:); 
%         end
%         subplot(4,4,dil)
%         plot(squeeze(mean_air(dil,1,:)),'*k')
%         title(['Dilution: ' num2str(odor_rat(dil)) '_' num2str(mix_rat(dil))])
%         ylabel('air flow (ml/min)')
%         xlabel('iter #')
%     end
% 
% 
%         data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',3,'w',10000);

end

%% Generate doses from both of the bottles
% record three trials and look at the statistics
if 1
clc;
clear all
    save_date = '2014_07_04_EA_dose';
%     odor_rat = [{[0, 3, 8, 20]} {[3, 8, 20]}];
%     odor_rat = [{[0, 1, 2, 3, 5, 8, 10, 15, 20]} {[]}];
    odor_rat = [{[0 1 3 5 10 15 20]} {[0,2,3,4,5,6,8,10,15,20]}];
    total_flow = 20;
    boost_width = .5;
    total_rec_time= 10;
    pulse_width = .5; %sec
    init_time = 5; %sec
    saveit = 1;
    load('2014_05_22_Boost_amplitudes')
    n_of_trials = 3;
    
    
 
    ControlParadigm = make_ctrl_prdgms_017(boost_amp(:,4),boost_width,odor_rat,total_flow,0,init_time,pulse_width,total_rec_time,[save_date 'airoff'],saveit);
    ControlParadigm = make_ctrl_prdgms_017(boost_amp(:,4),boost_width,odor_rat,total_flow,1,init_time,pulse_width,total_rec_time,[save_date 'airon'],saveit);

    clear all
    %         data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',10000);
%         figure; hold on
        
%     for dil = 1: length(odor_rat)
%         disp(['Dilution: ' num2str(odor_rat(dil)) '_' num2str(mix_rat(dil))])
%         odor_volt = conv_factor*odor_rat(dil)/(odor_rat(dil)+mix_rat(dil));   % voltage for odor MFC. This will be written to AO0
%         set_flow  = odor_volt*200;
% 
% 
%         for i = 1:length(boost_width)
%             
%                 goon = 1;
%                 iter = 1;
%                 iter1 = 0;
%                 
%                 if dil == 2
%                    boost_volt = 3;
%                 end
% 
%                 boost_voltmax = 3;
%                 boost_voltmin = 0;
% 
% 
% 
%             while goon
% 
%                 ControlParadigm = make_ctrl_prdgms_013(boost_volt,boost_width(i),odor_rat(dil),mix_rat(dil),init_time,pulse_width,total_rec_time,odor_gate_on,[save_date ''],saveit);
% 
%                 disp(['Doing the measurement of boost width: ' num2str(boost_width(i)) ' boost_volt: ' num2str(boost_volt)])
%                 data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',2,'w',10000);
% 
% 
%                 ind = find(data(2).Air_Odor(1,:)>odor_volt);
% %                 plot(timevec,data(2).Air_Odor(1,:),'k',timevec(ind),data(2).Air_Odor(1,ind),'*r')
%                 mean_air(dil,i,iter:end) = (mean(data(2).Air_Odor(1,(init_time-2)*10000+1:(init_time)*10000+1)))*200;
%                 max_air(dil,i,iter:end) = max(data(2).Air_Odor(1,:)*200);
%                 
%                 disp(['Set Flow: ' num2str(set_flow)])
%                 disp(['This iteration is ' num2str(iter)])
% 
%                 if mean_air(dil,i,iter) > set_flow+.5
% 
%                     boost_voltmax = boost_volt;
%                     boost_voltdif = (boost_voltmax-boost_voltmin)/2;
%                     boost_volt = boost_voltdif + boost_voltmin;
%                     disp(['Set Flow < measured Flow: ' num2str(mean_air(dil,i,iter)) ' and also max flow was: ' num2str(max_air(dil,i,iter))])
%                     disp('Decrease the boost voltage')
%                     disp(['Set max volt to ' num2str(boost_voltmax)])
%                     disp(['min volt = ' num2str(boost_voltmin)])
%                     disp(['boost volt = ' num2str(boost_volt)])
%                     iter1 = 0;
% 
% 
%                 elseif mean_air(dil,i,iter) < set_flow-.5
% 
%                     boost_voltmin = boost_volt;
%                     boost_voltdif = (boost_voltmax-boost_voltmin)/2;
%                     boost_volt = boost_voltdif + boost_voltmin;
%                     disp(['Set Flow > measured Flow: ' num2str(mean_air(dil,i,iter)) ' and also max flow was: ' num2str(max_air(dil,i,iter))])
%                     disp('Increase the boost voltage')
%                     disp(['Set min volt to ' num2str(boost_voltmin)])
%                     disp(['max volt = ' num2str(boost_voltmax)])
%                     disp(['boost volt = ' num2str(boost_volt)])
%                     iter1 = 0;
%                     
%                     if boost_voltdif == 0
%                         goon = 0;
%                     disp('This width is too less')
%                 end
%                 else
%                     if iter1 == 0
%                     disp('Air flow value is within 0.5 error range')
%                     disp('Repeat thwo more times and accept if stable')
%                     elseif iter1 == 2
%                         goon = 0;
%                         disp(['Stability Test - Iteration: ' num2str(iter1)])
%                         disp('This value is stable, will be accepted')
%                     else
%                         disp(['Stability Test - Iteration: ' num2str(iter1)])
%                         disp('Seems good...')
%                     end
%                     iter1 = iter1 +1;
%                     
%                 end
% 
%                 iter = iter + 1;
%                 if iter == 10
%                     goon = 0;
%                     disp('Max number of iterations is reached')
%                 end
%                 boost_val(dil,i,iter:end) = boost_volt;
%                 pause(15)
% 
%             end
%             air_flow_mat(dil,:) = data(2).Air_Odor(1,:); 
%         end
%         subplot(4,4,dil)
%         plot(squeeze(mean_air(dil,1,:)),'*k')
%         title(['Dilution: ' num2str(odor_rat(dil)) '_' num2str(mix_rat(dil))])
%         ylabel('air flow (ml/min)')
%         xlabel('iter #')
%     end
% 
% 
%         data = Kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',3,'w',10000);

end