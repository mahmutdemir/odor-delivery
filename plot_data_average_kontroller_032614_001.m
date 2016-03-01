clc;
if 1
CPrc = recombine_cparad(ControlParadigm);
data_temp = recombine_data(data);
data_orig = data;
data = data_temp;
CP_orig = ControlParadigm;
ControlParadigm = CPrc;
end
nofparad = length(data); % total number of paradigms including start & end
t_txt = '1pentan_ORN_PID_2014-04-09';
s_txt = '1pentan_1';
% inp_names = fieldnames(data);   % get input channel names
inp_names = {'PID'};   % set the input for the analysis manually
pri_time = 0.5; % seconds before the pulse applied
yliml = zeros(nofparad,2);  % array to determine the ylimit
ylimlb = zeros(nofparad,2);  % array to determine the ylimit
bck_val = zeros(nofparad,1);    % array to determine background values
meas_ti = 0.7;%sec
meas_te = 1;%sec
ave_pulse_amp = zeros(nofparad,1);    % average pulse amplitude measured between meas_ti and meas_te
ave_pulse_ampb = zeros(nofparad,1);    % average pulse amplitude measured between meas_ti and meas_te, background removed
err_pulse_amp = zeros(nofparad,1);    % average pulse amplitude measured between meas_ti and meas_te, standard deviation
% err_pulse_ampb = zeros(nofparad,1);    % standard deviation, average pulse amplitude measured between meas_ti and meas_te, background removed
err_bck_val = zeros(nofparad,1);    % average pulse amplitude measured between meas_ti and meas_te, standard deviation

% plots inputs channel by channel
for inp = 1:length(inp_names)
    figure('units','normalized','outerposition',[0 0 .8 .9]);
    % first figure out the y limits
    for prdgm = 1:nofparad
        yliml(prdgm,:)=([min(min(data(prdgm).(inp_names{inp}))) max(max(data(prdgm).(inp_names{inp})))]);
        bck_val(prdgm) = mean(mean(data(prdgm).(inp_names{inp})(2:end,1:pri_time*SamplingRate-1),1));
        err_bck_val(prdgm) = std(mean(data(prdgm).(inp_names{inp})(2:end,1:pri_time*SamplingRate-1),2),0,1);
        ylimlb(prdgm,:)= ([min(min(data(prdgm).(inp_names{inp}))) max(max(data(prdgm).(inp_names{inp})))])-bck_val(prdgm);
        ave_pulse_amp(prdgm) = mean(mean(data(prdgm).(inp_names{inp})(2:end,meas_ti*SamplingRate:meas_te*SamplingRate),1));
        err_pulse_amp(prdgm) = std(mean(data(prdgm).(inp_names{inp})(2:end,meas_ti*SamplingRate:meas_te*SamplingRate),2),0,1);
        ave_pulse_ampb(prdgm) = ave_pulse_amp(prdgm)-bck_val(prdgm);
    end
    ylimm = [99/100*min(yliml(:,1)) 101/100*max(yliml(:,2))];
    for prdgm = 1:nofparad
        [noftrial, ldata] = size(data(prdgm).(inp_names{inp}));
        % plot all curves together for different paradigms
        subplot(3,nofparad,prdgm)
        plot((1:ldata)/SamplingRate-pri_time,data(prdgm).(inp_names{inp}))
        del_ind = strfind(ControlParadigm(prdgm).Name,'_');
        if isempty(del_ind)
            ttxt = ControlParadigm(prdgm).Name;
        else
            ttxt = ControlParadigm(prdgm).Name(del_ind(1)+1:end);
        end
        title(ttxt,'Interpret','None')
%         xlabel('time (sec)')
        xlim([-pri_time ldata/SamplingRate-pri_time])
        ylim(ylimm)
        if prdgm==1
            ylabel([inp_names{inp} ' - All Curves'],'Interpret','None')
        end
        % now plot averages, throw the first one
        subplot(3,nofparad,nofparad+prdgm)
        plot((1:ldata)/SamplingRate-pri_time,mean(data(prdgm).(inp_names{inp})(2:end,:),1))

%         del_ind = strfind(ControlParadigm(prdgm).Name,'_');
%         if isempty(del_ind)
%             ttxt = ControlParadigm(prdgm).Name;
%         else
%             ttxt = ControlParadigm(prdgm).Name(del_ind(1)+1:end);
%         end
%         title(ttxt,'Interpret','None')
%         xlabel('time (sec)')
        xlim([-pri_time ldata/SamplingRate-pri_time])
        ylim(ylimm)
        if prdgm==1
            ylabel([inp_names{inp} ' - Average'],'Interpret','None')
        end
        % now plot averages minus the background
        subplot(3,nofparad,2*nofparad+prdgm)
        plot((1:ldata)/SamplingRate-pri_time,mean(data(prdgm).(inp_names{inp})(2:end,:),1)-bck_val(prdgm))
        bbot = min(mean(data(prdgm).(inp_names{inp})(2:end,meas_ti*SamplingRate:meas_te*SamplingRate),1))-bck_val(prdgm);
        bhgt = max(mean(data(prdgm).(inp_names{inp})(2:end,meas_ti*SamplingRate:meas_te*SamplingRate),1))-bck_val(prdgm)-bbot;
        if isnan(bbot) bbot = 0.1; end
        if isnan(bhgt) bhgt = 0.1; end
        rectangle('position',[.2 9.9/10*bbot .3 15/10*bhgt],'edgecolor','r','LineWidth',2)
%         rectangle('position',[.2 0 .3 1],'edgecolor','r','LineWidth',2)
%         del_ind = strfind(ControlParadigm(prdgm).Name,'_');
%         if isempty(del_ind)
%             ttxt = ControlParadigm(prdgm).Name;
%         else
%             ttxt = ControlParadigm(prdgm).Name(del_ind(1)+1:end);
%         end
%         title(ttxt,'Interpret','None')
        ylimmb = [1001/1000*min(ylimlb(:,1)) 1001/1000*max(ylimlb(:,2))];
        xlabel('time (sec)')
        xlim([-pri_time ldata/SamplingRate-pri_time])
        ylim(ylimmb)
        if prdgm==1
            ylabel([inp_names{inp} ' - Average - Bckg Rem.'],'Interpret','None')
        end
    end
    figname=[t_txt ' - plots'];saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r300',figname)
    figure('units','normalized','outerposition',[0 0 .8 .9]);
    subplot(1,2,1)
    errorbar(bck_val,err_bck_val,'or')
    hold on
    errorbar(ave_pulse_amp,err_pulse_amp,'*k')
    title([inp_names{inp} ' - Average Pulse Amplitude'],'Interpret','None')
    legend('background','amplitude')
    subplot(1,2,2)
    errorbar(ave_pulse_ampb,err_pulse_amp,'*k')
    title([inp_names{inp} ' - Average Pulse Amplitude - Background Removed'],'Interpret','None')
    c_odor = zeros(1,nofparad-2);
    for i=2:nofparad-1
        c = strsplit(ControlParadigm(i).Name,'_');
        if  strcmp(c(2),'pure');
            c_odor(i-1) = 100;
        else            
            c_odor(i-1) = str2double(c(2))/(str2double(c(2))+str2double(c(3)))*100;
        end
    end
    figname=[t_txt ' - amplitudes'];saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r300',figname)
    figure('units','normalized','outerposition',[0 0 .5 .7]);
    errorbar(c_odor,ave_pulse_ampb(2:end-1),err_pulse_amp(2:end-1),'*k')
    xlabel('odor %')
    ylabel('PID Amp. (volt)')
    xlim([-2 102])
%     ylim([-101/100*abs(ave_pulse_ampb(2)-err_pulse_amp(2)) 101/100*(ave_pulse_ampb(end-1)+err_pulse_amp(end-1))])
    ylim([min(ave_pulse_ampb)-999/1000*min(ave_pulse_ampb) 1001/1000*(max(ave_pulse_ampb)+max(err_pulse_amp))])
    title(t_txt,'Interpret','None')
    figname=[t_txt ' - Dose vs Odor Conc'];saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r300',figname)
            
end
%rename analyzed values
eval(['bck_val_' s_txt '= bck_val;']);
eval(['err_bck_val_' s_txt '= err_bck_val;']);
eval(['ave_pulse_amp_' s_txt '= ave_pulse_amp;']);
eval(['ave_pulse_ampb_' s_txt '= ave_pulse_ampb;']);
eval(['err_pulse_amp_' s_txt '= err_pulse_amp;']);
eval(['c_odor_' s_txt '= c_odor;']);
