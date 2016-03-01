% this code uses max to get the amplitude
clc;
auto_load = 0;
if auto_load ==1;
    clear all;
    load('2014_05_02_1_octen_3_ol_3.mat');
    if length(ControlParadigm)==length(data)
        there_is_end = 1;
    else
        there_is_end = 0;
    end
end
t_txt = 'Geranyl_butyrate_Clean_Test';
savefigure = 0;
recombine = 0;
setczero=1;
czeroval = .1;
there_is_end = 1;
plot_box = 0;
box_xloc = [1.3 1.8];
plot_area = 1;
plot_aoi = 1;
aoival = [0.00653 0.06551 3.117];
aoitxt = [{'-4'},{'-3'},{'-1'}];
ylim_auto = 0; % 0: no, 1: yes, 3: you set
ylimset = [-.002 .007];
xlimm = [1 -1 2];
stimon = 5;
pulse_dur = 0.5;
pri_time = 1; % seconds before the pulse applied
meas_ti = box_xloc(1);%sec
meas_te = box_xloc(2);%sec

% s_txt = '2014_04_11_csm1_ab3_1_meth_but_2';
s_txt = t_txt;
if recombine
CPrc = recombine_cparad(ControlParadigm);
data_temp = recombine_data(data);
data_orig = data;
data = data_temp;
CP_orig = ControlParadigm;
ControlParadigm = CPrc;
end
nofparad = length(data); % total number of paradigms including start & end
% inp_names = fieldnames(data);   % get input channel names
inp_names = {'PID'};   % set the input for the analysis manually

yliml = zeros(nofparad,2);  % array to determine the ylimit
ylimlb = zeros(nofparad,2);  % array to determine the ylimit
bck_val = zeros(nofparad,1);    % array to determine background values
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

        if length(data(prdgm).(inp_names{inp})(1,:))<meas_te*SamplingRate
        yliml(prdgm,:)=NaN;
        bck_val(prdgm) = NaN;
        err_bck_val(prdgm) = NaN;
%         ylimlb(prdgm,:)= ([min(min(data(prdgm).(inp_names{inp}))) max(max(data(prdgm).(inp_names{inp})))])-bck_val(prdgm);
        ylimlb(prdgm,:)= yliml(prdgm,:)-bck_val(prdgm);
            ave_pulse_amp(prdgm) = mean(mean(data(prdgm).(inp_names{inp})(1:end,:),1));
            err_pulse_amp(prdgm) = std(mean(data(prdgm).(inp_names{inp})(1:end,:),2),0,1);
        else
        yliml(prdgm,:)=([min(min(data(prdgm).(inp_names{inp})(1:end,meas_ti*SamplingRate:meas_te*SamplingRate))) max(max(data(prdgm).(inp_names{inp})(1:end,meas_ti*SamplingRate:meas_te*SamplingRate)))]);
        bck_val(prdgm) = mean(mean(data(prdgm).(inp_names{inp})(1:end,1:pri_time*SamplingRate-1),1));
        err_bck_val(prdgm) = std(mean(data(prdgm).(inp_names{inp})(1:end,1:pri_time*SamplingRate-1),2),0,1);
%         ylimlb(prdgm,:)= ([min(min(data(prdgm).(inp_names{inp}))) max(max(data(prdgm).(inp_names{inp})))])-bck_val(prdgm);
        ylimlb(prdgm,:)= yliml(prdgm,:)-bck_val(prdgm);
            ave_pulse_amp(prdgm) = mean(max(data(prdgm).(inp_names{inp})(1:end,meas_ti*SamplingRate:meas_te*SamplingRate),[],2));
            err_pulse_amp(prdgm) = std(max(data(prdgm).(inp_names{inp})(1:end,meas_ti*SamplingRate:meas_te*SamplingRate),[],2),0,1);
        end
        ave_pulse_ampb(prdgm) = ave_pulse_amp(prdgm)-bck_val(prdgm);
    end
    ylimm = [99/100*min(yliml(:,1)) 101/100*max(yliml(:,2))];
    for prdgm = 1:nofparad
        [noftrial, ldata] = size(data(prdgm).(inp_names{inp}));
        % plot all curves together for different paradigms
        subplot(3,nofparad,prdgm)
        if plot_area ==1
            area([stimon-pri_time stimon-pri_time+pulse_dur], [ylimm(2) ylimm(2)], 'faceColor', [.85 .85 .85], 'edgeColor', [.85 .85 .85])
            hold on
        end
        plot((1:ldata)/SamplingRate-pri_time,data(prdgm).(inp_names{inp}))
        del_ind = strfind(ControlParadigm(prdgm).Name,'_');
        if isempty(del_ind)
            ttxt = ControlParadigm(prdgm).Name;
        else
            ttxt = ControlParadigm(prdgm).Name(del_ind(1)+1:end);
        end
        title(ttxt,'Interpret','None')
%         xlabel('time (sec)')
        if xlimm(1) ==1
            xlim(xlimm(2:3));
        else
            xlim([-pri_time ldata/SamplingRate-pri_time]);
        end
        if ylim_auto==0
            ylim(ylimm);
        elseif ylim_auto==3;
            ylim(ylimset);
        end
        if prdgm==1
            ylabel([inp_names{inp} ' - All Curves'],'Interpret','None')
        else set(gca,'YTickLabel',[])
        end
        % now plot averages, throw the first one
        subplot(3,nofparad,nofparad+prdgm)
        if plot_area ==1
            area([stimon-pri_time stimon-pri_time+pulse_dur], [ylimm(2) ylimm(2)], 'faceColor', [.85 .85 .85], 'edgeColor', [.85 .85 .85])
            hold on
        end
        plot((1:ldata)/SamplingRate-pri_time,mean(data(prdgm).(inp_names{inp})(1:end,:),1))

%         del_ind = strfind(ControlParadigm(prdgm).Name,'_');
%         if isempty(del_ind)
%             ttxt = ControlParadigm(prdgm).Name;
%         else
%             ttxt = ControlParadigm(prdgm).Name(del_ind(1)+1:end);
%         end
%         title(ttxt,'Interpret','None')
%         xlabel('time (sec)')
        if xlimm(1) ==1
            xlim(xlimm(2:3));
        else
            xlim([-pri_time ldata/SamplingRate-pri_time]);
        end
        if ylim_auto==0
            ylim(ylimm)
        elseif ylim_auto==3;
            ylimm = ylimset;
        end
        if prdgm==1
            ylabel([inp_names{inp} ' - Average'],'Interpret','None')
        else set(gca,'YTickLabel',[])
        end
        % now plot averages minus the background
        subplot(3,nofparad,2*nofparad+prdgm)
        if plot_area ==1
            area([stimon-pri_time stimon-pri_time+pulse_dur], [ylimm(2) ylimm(2)], 'faceColor', [.85 .85 .85], 'edgeColor', [.85 .85 .85])
            hold on
        end
        plot((1:ldata)/SamplingRate-pri_time,mean(data(prdgm).(inp_names{inp})(1:end,:),1)-bck_val(prdgm))
        if length(data(prdgm).(inp_names{inp})(1,:))<meas_te*SamplingRate
            bbot = min(mean(data(prdgm).(inp_names{inp})(1:end,:),1))-bck_val(prdgm);
            bhgt = max(mean(data(prdgm).(inp_names{inp})(1:end,:),1))-bck_val(prdgm)-bbot;
        else
            bbot = min(mean(data(prdgm).(inp_names{inp})(1:end,meas_ti*SamplingRate:meas_te*SamplingRate),1))-bck_val(prdgm);
            bhgt = max(mean(data(prdgm).(inp_names{inp})(1:end,meas_ti*SamplingRate:meas_te*SamplingRate),1))-bck_val(prdgm)-bbot;
        end

        if isnan(bbot) bbot = 0.1; end
        if isnan(bhgt) bhgt = 0.1; end
        rectangle('position',[meas_ti-pri_time 9.9/10*bbot meas_te-meas_ti 15/10*bhgt],'edgecolor','r','LineWidth',2)
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
        if xlimm(1) ==1
            xlim(xlimm(2:3));
        else
            xlim([-pri_time ldata/SamplingRate-pri_time]);
        end
        if ylim_auto==0
            ylim(ylimmb)
        elseif ylim_auto==3;
            ylim(ylimset);
        end
        if prdgm==1
            ylabel([inp_names{inp} ' - Average - Bckg Rem.'],'Interpret','None')
        else set(gca,'YTickLabel',[])
        end
    end
    if savefigure
    figname=[t_txt ' - plots'];saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r300',figname)
    end
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
    c_odor = zeros(1,nofparad-1);
    for i=2:nofparad
        c = strsplit(ControlParadigm(i).Name,'_');
        if length(c)>1
            if  strcmp(c(2),'pure');
                c_odor(i-1) = 100;
            else            
                c_odor(i-1) = str2double(c(2))/(str2double(c(2))+str2double(c(3)))*100;
            end
            if setczero == 1
                if c_odor(i-1) ==0
                    c_odor(i-1) = czeroval;
                end
            end
        end
    end
    if 1;%savefigure
    figname=[t_txt ' - amplitudes'];saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r300',figname)
    end
    
    
    
    figure('units','normalized','outerposition',[0 0 .5 .7]);
    errorbar(c_odor/100,abs(ave_pulse_ampb(2:end)),err_pulse_amp(2:end),'*k')
    xlabel('dilution')
    ylabel('PID Amp. (volt)')
%     xlim([-2 102])
%     ylim([min(ave_pulse_ampb-err_pulse_amp) max(ave_pulse_ampb+err_pulse_amp)])
%     ylim([min(ave_pulse_ampb)-8/10*min(ave_pulse_ampb) 12/10*(max(ave_pulse_ampb)+max(err_pulse_amp))])
    title(t_txt,'Interpret','None')
    if 0;%savefigure
    figname=[t_txt ' - Dose vs Odor Conc'];saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r300',figname)
    end
        
    figure('units','normalized','outerposition',[0 0 .7 .9]);
    % first figure out the y limits
    if there_is_end ==1
        doit = nofparad - 1;
    else doit = nofparad;
    end
    for prdgm = 2:doit
        [noftrial, ldata] = size(data(prdgm).(inp_names{inp}));
        % plot all curves together for different paradigms
        subplot(2,doit-1,prdgm-1)
        if plot_area ==1
            area([stimon-pri_time stimon-pri_time+pulse_dur], [ylimm(2) ylimm(2)], 'faceColor', [.85 .85 .85], 'edgeColor', [.85 .85 .85])
            hold on
        end
        plot((1:ldata)/SamplingRate-pri_time,mean(data(prdgm).(inp_names{inp})(1:end,:),1)-bck_val(prdgm))
        del_ind = strfind(ControlParadigm(prdgm).Name,'_');
        if isempty(del_ind)
            ttxt = ControlParadigm(prdgm).Name;
        else
            ttxt = ControlParadigm(prdgm).Name(del_ind(1)+1:end);
        end
        title(ttxt,'Interpret','None')
        if plot_box ==1
            if length(data(prdgm).(inp_names{inp})(1,:))<meas_te*SamplingRate
                bbot = min(mean(data(prdgm).(inp_names{inp})(1:end,:),1))-bck_val(prdgm);
                bhgt = max(mean(data(prdgm).(inp_names{inp})(1:end,:),1))-bck_val(prdgm)-bbot;
            else
                bbot = min(mean(data(prdgm).(inp_names{inp})(1:end,meas_ti*SamplingRate:meas_te*SamplingRate),1))-bck_val(prdgm);
                bhgt = max(mean(data(prdgm).(inp_names{inp})(1:end,meas_ti*SamplingRate:meas_te*SamplingRate),1))-bck_val(prdgm)-bbot;
            end
            if isnan(bbot) bbot = 0.1; end
            if isnan(bhgt) bhgt = 0.1; end
            rectangle('position',[meas_ti-pri_time 9.9/10*bbot meas_te-meas_ti 15/10*bhgt],'edgecolor','r','LineWidth',2)
        end
        if xlimm(1) ==1
            xlim(xlimm(2:3));
        else
        xlim([-pri_time ldata/SamplingRate-pri_time-2])
        end
        ylimmb = [1001/1000*min(ylimlb(:,1)) 1001/1000*max(ylimlb(:,2))];
        if ylim_auto==0
            ylim(ylimmb)
        elseif ylim_auto==3;
            ylim(ylimset);
        end


        if prdgm ==2
            ylabel('PID Amp. (volt)')
            xlabel('time (sec)')
        else set(gca,'YTickLabel',[])
        end
    end
    subplot(2,nofparad-2,nofparad-1:2*nofparad-4)
    errorbar(c_odor/100,abs(ave_pulse_ampb(2:end)),err_pulse_amp(2:end),'*k')
    if plot_aoi ==1
        aoimat = zeros(length(aoival),length(c_odor));
        for kl = 1:length(aoival)
            aoimat(kl,:) = aoival(kl);
            text((c_odor(kl)+c_odor(kl+1))/200,aoival(kl),aoitxt(kl),...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left',...
        'FontSize',8,'BackgroundColor','cyan')
        end
        hold on
        plot(c_odor/100,aoimat)
        hold off
    end
    xlabel('odor dilution')
    ylabel('PID Amp. (volt)')
%     xlim([-2 102])
%     ylim([min(ave_pulse_ampb-err_pulse_amp) max(ave_pulse_ampb+err_pulse_amp)])
%     ylim([min(ave_pulse_ampb)-8/10*min(ave_pulse_ampb) 12/10*(max(ave_pulse_ampb)+max(err_pulse_amp))])
    title(t_txt,'Interpret','None')
    if savefigure
    figname=[t_txt ' - PID Curves and Doses'];saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r300',figname)
    end
    
    
    
    
end
%rename analyzed values
eval(['bck_val_' s_txt '= bck_val;']);
eval(['err_bck_val_' s_txt '= err_bck_val;']);
eval(['ave_pulse_amp_' s_txt '= ave_pulse_amp;']);
eval(['ave_pulse_ampb_' s_txt '= ave_pulse_ampb;']);
eval(['err_pulse_amp_' s_txt '= err_pulse_amp;']);
eval(['c_odor_' s_txt '= c_odor;']);



%% plot two experiments together
%
if 0
    
clc;    
figure('units','normalized','outerposition',[0 0 .7 .9]);
    % first figure out the y limits
savefigure = 0;
plot_box = 0;
n_of_prd = length(Exp1.data);

    
    for prdgm = 2:n_of_prd-1
        [noftrial, ldata] = size(Exp1.data(prdgm).(Exp1.inp_names{Exp1.inp}));
        % plot all curves together for different paradigms
        subplot(2,n_of_prd-2,prdgm-1)
        plot((1:ldata)/Exp1.SamplingRate-Exp1.pri_time,mean(Exp1.data(prdgm).(Exp1.inp_names{Exp1.inp})(1:end,:),1)-Exp1.bck_val(prdgm),'k')
        hold on
        plot((1:ldata)/Exp2.SamplingRate-Exp2.pri_time,mean(Exp2.data(prdgm).(Exp2.inp_names{Exp2.inp})(1:end,:),1)-Exp2.bck_val(prdgm),'r')
        plot((1:ldata)/Exp3.SamplingRate-Exp3.pri_time,mean(Exp3.data(prdgm).(Exp3.inp_names{Exp3.inp})(1:end,:),1)-Exp3.bck_val(prdgm),'b')
        hold off
        del_ind = strfind(Exp1.ControlParadigm(prdgm).Name,'_');
        if isempty(del_ind)
            ttxt = Exp1.ControlParadigm(prdgm).Name;
        else
            ttxt = Exp1.ControlParadigm(prdgm).Name(del_ind(1)+1:end);
        end
        title(ttxt,'Interpret','None')
        if plot_box ==1
            if length(Exp1.data(prdgm).(Exp1.inp_names{Exp1.inp})(1,:))<meas_te*Exp1.SamplingRate
                bbot = min(mean(Exp1.data(prdgm).(Exp1.inp_names{Exp1.inp})(1:end,:),1))-Exp1.bck_val(prdgm);
                bhgt = max(mean(Exp1.data(prdgm).(Exp1.inp_names{Exp1.inp})(1:end,:),1))-Exp1.bck_val(prdgm)-bbot;
            else
                bbot = min(mean(Exp1.data(prdgm).(Exp1.inp_names{Exp1.inp})(1:end,meas_ti*Exp1.SamplingRate:meas_te*Exp1.SamplingRate),1))-Exp1.bck_val(prdgm);
                bhgt = max(mean(Exp1.data(prdgm).(Exp1.inp_names{Exp1.inp})(1:end,meas_ti*Exp1.SamplingRate:meas_te*Exp1.SamplingRate),1))-Exp1.bck_val(prdgm)-bbot;
            end
            if isnan(bbot) bbot = 0.1; end
            if isnan(bhgt) bhgt = 0.1; end
            rectangle('position',[meas_ti-pri_time 9.9/10*bbot meas_te-meas_ti 15/10*bhgt],'edgecolor','r','LineWidth',2)
        end
        ylimmb = [1001/1000*min(Exp1.ylimlb(:,1)) 1001/1000*max(Exp1.ylimlb(:,2))];
%         xlabel('time (sec)')
        if prdgm ==2
            ylabel('PID Amp. (volt)')
        end
        xlim([-Exp1.pri_time ldata/Exp1.SamplingRate-Exp1.pri_time-2])
        ylim(ylimmb)
        if prdgm==1
            ylabel([Exp1.inp_names{Exp1.inp} ' - Average - Bckg Rem.'],'Interpret','None')
        end
    end
    subplot(2,13,14:26)
    errorbar(Exp1.c_odor,Exp1.ave_pulse_ampb(2:end),Exp1.err_pulse_amp(2:end),'*k')
    hold on;
    errorbar(Exp2.c_odor,Exp2.ave_pulse_ampb(2:end),Exp2.err_pulse_amp(2:end),'*r')
    errorbar(Exp3.c_odor,Exp3.ave_pulse_ampb(2:end),Exp3.err_pulse_amp(2:end),'*b')
    hold off
    xlabel('odor %')
    ylabel('PID Amp. (volt)')
    xlim([-2 102])
    hleg1=legend('Exp-1','Exp-2','Exp-3');
    set(hleg1,'Location','NorthWest')
    set(hleg1,'Interpreter','none')
%     ylim([min(ave_pulse_ampb-err_pulse_amp) max(ave_pulse_ampb+err_pulse_amp)])
%     ylim([min(ave_pulse_ampb)-8/10*min(ave_pulse_ampb) 12/10*(max(ave_pulse_ampb)+max(err_pulse_amp))])
    title(Exp1.t_txt,'Interpret','None')
    if savefigure
    figname=[Exp1.t_txt ' - PID Curves and Doses'];saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r300',figname)
    end
    
    
end