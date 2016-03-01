%%
for i = 2:34
    time = (1:length(data(i).PID(1,:)))/10000;
    subplot(5,7,i-1)
    plot(time,data(i).PID')
    xlim([1.5 8])
    ylim([0 4])
    C = strsplit(ControlParadigm(i).Name,'_');
    title(['p1:' C{2} '/' C{3} ' - p2:' C{5} '/' C{6}],'Interpret', 'None')
end
%%
% for 2 sec pulse
clear
cd('C:\Users\md762\Desktop\Data\EP Rig\2014_11_04');
figure('units','normalized','outerposition',[0 0 .9 .9]);
load('2014_11_04_ethyl_butyrate_pulse_addition_2sec_pulse_1.mat');
for i = 2:34
    time = (1:length(data(i).PID(1,:)))/10000;
    subplot(5,7,i-1)
    plot(time,data(i).PID')
    xlim([1.5 8])
    ylim([0 4])
    C = strsplit(ControlParadigm(i).Name,'_');
    title(['[' C{2} ' + ' C{5} ']  lag:' num2str(str2double(C{8})/1000) 's'],'Interpret', 'None')
end
%%
% for 2 sec pulse
clear
cd('C:\Users\md762\Desktop\Data\EP Rig\2014_11_04');
figure('units','normalized','outerposition',[0 0 .9 .9]);
load('2014_11_04_ethyl_butyrate_pulse_addition_2sec_pulse_1.mat');
for i = 2:34
    time = (1:length(data(i).PID(1,:)))/10000;
    subplot(5,7,i-1)
    plot(time,mean(data(i).PID,1),'k')
    xlim([1.5 8])
    ylim([0 4])
    C = strsplit(ControlParadigm(i).Name,'_');
    title(['[' C{2} ' + ' C{5} ']  lag:' num2str(str2double(C{8})/1000) 's'],'Interpret', 'None')
end
figname='eth_but addition - 2 sec pulse - 2014_11_04';saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r900',figname)
%%
% for .5 sec pulse
clear
cd('C:\Users\md762\Desktop\Data\EP Rig\2014_11_04');
figure('units','normalized','outerposition',[0 0 1 1]);
load('2014_11_04_ethyl_butyrate_pulse_addition_1.mat');
for i = 2:34
    time = (1:length(data(i).PID(1,:)))/10000;
    subplot(5,7,i-1)
    plot(time,mean(data(i).PID,1),'k')
    xlim([4.5 8])
    ylim([0 4])
    C = strsplit(ControlParadigm(i).Name,'_');
    title(['[' C{2} ' + ' C{5} ']  lag:' num2str(str2double(C{8})/1000) 's'],'Interpret', 'None')
end
figname='eth_but addition - half sec pulse - 2014_11_04';saveas(gcf,figname,'fig');set(gcf,'PaperPositionMode','auto');print(gcf,'-djpeg','-r900',figname)
%%
figure;
for i = 2:16
    subplot(3,5,i-1)
    plot(time,mean(data(i).PID,1))
    xlim([4.5 6])
    ylim([0 6])
    C = strsplit(ControlParadigm(i).Name,'_');
    title(['p1:' C{2} '/' C{3} ' - p2:' C{5} '/' C{6}],'Interpret', 'None')
end
%%
figure;

% title(['Hello {\color{blue}World}']);
for i = 2:16
    subplot(3,5,i-1)
    plot(time,data(i).odor'/5*500,'*k')
    hold on;
    plot(time,data(i).bck'/5*200,'r')
    plot(time,data(i).odor_dil'/5*1000,'b')
    plot(time,data(i).bck_dil'/5*500,'g')
    hold off
    xlim([4.5 6])
    ylim([0 250])
    C = strsplit(ControlParadigm(i).Name,'_');
    title(['v1:' [C{2}] '/' C{3} ' - v2:' C{5} '/' C{6}],'Interpret', 'None')
    legend(C{2}, C{3}, C{5}, C{6})
end
%%
figure;

for i = 2:16
    subplot(3,5,i-1)
    plot(time,mean(data(i).odor,1)/5*500,'k')
    hold on;
    plot(time,mean(data(i).bck,1)/5*200,'r')
    plot(time,mean(data(i).odor_dil,1)/5*1000,'b')
    plot(time,mean(data(i).bck_dil,1)/5*500,'g')
    hold off
    xlim([4.5 6])
    ylim([0 250])
    C = strsplit(ControlParadigm(i).Name,'_');
    title(['v1:' [C{2}] '/' C{3} ' - v2:' C{5} '/' C{6}],'Interpret', 'None')
    legend(C{2}, C{3}, C{5}, C{6})
end

%%
%
figure;
for i = 2:16
    subplot(3,5,i-1)
    odr = mean(data(i).odor(:,40000:60000),2);
    plot(1,odr/5*500,'*k')
    hold on;
    bck = mean(data(i).bck(:,40000:60000),2);
    plot(2,bck/5*200,'r')
    odr_dil = mean(data(i).odor_dil(:,40000:60000),2);
    plot(3,odr_dil/5*1000,'b')
    bck_dil = mean(data(i).odor(:,40000:60000),2);
    plot(4,bck_dil/5*500,'g')
    hold off
    xlim([0 5])
    ylim([-50 250])
    C = strsplit(ControlParadigm(i).Name,'_');
    title(['v1:' [C{2}] '/' C{3} ' - v2:' C{5} '/' C{6}],'Interpret', 'None')
    legend(C{2}, C{3}, C{5}, C{6})
end