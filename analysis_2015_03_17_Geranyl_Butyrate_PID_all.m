%% PID measurements of Geranyl Butyrate
% 03/18/2015, Mahmut Demir
%% Experiment and Data
% 3 ml Pure Geranyl Butyrate was loaded in glass bottles. Odor was 
% diluted in gas phase. The fow through the odor was on all times. (200
% ml/min)
% 
%  Data is located in newton>mahmut>EP Rig>2014_03_(16-17) Geranyl Butyrate PID
%
clear;
clc;
cd('C:\Users\md762\Desktop\Data\EP Rig\2015_03_16 Geranyl Butyrate PID');
d1 = load('2015_03_16_OBP_Geranyl_Butyrate_PID_1.mat');
d2 = load('2015_03_16_OBP_Geranyl_Butyrate_PID_2.mat');
cd('C:\Users\md762\Desktop\Data\EP Rig\2015_03_17 Geranyl Butyrate PID');
d3 = load('2015_03_17_OBP_Granyl_Butyrate_PID_1.mat');
d4 = load('2015_03_17_OBP_Granyl_Butyrate_PID_2.mat'); % has long recordings 11-12-13 300 sec long
d5 = load('2015_03_17_OBP_Granyl_Butyrate_PID_3.mat'); % delete 300sec long empty lines
d6 = load('2015_03_18_OBP_Geranyl_Butyrate_PID_1.mat'); % 10 min long step
d5.data(11:13) = [];
d5.ControlParadigm(11:13) = [];

%% Conclusion
% Geranyl Butyrate showed variable responses in PID. This odorant is a slow
% odorant and stabilizes in more than 10 minutes. Because of this slow
% stabilizing time, it is very diffucult to generate repeatable puffs of gas
% dilutions and flickering stimulus within reasonable duratiosn to perform experiments.
%
% We need to either change the odor or try again with a newer Geranyl
% Butyrate. Other options could be:
%%
% 
% * Ethyl Acetate (fast - weak response)
% * 2-3 Butanediol (slow - strong response)
% 


%% Results


%% 50 ms Pulses
%
% The first trial is higher than the rest. It becomes more repeatable in
% later trials
%  
clc;
flag = 5; % odor_valve
ControlParadigm = d2.ControlParadigm;
sr = d1.SamplingRate;
der_data = diff(ControlParadigm(2).Outputs(flag,:));
ind_tonoff = find(der_data~=0);
box_xloc = ind_tonoff/sr;
pri_time = ind_tonoff(1)/sr; % seconds before the pulse applied
stimon = pri_time;
if length(ind_tonoff)>2
    pulse_dur = (ind_tonoff(end)-ind_tonoff(1)-1)/sr;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(end);%sec
else
    pulse_dur = (ind_tonoff(2)-ind_tonoff(1)-1)/sr;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(2);%sec
end

time = (1:length(d1.data(2).PID(1,:)))/sr;
figure('units','normalized','outerposition',[0 0 .7 .8]);
prdglist = 2:4;
gmax = 0;
tempdata = [d1.data(prdglist(1)).PID ; d2.data(prdglist(1)).PID ; ...
        d3.data(prdglist(1)).PID ; d5.data(prdglist(1)).PID];
dil_val = zeros(length(prdglist),1);
lgnd_list = cell(length(prdglist),1);
pid_maxval = zeros(length(prdglist),size(tempdata,1));
for prdgms = prdglist
    data_temp = [d1.data(prdgms).PID ; d2.data(prdgms).PID ; ...
        d3.data(prdgms).PID ; d5.data(prdgms).PID];
    data_br = zeros(size(data_temp));
    max_val = zeros(1,size(data_temp,1));
    for j = 1:size(data_temp,1)
        data_br(j,:) = data_temp(j,:) - mean(data_temp(j,20000:pri_time*sr));
        max_val(j) = max(data_br(j,:));
    end
    subplot(1,3,prdgms-prdglist(1)+1)
    maxp = max(max_val);
    gmax = max(gmax,maxp);
    area([stimon stimon+pulse_dur], 1.1*[maxp maxp], 'faceColor', [.85 .85 .85], 'edgeColor', [.85 .85 .85]);
    hold on
    plot(time,data_br)
    extra_window = .1;
    xlim([pri_time-extra_window pri_time+pulse_dur+extra_window+.5])
    ylim([-.005 .8])
    nameparts = strsplit(ControlParadigm(prdgms).Name,'_');
    csub = strsplit(nameparts{end},'/');
    dil_val(prdgms-prdglist(1)+1) = str2double(csub{1})/(str2double(csub{1})+str2double(csub{2}));
    pid_maxval(prdgms-prdglist(1)+1,:) = max_val;
    title(nameparts{end},'Interpret','None')
    lgnd_list{prdgms-prdglist(1)+1} = nameparts{end};
    if prdgms-prdglist(1)+1==1
    xlabel('time (sec)')
    ylabel('voltage (V)')
    end
    PrettyFig

    
end
%%
% The dose values stabilize after 15 trials (note the 5 min approximate time
% difference between trials)
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
semilogy(1:size(tempdata,1),pid_maxval,'--o')
xlabel('Trial #')
ylabel('Max PID (V)')
title([num2str(pulse_dur),' sec pulses'])
legend(lgnd_list)
PrettyFig

%% 500 ms Pulses
%
% The shape of the lowest dose is not smooth probably. 
%  
clc;
flag = 5; % odor_valve
prdglist = 5:7;
der_data = diff(ControlParadigm(prdglist(1)).Outputs(flag,:));
ind_tonoff = find(der_data~=0);
box_xloc = ind_tonoff/sr;
pri_time = ind_tonoff(1)/sr; % seconds before the pulse applied
stimon = pri_time;
if length(ind_tonoff)>2
    pulse_dur = (ind_tonoff(end)-ind_tonoff(1)-1)/sr;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(end);%sec
else
    pulse_dur = (ind_tonoff(2)-ind_tonoff(1)-1)/sr;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(2);%sec
end

time = (1:length(d1.data(prdglist(1)).PID(1,:)))/sr;
figure('units','normalized','outerposition',[0 0 .7 .8]);
gmax = 0;
tempdata = [d1.data(prdglist(1)).PID ; d2.data(prdglist(1)).PID ; ...
        d3.data(prdglist(1)).PID ; d5.data(prdglist(1)).PID];
dil_val = zeros(length(prdglist),1);
lgnd_list = cell(length(prdglist),1);
pid_maxval = zeros(length(prdglist),size(tempdata,1));
for prdgms = prdglist
    data_temp = [d1.data(prdgms).PID ; d2.data(prdgms).PID ; ...
        d3.data(prdgms).PID ; d5.data(prdgms).PID];
    data_br = zeros(size(data_temp));
    max_val = zeros(1,size(data_temp,1));
    for j = 1:size(data_temp,1)
        data_br(j,:) = data_temp(j,:) - mean(data_temp(j,20000:pri_time*sr));
        max_val(j) = max(data_br(j,:));
    end
    subplot(1,3,prdgms-prdglist(1)+1)
    maxp = max(max_val);
    gmax = max(gmax,maxp);
    area([stimon stimon+pulse_dur], 1.1*[maxp maxp], 'faceColor', [.85 .85 .85], 'edgeColor', [.85 .85 .85]);
    hold on
    plot(time,data_br)
    extra_window = .1;
    xlim([pri_time-extra_window pri_time+pulse_dur+extra_window+.5])
    ylim([-.005 .7])
    nameparts = strsplit(ControlParadigm(prdgms).Name,'_');
    csub = strsplit(nameparts{end},'/');
    dil_val(prdgms-prdglist(1)+1) = str2double(csub{1})/(str2double(csub{1})+str2double(csub{2}));
    pid_maxval(prdgms-prdglist(1)+1,:) = max_val;
    title(nameparts{end},'Interpret','None')
    lgnd_list{prdgms-prdglist(1)+1} = nameparts{end};
    if prdgms-prdglist(1)+1==1
    xlabel('time (sec)')
    ylabel('voltage (V)')
    end
    PrettyFig


    
end
%%
% The doses become stable after 15 trials (~30 min) 
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
semilogy(1:size(tempdata,1),pid_maxval,'--o')
xlabel('Trial #')
ylabel('Max PID (V)')
title([num2str(pulse_dur),' sec pulses'])
legend(lgnd_list)
PrettyFig

%% 60 sec Pulses
%
% Middle dose decreases more sharply then others. 
%  
clc;
flag = 5; % odor_valve
prdglist = 8:10;
der_data = diff(ControlParadigm(prdglist(1)).Outputs(flag,:));
ind_tonoff = find(der_data~=0);
box_xloc = ind_tonoff/sr;
pri_time = ind_tonoff(1)/sr; % seconds before the pulse applied
stimon = pri_time;
if length(ind_tonoff)>2
    pulse_dur = (ind_tonoff(end)-ind_tonoff(1)-1)/sr;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(end);%sec
else
    pulse_dur = (ind_tonoff(2)-ind_tonoff(1)-1)/sr;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(2);%sec
end

time = (1:length(d2.data(prdglist(1)).PID(1,:)))/sr;
figure('units','normalized','outerposition',[0 0 .7 .8]);
gmax = 0;
tempdata = [d2.data(prdglist(1)).PID ; ...
        d3.data(prdglist(1)).PID ; d5.data(prdglist(1)).PID];
dil_val = zeros(length(prdglist),1);
lgnd_list = cell(length(prdglist),1);
pid_maxval = zeros(length(prdglist),size(tempdata,1));
for prdgms = prdglist
    data_temp = [d2.data(prdgms).PID ; ...
        d3.data(prdgms).PID ; d5.data(prdgms).PID];
    data_br = zeros(size(data_temp));
    max_val = zeros(1,size(data_temp,1));
    for j = 1:size(data_temp,1)
        data_br(j,:) = data_temp(j,:) - mean(data_temp(j,20000:pri_time*sr));
        max_val(j) = max(data_br(j,:));
    end
    subplot(1,3,prdgms-prdglist(1)+1)
    maxp = max(max_val);
    gmax = max(gmax,maxp);
    area([stimon stimon+pulse_dur], 1.1*[maxp maxp], 'faceColor', [.85 .85 .85], 'edgeColor', [.85 .85 .85]);
    hold on
    plot(time,data_br)
    extra_window = 5;
    xlim([pri_time-extra_window pri_time+pulse_dur+extra_window+.5])
    ylim([-.005 .35])
    nameparts = strsplit(ControlParadigm(prdgms).Name,'_');
    csub = strsplit(nameparts{end},'/');
    dil_val(prdgms-prdglist(1)+1) = str2double(csub{1})/(str2double(csub{1})+str2double(csub{2}));
    pid_maxval(prdgms-prdglist(1)+1,1:length(max_val)) = max_val;
    title(nameparts{end},'Interpret','None')
    lgnd_list{prdgms-prdglist(1)+1} = nameparts{end};
    if prdgms-prdglist(1)+1==1
    xlabel('time (sec)')
    ylabel('voltage (V)')
    end
    PrettyFig


    
end
%%
% The doses become stable after 15 trials (~30 min)
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
semilogy(1:size(tempdata,1),pid_maxval,'--o')
xlabel('Trial #')
ylabel('Max PID (V)')
title([num2str(pulse_dur),' sec pulses'])
legend(lgnd_list)
PrettyFig

%% Gaussian Noise
%
% Guassion noise is generated with 100 ms correlation length.
% the MFC flow is pretty stable
%  
clc;
prdglist = 11;
figure('units','normalized','outerposition',[0 0 .7 .8]);
time = (1:length(d2.data(prdglist).PID(1,:)))/sr;
data_temp_pid = [d2.data(prdglist(1)).PID ; ...
        d3.data(prdglist(1)).PID ; d5.data(prdglist(1)).PID];
data_temp_mfc = [d2.data(prdglist(1)).mfc500 ; ...
        d3.data(prdglist(1)).mfc500 ; d5.data(prdglist(1)).mfc500];
corr_mat = zeros(size(data_temp_pid,1),size(data_temp_pid,1));
for prdgms = prdglist
    for j = 1:size(data_temp_mfc,1)
        for k = 1:size(data_temp_mfc,1)
            temp_mat = corrcoef(data_temp_mfc(j,:),data_temp_mfc(k,:));
            corr_mat(j,k) = temp_mat(2);
        end
    end
end
subplot(3,3,[1:2,4:5,7:8])
plot(time,data_temp_mfc(1:end,:))
xlabel('time (sec)')
ylabel('MFC Signal (V)')
title('MFC Flow Signal')
xlim([30 35])
subplot(3,3,6)
imagesc(corr_mat)
colorbar
title('Corr. Coef.')

data_temp = data_temp_pid;

corr_mat = zeros(size(data_temp,1),size(data_temp,1));
for prdgms = prdglist
    for j = 1:size(data_temp,1)
        for k = 1:size(data_temp,1)
            temp_mat = corrcoef(data_temp(j,:),data_temp(k,:));
            corr_mat(j,k) = temp_mat(2);
        end
    end
end
%%
% The first trial is huge in the beginning and becomes stable in seconds.
% Trials after the first one have high correlation
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
subplot(3,3,[1:2,4:5,7:8])
plot(time,data_temp_pid(1:end,:))
xlabel('time (sec)')
ylabel('PID Voltage (V)')
title('Gaussian Noise')
subplot(3,3,6)
imagesc(corr_mat)
colorbar
title('Corr. Coef.')
%%
% A more close look at the data. Notice that the mean decreases at each
% trial
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);

plot(time,data_temp_pid(1:end,:))
xlabel('time (sec)')
ylabel('PID Voltage (V)')
title('Gaussian Noise - Detailed')
xlim([30 32])
ylim([0.15 .8])
PrettyFig
%%
% Stimulus distributions are not gaussian probably because this noise is
% not optimized for this odorant.
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
hold on
for prdgms = prdglist
    hist_val = zeros(200,size(data_temp_pid,1));
    x_val = zeros(200,size(data_temp_pid,1));
    for j = 1:size(data_temp,1)
%         hist(data_temp(j,:),200);
        [hist_val(:,j), x_val(:,j)] = hist(data_temp_pid(j,:),200);
        plot(x_val(:,j),hist_val(:,j));
        for k = j:size(data_temp,1)
            temp_mat = corrcoef(data_temp_pid(j,:),data_temp_pid(k,:));
            corr_mat(j,k) = temp_mat(2);
        end
    end
       
end
xlabel('PID (V)')
ylabel('count')
title('PID Histograms')
hold off

%% 10 min long step
%
% 10 min long step is recorded. It has a sharp transient maxiumum at the
% onset and decreases sharply after the onset. A slow increase follows that
% transient and then the PID value decreases up to 10 minutes with a slower
% rate.
% 
clc;
figure('units','normalized','outerposition',[0 0 .7 .8]);
plot((1:length(d6.data(4).PID))/sr,d6.data(4).PID,'k')
xlabel('time (sec)')
ylabel('PID (V)')
title('Geranyl Butyrate - 10 min step')
PrettyFig

%% Odor bottle content
% According to the label on the bottle this odor is not a monomolecular
% odor. It has other components as well and if the chemicals decomposed in
% time, we have more than two components in the bottle

botim = imread('C:\Users\md762\Desktop\Data\EP Rig\2015_03_17 Geranyl Butyrate PID\2015-03-18 13.32.21.jpg');
imshow(botim,'InitialMagnification',30)