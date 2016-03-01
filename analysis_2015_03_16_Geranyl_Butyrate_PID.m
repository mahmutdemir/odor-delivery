%% PID measurements of Geranyl Butyrate
% 03/17/2015, Mahmut Demir
%% Experiment and Data
% 3 ml Pure Geranyl Butyrate was loaded in glass bottles. Odor was 
% diluted in gas phase. The fow through the odor was on all times. (200
% ml/min)
% 
%  Data is located in newton>mahmut>EP Rig>2014_03_16 Geranyl Butyrate PID
% clear;
clc;
% cd('C:\Users\md762\Desktop\Data\EP Rig\2015_03_16 Geranyl Butyrate PID');
% load('2015_03_17_OBP_Granyl_Butyrate_PID_1.mat');
%% Conclusion
% Geranyl Butyrate showed variable responses in PID. The stimulus decreases
% in consecutive trials. The variablity among identical puffs seems to be
% high. Although the PID readout decreases in ceonsecutive trials the
% correscponding pid vs. time curves are similar with correlation
% coefficient above 95% excluding the first trial. The age of odor might
% be a reason for this outcome.

%% Results


%% 50 ms Pulses
%
% The first trial is higher than the rest. It becomes more repeatable in
% later trials
%  
clc;
flag = 5; % odor_valve
der_data = diff(ControlParadigm(2).Outputs(flag,:));
ind_tonoff = find(der_data~=0);
box_xloc = ind_tonoff/SamplingRate;
pri_time = ind_tonoff(1)/SamplingRate; % seconds before the pulse applied
stimon = pri_time;
sr = SamplingRate;
if length(ind_tonoff)>2
    pulse_dur = (ind_tonoff(end)-ind_tonoff(1)-1)/SamplingRate;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(end);%sec
else
    pulse_dur = (ind_tonoff(2)-ind_tonoff(1)-1)/SamplingRate;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(2);%sec
end

time = (1:length(data(2).PID(1,:)))/10000;
figure('units','normalized','outerposition',[0 0 .7 .8]);
prdglist = 2:4;
gmax = 0;
dil_val = zeros(3,1);
pid_val = zeros(3,1);
pid_err = zeros(3,1);
for prdgms = prdglist
    data_temp = data(prdgms).PID;
    data_br = zeros(size(data_temp));
    max_val = zeros(1,size(data_temp,1));
    for j = 1:size(data_temp,1)
        data_br(j,:) = data_temp(j,:) - mean(data_temp(j,1:pri_time*sr));
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
    ylim([-.005 .025])
    nameparts = strsplit(ControlParadigm(prdgms).Name,'_');
    csub = strsplit(nameparts{end},'/');
    dil_val(prdgms-prdglist(1)+1) = str2double(csub{1})/(str2double(csub{1})+str2double(csub{2}));
    pid_val(prdgms-prdglist(1)+1) = mean(max_val);
    pid_err(prdgms-prdglist(1)+1) = std(max_val);
    title(nameparts{end},'Interpret','None')
    if prdgms-prdglist(1)+1==1
    xlabel('time (sec)')
    ylabel('voltage (V)')
    end
    PrettyFig

    
end
%%
% Different doses are clearly different from each other on average
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
errorbar(dil_val*100,pid_val,pid_err,'ok')
xlabel('Dilution (%)')
ylabel('Average PID Voltage (V)')
title([num2str(pulse_dur),' sec pulses'])
PrettyFig

%% 500 ms Pulses
%
% The shape of the lowest dose is not smooth probably. 
%  
clc;
flag = 5; % odor_valve
der_data = diff(ControlParadigm(5).Outputs(flag,:));
ind_tonoff = find(der_data~=0);
box_xloc = ind_tonoff/SamplingRate;
pri_time = ind_tonoff(1)/SamplingRate; % seconds before the pulse applied
stimon = pri_time;
sr = SamplingRate;
if length(ind_tonoff)>2
    pulse_dur = (ind_tonoff(end)-ind_tonoff(1)-1)/SamplingRate;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(end);%sec
else
    pulse_dur = (ind_tonoff(2)-ind_tonoff(1)-1)/SamplingRate;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(2);%sec
end

time = (1:length(data(5).PID(1,:)))/10000;
figure('units','normalized','outerposition',[0 0 .7 .8]);
prdglist = 5:7;
gmax = 0;
dil_val = zeros(3,1);
pid_val = zeros(3,1);
pid_err = zeros(3,1);
for prdgms = prdglist
    data_temp = data(prdgms).PID;
    data_br = zeros(size(data_temp));
    max_val = zeros(1,size(data_temp,1));
    for j = 1:size(data_temp,1)
        data_br(j,:) = data_temp(j,:) - mean(data_temp(j,1:pri_time*sr));
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
    ylim([-.005 .025])
    nameparts = strsplit(ControlParadigm(prdgms).Name,'_');
    csub = strsplit(nameparts{end},'/');
    dil_val(prdgms-prdglist(1)+1) = str2double(csub{1})/(str2double(csub{1})+str2double(csub{2}));
    pid_val(prdgms-prdglist(1)+1) = mean(max_val);
    pid_err(prdgms-prdglist(1)+1) = std(max_val);
    title(nameparts{end},'Interpret','None')
    if prdgms-prdglist(1)+1==1
    xlabel('time (sec)')
    ylabel('voltage (V)')
    end
    PrettyFig

    
end
%%
% These are the averages of the doses and standard deviations 
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
errorbar(dil_val*100,pid_val,pid_err,'ok')
xlabel('Dilution (%)')
ylabel('Average PID Voltage (V)')
title([num2str(pulse_dur),' sec pulses'])
PrettyFig

%% 30 sec Pulses
%
% Middle dose decreases more sharply then aothers. 
%  
clc;
flag = 5; % odor_valve
der_data = diff(ControlParadigm(8).Outputs(flag,:));
ind_tonoff = find(der_data~=0);
box_xloc = ind_tonoff/SamplingRate;
pri_time = ind_tonoff(1)/SamplingRate; % seconds before the pulse applied
stimon = pri_time;
sr = SamplingRate;
if length(ind_tonoff)>2
    pulse_dur = (ind_tonoff(end)-ind_tonoff(1)-1)/SamplingRate;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(end);%sec
else
    pulse_dur = (ind_tonoff(2)-ind_tonoff(1)-1)/SamplingRate;
    meas_ti = box_xloc(1);%sec
    meas_te = box_xloc(2);%sec
end

time = (1:length(data(8).PID(1,:)))/10000;
figure('units','normalized','outerposition',[0 0 .7 .8]);
prdglist = 8:10;
gmax = 0;
dil_val = zeros(3,1);
pid_val = zeros(3,1);
pid_err = zeros(3,1);
for prdgms = prdglist
    data_temp = data(prdgms).PID;
    data_br = zeros(size(data_temp));
    max_val = zeros(1,size(data_temp,1));
    final_val = zeros(1,size(data_temp,1));
    for j = 1:size(data_temp,1)
        data_br(j,:) = data_temp(j,:) - mean(data_temp(j,1:pri_time*sr));
        max_val(j) = max(data_br(j,:));
        final_val(j) = mean(data_br(j,60*sr+1:65*sr+1));
    end
    subplot(1,3,prdgms-prdglist(1)+1)
    maxp = max(max_val);
    gmax = max(gmax,maxp);
    area([stimon stimon+pulse_dur], 1.1*[maxp maxp], 'faceColor', [.85 .85 .85], 'edgeColor', [.85 .85 .85]);
    hold on
    plot(time,data_br)
    extra_window = 10;
    xlim([pri_time-extra_window pri_time+pulse_dur+extra_window+.5])
    ylim([-.015 .05])
    nameparts = strsplit(ControlParadigm(prdgms).Name,'_');
    csub = strsplit(nameparts{end},'/');
    dil_val(prdgms-prdglist(1)+1) = str2double(csub{1})/(str2double(csub{1})+str2double(csub{2}));
    pid_val(prdgms-prdglist(1)+1) = mean(final_val);
    pid_err(prdgms-prdglist(1)+1) = std(final_val);
    title(nameparts{end},'Interpret','None')
    if prdgms-prdglist(1)+1==1
    xlabel('time (sec)')
    ylabel('voltage (V)')
    end
    PrettyFig

    
end
%%
% The variance increases as dose increases
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
errorbar(dil_val*100,pid_val,pid_err,'ok')
xlabel('Dilution (%)')
ylabel('Average PID Voltage (V)')
title([num2str(pulse_dur),' sec pulses'])
PrettyFig

%% Gaussian Noise
%
% Guassion noise is generated with 100 ms correlation length.
% the MFC flow is pretty stable
%  
clc;
prdglist = 14;
figure('units','normalized','outerposition',[0 0 .7 .8]);
time = (1:length(data(prdglist).PID(1,:)))/10000;
data_temp = data(prdglist).mfc500;
corr_mat = zeros(size(data_temp,1),size(data_temp,1));
for prdgms = prdglist
    data_temp = data(prdgms).mfc500;
    for j = 1:size(data_temp,1)
        for k = 1:size(data_temp,1)
            temp_mat = corrcoef(data_temp(j,:),data_temp(k,:));
            corr_mat(j,k) = temp_mat(2);
        end
    end
end
subplot(3,3,[1:2,4:5,7:8])
plot(time,data_temp(1:end,:))
xlabel('time (sec)')
ylabel('MFC Signal (V)')
title('MFC Flow Signal')
xlim([30 35])
subplot(3,3,6)
imagesc(corr_mat)
colorbar
title('Corr. Coef.')

data_temp = data(11).PID;

corr_mat = zeros(size(data_temp,1),size(data_temp,1));
for prdgms = prdglist
    data_temp = data(prdgms).PID;
    for j = 1:size(data_temp,1)
        for k = 1:size(data_temp,1)
            temp_mat = corrcoef(data_temp(j,:),data_temp(k,:));
            corr_mat(j,k) = temp_mat(2);
        end
    end
end
%%
% The firts trial is huge in the beginning and becomes stable in seconds.
% Trials after the fisrt one have high correlation
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
subplot(3,3,[1:2,4:5,7:8])
plot(time,data_temp(1:end,:))
xlabel('time (sec)')
ylabel('PID Voltage (V)')
title('Gaussian Noise')
subplot(3,3,6)
imagesc(corr_mat)
colorbar
title('Corr. Coef.')
%%
% A more close look at the data
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);

plot(time,data_temp(1:end,:))
xlabel('time (sec)')
ylabel('PID Voltage (V)')
title('Gaussian Noise - Detailed')
xlim([30 35])
ylim([0.15 .45])
PrettyFig
%%
% Stimulus distributions are not gaussian probably because this noise is
% not optimized for this odorant.
%  
figure('units','normalized','outerposition',[0 0 .7 .8]);
hold on
for prdgms = prdglist
    data_temp = data(prdgms).PID;
    hist_val = zeros(200,size(data_temp,1));
    x_val = zeros(200,size(data_temp,1));
    for j = 1:size(data_temp,1)
%         hist(data_temp(j,:),200);
        [hist_val(:,j), x_val(:,j)] = hist(data_temp(j,:),200);
        plot(x_val(:,j),hist_val(:,j));
        for k = j:size(data_temp,1)
            temp_mat = corrcoef(data_temp(j,:),data_temp(k,:));
            corr_mat(j,k) = temp_mat(2);
        end
    end
       
end
xlabel('PID (V)')
ylabel('count')
title('PID Histograms')
hold off



