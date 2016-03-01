%% Odorant onset time detection

winsdata = 1;
winsdiff = 1;
stdthreshold = 5;
txt ='_w11_st5';
percent_mat = [{'0%'} {'20%'} {'40%'} {'60%'} {'80%'} {'100%'}];
hpoint = zeros(6,5);
hpoint_std = zeros(6,5);

wsd =  [1 5 10 20 50 100 1000  1   5   10];
wsdf = [1 5 10 20 50 100 1000 100 100 100];
thr = [3 5 7 10];
meanton = zeros(length(wsd),4);
stdton = zeros(length(wsd),4);

for thrind = 1:1%length(thr)
    stdthreshold = thr(thrind);
    for lind = 1:1%length(wsd)
        winsdata = wsd(lind);
        winsdiff = wsdf(lind);
        txt = ['_w' num2str(winsdata) num2str(winsdiff) '_st' num2str(stdthreshold)];
                hpoint_std = zeros(6,5);

        for dil=1:6
            if dil==1
                figure('units','normalized','outerposition',[0 0 .9 .9]);
                for trial=1:5
                    subplot(2,3,trial)
                    std_noise = std(box_ave(diff(box_ave(PID(dil,trial,1:9000),100)),100));
                    indm = find(box_ave(diff(box_ave(PID(dil,trial,:),100)),100)>stdthreshold*std_noise);
                    indm(indm<10000)=[];
                    indm(indm>11000)=[];
                    if isempty(indm)
                        indm = 10000;
                    end
                    hpoint_std(dil,trial)= min(indm);
                    derdata = box_ave(diff(box_ave(squeeze(PID(dil,trial,:)),100)),100);
                    plot(timel(1:end-1),derdata,timel(indm),derdata(indm),'*r')
                    title('derivative')
                    xlim([1.04 1.13])
                end
                subplot(2,3,6)
                text(.1,.6,['\fontsize{16}windsize =' num2str(winsdata) num2str(winsdiff)])
                text(.1,.5,['\fontsize{16}std thresh = ' num2str(stdthreshold)])
                text(.1,.4,['\fontsize{16}dose = ' percent_mat{dil}])
                    figure('units','normalized','outerposition',[0 0 .9 .9]);
                for trial=1:5
                    subplot(2,3,trial)
                    std_noise = std(box_ave(diff(box_ave(PID(dil,trial,1:9000),100)),100));
                    indm = find(box_ave(diff(box_ave(PID(dil,trial,:),100)),100)>stdthreshold*std_noise);
                    indm(indm<10000)=[];
                    indm(indm>11000)=[];
                    if isempty(indm)
                        indm = 10000;
                    end
                    hpoint_std(dil,trial)= min(indm);
            %                 plot(timel,squeeze(PID(dil,trial,:)),timel(min(indm)+winsdata+winsdiff),PID(dil,trial,min(indm)+winsdata+winsdiff),'*r')
                    plot(timel,squeeze(PID(dil,trial,:)),timel(min(indm)),PID(dil,trial,min(indm)),'*r')
                    title('std calc')
                    xlim([1.04 1.13])
                end
                subplot(2,3,6)
                text(.1,.6,['\fontsize{16}windsize =' num2str(winsdata) num2str(winsdiff)])
                text(.1,.5,['\fontsize{16}std thresh = ' num2str(stdthreshold)])
                text(.1,.4,['\fontsize{16}dose = ' percent_mat{dil}])

            else

                figure('units','normalized','outerposition',[0 0 .9 .9]);
                for trial=1:5
                    subplot(2,3,trial)
                    std_noise = std(box_ave(diff(box_ave(PID(dil,trial,1:9000),winsdata)),winsdiff));
                    indm = find(box_ave(diff(box_ave(PID(dil,trial,:),winsdata)),winsdiff)>stdthreshold*std_noise);
                    indm(indm<10000)=[];
                    indm(indm>11000)=[];
                    if isempty(indm)
                        indm = 10000;
                    end
                    hpoint_std(dil,trial)= min(indm);
                    derdata = box_ave(diff(box_ave(squeeze(PID(dil,trial,:)),winsdata)),winsdiff);
                    plot(timel(1:end-1),derdata,timel(indm),derdata(indm),'*r')
                    title('derivative')
                    xlim([1.04 1.13])
                end
                subplot(2,3,6)
                text(.1,.6,['\fontsize{16}windsize =' num2str(winsdata) num2str(winsdiff)])
                text(.1,.5,['\fontsize{16}std thresh = ' num2str(stdthreshold)])
                text(.1,.4,['\fontsize{16}dose = ' percent_mat{dil}])
                    figure('units','normalized','outerposition',[0 0 .9 .9]);
                for trial=1:5
                    subplot(2,3,trial)
                    std_noise = std(box_ave(diff(box_ave(PID(dil,trial,1:9000),winsdata)),winsdiff));
                    indm = find(box_ave(diff(box_ave(PID(dil,trial,:),winsdata)),winsdiff)>stdthreshold*std_noise);
                    indm(indm<10000)=[];
                    indm(indm>11000)=[];
                    if isempty(indm)
                        indm = 10000;
                    end
                    hpoint_std(dil,trial)= min(indm);
            %                 plot(timel,squeeze(PID(dil,trial,:)),timel(min(indm)+winsdata+winsdiff),PID(dil,trial,min(indm)+winsdata+winsdiff),'*r')
                    plot(timel,squeeze(PID(dil,trial,:)),timel(min(indm)),PID(dil,trial,min(indm)),'*r')
                    title('std calc')
                    xlim([1.04 1.13])
                end
                subplot(2,3,6)
                text(.1,.6,['\fontsize{16}windsize =' num2str(winsdata) num2str(winsdiff)])
                text(.1,.5,['\fontsize{16}std thresh = ' num2str(stdthreshold)])
                text(.1,.4,['\fontsize{16}dose = ' percent_mat{dil}])
            end

        end
        t_hpoint = hpoint_std/10;
        t_hpoint(:,6) = mean(t_hpoint(:,1:5),2);
        t_hpoint(:,7) = std(t_hpoint(:,1:5),0,2);
        figure('units','normalized','outerposition',[0 0 .6 .7]);
        errorbar([0 20 40 60 80 100],t_hpoint(:,6),t_hpoint(:,7),'ko')
        xlabel('dose')
        grid on
        ylabel('odor onset time (ms)')
        % ylim([1040 1100])

        eval(['t_hpoint' txt '=t_hpoint;']);
    end
end


        % calculate the time from PID itself
        % for dil=1:6
        %     figure('units','normalized','outerposition',[0 0 .9 .9]);
        %     for trial=1:5
        %         subplot(2,3,trial)
        %         std_noise = std((box_ave(PID(dil,trial,1:9000),winsdata)));
        %         indm = find((box_ave(PID(dil,trial,:),winsdata))>stdthreshold*std_noise);
        %         indm(indm<10000)=[];
        %         indm(indm>11000)=[];
        %         if isempty(indm)
        %             indm = 10600;
        %         end
        %         hpoint_std(dil,trial)= min(indm);
        %         data = box_ave(squeeze(PID(dil,trial,:)),winsdata);
        %         plot(timel,data,timel(indm),data(indm),'g',timel(min(indm)),data(min(indm)),'*r')
        %         title('data')
        %         xlim([1.04 1.13])
        %     end
        %         
        % end
        % t_hpointd = hpoint_std/10;
        % t_hpointd(:,6) = mean(t_hpointd(:,1:5),2);
        % t_hpointd(:,7) = std(t_hpointd(:,1:5),0,2);
        % figure('units','normalized','outerposition',[0 0 .9 .9]);
        % errorbar([0 20 40 60 80 100],t_hpointd(:,6),t_hpointd(:,7),'ko')
        % xlabel('dose')
        % grid on
        % ylabel('time_hit_data (sec)')
        % ylim([1040 1100])
        
 for thrind = 1:length(thr)
    stdthreshold = thr(thrind);
    for lind = 1:length(wsd)
        winsdata = wsd(lind);
        winsdiff = wsdf(lind);
        txt = ['_w' num2str(winsdata) num2str(winsdiff) '_st' num2str(stdthreshold)];
                hpoint_std = zeros(6,5);

        for dil=1:6
            if dil==1
%                 figure('units','normalized','outerposition',[0 0 .9 .9]);
                for trial=1:5
%                     subplot(2,3,trial)
                    std_noise = std(box_ave(diff(box_ave(PID(dil,trial,1:9000),100)),100));
                    indm = find(box_ave(diff(box_ave(PID(dil,trial,:),100)),100)>stdthreshold*std_noise);
                    indm(indm<10000)=[];
                    indm(indm>11000)=[];
                    if isempty(indm)
                        indm = 10000;
                    end
                    hpoint_std(dil,trial)= min(indm);
                    derdata = box_ave(diff(box_ave(squeeze(PID(dil,trial,:)),100)),100);
%                     plot(timel(1:end-1),derdata,timel(indm),derdata(indm),'*r')
%                     title('derivative')
%                     xlim([1.04 1.13])
                end
%                 subplot(2,3,6)
%                 text(.1,.6,['\fontsize{16}windsize =' num2str(winsdata) num2str(winsdiff)])
%                 text(.1,.5,['\fontsize{16}std thresh = ' num2str(stdthreshold)])
%                 text(.1,.4,['\fontsize{16}dose = ' percent_mat{dil}])
%                     figure('units','normalized','outerposition',[0 0 .9 .9]);
                for trial=1:5
%                     subplot(2,3,trial)
%                     std_noise = std(box_ave(diff(box_ave(PID(dil,trial,1:9000),100)),100));
%                     indm = find(box_ave(diff(box_ave(PID(dil,trial,:),100)),100)>stdthreshold*std_noise);
                    indm(indm<10000)=[];
                    indm(indm>11000)=[];
                    if isempty(indm)
                        indm = 10000;
                    end
                    hpoint_std(dil,trial)= min(indm);
            %                 plot(timel,squeeze(PID(dil,trial,:)),timel(min(indm)+winsdata+winsdiff),PID(dil,trial,min(indm)+winsdata+winsdiff),'*r')
%                     plot(timel,squeeze(PID(dil,trial,:)),timel(min(indm)),PID(dil,trial,min(indm)),'*r')
%                     title('std calc')
%                     xlim([1.04 1.13])
                end
%                 subplot(2,3,6)
%                 text(.1,.6,['\fontsize{16}windsize =' num2str(winsdata) num2str(winsdiff)])
%                 text(.1,.5,['\fontsize{16}std thresh = ' num2str(stdthreshold)])
%                 text(.1,.4,['\fontsize{16}dose = ' percent_mat{dil}])

            else

%                 figure('units','normalized','outerposition',[0 0 .9 .9]);
                for trial=1:5
%                     subplot(2,3,trial)
                    std_noise = std(box_ave(diff(box_ave(PID(dil,trial,1:9000),winsdata)),winsdiff));
                    indm = find(box_ave(diff(box_ave(PID(dil,trial,:),winsdata)),winsdiff)>stdthreshold*std_noise);
                    indm(indm<10000)=[];
                    indm(indm>11000)=[];
                    if isempty(indm)
                        indm = 10000;
                    end
                    hpoint_std(dil,trial)= min(indm);
                    derdata = box_ave(diff(box_ave(squeeze(PID(dil,trial,:)),winsdata)),winsdiff);
%                     plot(timel(1:end-1),derdata,timel(indm),derdata(indm),'*r')
%                     title('derivative')
%                     xlim([1.04 1.13])
                end
%                 subplot(2,3,6)
%                 text(.1,.6,['\fontsize{16}windsize =' num2str(winsdata) num2str(winsdiff)])
%                 text(.1,.5,['\fontsize{16}std thresh = ' num2str(stdthreshold)])
%                 text(.1,.4,['\fontsize{16}dose = ' percent_mat{dil}])
%                     figure('units','normalized','outerposition',[0 0 .9 .9]);
                for trial=1:5
%                     subplot(2,3,trial)
                    std_noise = std(box_ave(diff(box_ave(PID(dil,trial,1:9000),winsdata)),winsdiff));
                    indm = find(box_ave(diff(box_ave(PID(dil,trial,:),winsdata)),winsdiff)>stdthreshold*std_noise);
                    indm(indm<10000)=[];
                    indm(indm>11000)=[];
                    if isempty(indm)
                        indm = 10000;
                    end
                    hpoint_std(dil,trial)= min(indm);
            %                 plot(timel,squeeze(PID(dil,trial,:)),timel(min(indm)+winsdata+winsdiff),PID(dil,trial,min(indm)+winsdata+winsdiff),'*r')
%                     plot(timel,squeeze(PID(dil,trial,:)),timel(min(indm)),PID(dil,trial,min(indm)),'*r')
%                     title('std calc')
%                     xlim([1.04 1.13])
                end
%                 subplot(2,3,6)
%                 text(.1,.6,['\fontsize{16}windsize =' num2str(winsdata) num2str(winsdiff)])
%                 text(.1,.5,['\fontsize{16}std thresh = ' num2str(stdthreshold)])
%                 text(.1,.4,['\fontsize{16}dose = ' percent_mat{dil}])
            end

        end
        t_hpoint = hpoint_std/10;
        t_hpoint(:,6) = mean(t_hpoint(:,1:5),2);
        t_hpoint(:,7) = std(t_hpoint(:,1:5),0,2);
        meanton(lind,thrind) = mean(t_hpoint(:,6));
        stdton(lind,thrind) = std(t_hpoint(:,6));
        figure('units','normalized','outerposition',[0 0 .6 .7]);
        errorbar([0 20 40 60 80 100],t_hpoint(:,6),t_hpoint(:,7),'ko')
        xlabel('dose')
        grid on
        ylabel('Odor onset time (ms)')
        % ylim([1040 1100])
        title(['windsize =' num2str(winsdata) ' - ' num2str(winsdiff) ', std thresh = ' num2str(stdthreshold)])
        eval(['t_hpoint' txt '=t_hpoint;']);
    end
 end

        figure('units','normalized','outerposition',[0 0 .8 .9]);
        errorbar((1:length(wsd)),meanton(:,1),stdton(:,1),'ko')
        hold on
        errorbar((1:length(wsd)),meanton(:,2),stdton(:,2),'ro')
        errorbar((1:length(wsd)),meanton(:,3),stdton(:,3),'bo')
        errorbar((1:length(wsd)),meanton(:,4),stdton(:,4),'go')
        xlabel('window size #')
        grid on
        ylabel('mean time (ms))')
        title('Mean odor onset times')
        legend('thr = 3','thr = 5','thr = 7','thr = 10')
        hold off

        % calculate the time from PID itself
        % for dil=1:6
        %     figure('units','normalized','outerposition',[0 0 .9 .9]);
        %     for trial=1:5
        %         subplot(2,3,trial)
        %         std_noise = std((box_ave(PID(dil,trial,1:9000),winsdata)));
        %         indm = find((box_ave(PID(dil,trial,:),winsdata))>stdthreshold*std_noise);
        %         indm(indm<10000)=[];
        %         indm(indm>11000)=[];
        %         if isempty(indm)
        %             indm = 10600;
        %         end
        %         hpoint_std(dil,trial)= min(indm);
        %         data = box_ave(squeeze(PID(dil,trial,:)),winsdata);
        %         plot(timel,data,timel(indm),data(indm),'g',timel(min(indm)),data(min(indm)),'*r')
        %         title('data')
        %         xlim([1.04 1.13])
        %     end
        %         
        % end
        % t_hpointd = hpoint_std/10;
        % t_hpointd(:,6) = mean(t_hpointd(:,1:5),2);
        % t_hpointd(:,7) = std(t_hpointd(:,1:5),0,2);
        % figure('units','normalized','outerposition',[0 0 .9 .9]);
        % errorbar([0 20 40 60 80 100],t_hpointd(:,6),t_hpointd(:,7),'ko')
        % xlabel('dose')
        % grid on
        % ylabel('time_hit_data (sec)')
        % ylim([1040 1100])