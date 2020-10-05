function plot_one(data, stimulus_times)

    max_time = stimulus_times(end) + stimulus_times(2)*3;
    [ts, state, Fused_im, Fused_m, Ca_sim] = stim_sim(stimulus_times, max_time);
    Fused = Fused_im + Fused_m;
    Fused = Fused/Fused(1);
    
    figure
    subplot(3,2,1)
    plot(Fused,'ko','MarkerSize', 5)
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 length(Fused)])
    set(gca,'ylim',[0 1])
    hold on
    
    hz_1_data = data(:,1);
    hz_10_data = data(:,2);
    hz_20_data = data(:,3);
    hz_50_data = data(:,4);
    hz_100_data = data(:,5);
    hz_200_data = data(:,6);

    labels = ["1 Hz data", "10 Hz data", "20 Hz data", "50 Hz data", "100 Hz data", "200 Hz data"];
    

    plot(hz_1_data,'b')
    plot(hz_10_data,'Color',[17 17 17]/255)
    plot(hz_20_data,'g')
    plot(hz_50_data,'m')
    plot(hz_100_data,'r')
    plot(hz_200_data,'y')

    legend('Simulation','1 Hz data','10 Hz data','20 Hz data','50 Hz data','100 Hz data','200 Hz data')
    hold off
    
    subplot(3,2,2)
    semilogy(linspace(1,max_time,length(Ca_sim)),Ca_sim)
    ylabel('Simulated Calcium Concentration')
    xlabel('Time (ms)')
    set(gca,'xlim',[ts(1) max_time])

    subplot(3,2,3)
    hold on

    plot(ts, state(:,2),'-.b')
    plot(ts, state(:,3),'-b')
    plot(ts, state(:,2)+state(:,3),'-r')
    xlabel('Time (ms)')
    ylabel('Occupied Sites (norm)')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])

    hold off
    legend('Occupied Low P Sites','Occupied High P Sites','Total Occupied Sites');
    
    subplot(3,2,4)
    hold on

    plot(ts, state(:,1),'-r')
    xlabel('Time (ms)')
    ylabel('Empty Sites (norm)')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])

    hold off

    subplot(3,2,5)
    plot(ts, state(:,4))
    xlabel('Time (ms)')
    ylabel('Syt3')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])

    subplot(3,2,6)
    plot(ts, state(:,5))
    xlabel('Time (ms)')
    ylabel('Syt7')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])
    
end