function plot_six(data,Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50, Fused_im_100, Fused_m_100, hz_100, Fused_im_200, Fused_m_200, hz_200)
    
    hz_1_data = data(:,1);
    hz_10_data = data(:,2);
    hz_20_data = data(:,3);
    hz_50_data = data(:,4);
    hz_100_data = data(:,5);
    hz_200_data = data(:,6);

    labels = ["1 Hz data", "10 Hz data", "20 Hz data", "50 Hz data", "100 Hz data", "200 Hz data"];
    
    Fused_1_norm = Fused_im_1(1) + Fused_m_1(1);
    Fused_10_norm = Fused_im_10(1) + Fused_m_10(1);
    Fused_20_norm = Fused_im_20(1) + Fused_m_20(1);
    Fused_50_norm = Fused_im_50(1) + Fused_m_50(1);
    Fused_100_norm = Fused_im_100(1) + Fused_m_100(1);
    Fused_200_norm = Fused_im_200(1) + Fused_m_200(1);
    
    figure 
    subplot(3,2,1)
    plot(hz_1_data,'-k')
    title('1 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_1,'ko','Markersize',5)
    plot(Fused_im_1/Fused_1_norm,'rv')
    plot(Fused_m_1/Fused_1_norm,'g^')
    legend(labels(1),'1 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(3,2,2)
    plot(hz_10_data,'-k')
    title('10 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_10,'ko','Markersize',5)
    plot(Fused_im_10/Fused_10_norm,'rv')
    plot(Fused_m_10/Fused_10_norm,'g^')
    %legend(labels(2),'10 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(3,2,3)
    plot(hz_20_data,'-k')
    title('20 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_20,'ko','Markersize',5)
    plot(Fused_im_20/Fused_20_norm,'rv')
    plot(Fused_m_20/Fused_20_norm,'g^')
    %legend(labels(3),'20 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(3,2,4)
    plot(hz_50_data,'-k')
    title('50 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_50,'ko','Markersize',5)
    plot(Fused_im_50/Fused_50_norm,'rv')
    plot(Fused_m_50/Fused_50_norm,'g^')
    %legend(labels(4),'50 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(3,2,5)
    plot(hz_100_data,'-k')
    title('100 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_100,'ko','Markersize',5)
    plot(Fused_im_100/Fused_100_norm,'rv')
    plot(Fused_m_100/Fused_100_norm,'g^')
    
    subplot(3,2,6)
    plot(hz_200_data,'-k')
    title('200 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_200,'ko','Markersize',5)
    plot(Fused_im_200/Fused_200_norm,'rv')
    plot(Fused_m_200/Fused_200_norm,'g^')
end