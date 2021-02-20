Ca_rest = 50e-9;
T_Ca_decay = 48;
Ca_residual = 250e-9;
Ca_spike = 22.5e-6;
FWHM = .5;
mu = 2*FWHM;
sigma = FWHM/2.35;
delta_t = 0.1;

%make 1hz
stimulus_times = linspace(0,1000*19,20);
max_time = stimulus_times(end) + stimulus_times(2)*10;

ts = linspace(0,max_time,max_time/delta_t + 1);
Ca_sim = zeros(1,length(ts));
Ca_sim = Ca_sim + Ca_rest;
   
for i = 1
    
    stimulus_times_2 = stimulus_times;
    for t = 1:length(stimulus_times_2) %simulate calcium influx

        spike_start_index = round(stimulus_times_2(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
        spike_peak_index = round((stimulus_times_2(t)+mu)/delta_t) + 1;
        
        Ca_sim(i,spike_peak_index:end) = (Ca_sim(i,spike_peak_index) + Ca_residual - Ca_rest)*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay)+Ca_rest; 
        
    end
%     for t = 1:length(stimulus_times_2) %simulate calcium influx
% 
%         spike_start_index = round(stimulus_times_2(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
%         spike_peak_index = round((stimulus_times_2(t)+mu)/delta_t) + 1;
% 
%         Ca_sim(i,spike_start_index:end) = Ca_sim(i,spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one
%         
%     end
end
semilogy(ts,Ca_sim(round(ts/.1)+1));
save('CF1HzCa.mat','Ca_sim');

%make 10hz
stimulus_times = linspace(0,100*19,20);
max_time = stimulus_times(end) + stimulus_times(2)*10;

ts = linspace(0,max_time,max_time/delta_t + 1);
Ca_sim = zeros(1,length(ts));
Ca_sim = Ca_sim + Ca_rest;
   
for i = 1
    
    stimulus_times_2 = stimulus_times;
    for t = 1:length(stimulus_times_2) %simulate calcium influx

        spike_start_index = round(stimulus_times_2(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
        spike_peak_index = round((stimulus_times_2(t)+mu)/delta_t) + 1;
        
        Ca_sim(i,spike_peak_index:end) = (Ca_sim(i,spike_peak_index) + Ca_residual - Ca_rest)*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay)+Ca_rest; 
        
    end
%     for t = 1:length(stimulus_times_2) %simulate calcium influx
% 
%         spike_start_index = round(stimulus_times_2(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
%         spike_peak_index = round((stimulus_times_2(t)+mu)/delta_t) + 1;
% 
%         Ca_sim(i,spike_start_index:end) = Ca_sim(i,spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one
%         
%     end
end
semilogy(ts,Ca_sim(round(ts/.1)+1));
save('CF10HzCa.mat','Ca_sim');

%make 20hz
stimulus_times = linspace(0,50*19,20);
max_time = stimulus_times(end) + stimulus_times(2)*10;

ts = linspace(0,max_time,max_time/delta_t + 1);
Ca_sim = zeros(1,length(ts));
Ca_sim = Ca_sim + Ca_rest;
   
for i = 1
    
    stimulus_times_2 = stimulus_times;
    for t = 1:length(stimulus_times_2) %simulate calcium influx

        spike_start_index = round(stimulus_times_2(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
        spike_peak_index = round((stimulus_times_2(t)+mu)/delta_t) + 1;
        
        Ca_sim(i,spike_peak_index:end) = (Ca_sim(i,spike_peak_index) + Ca_residual - Ca_rest)*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay)+Ca_rest; 
        
    end
%     for t = 1:length(stimulus_times_2) %simulate calcium influx
% 
%         spike_start_index = round(stimulus_times_2(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
%         spike_peak_index = round((stimulus_times_2(t)+mu)/delta_t) + 1;
% 
%         Ca_sim(i,spike_start_index:end) = Ca_sim(i,spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one
%         
%     end
end
semilogy(ts,Ca_sim(round(ts/.1)+1));
save('CF20HzCa.mat','Ca_sim');

%make 50hz w/50hz rec
stimulus_times = linspace(0,20*19,20);
Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
max_time = stimulus_times(end) + stimulus_times(2)*10 + Recovery(end); 

ts = linspace(0,max_time,max_time/delta_t + 1);
Ca_sim = zeros(length(Recovery),length(ts));
Ca_sim = Ca_sim + Ca_rest;
   
for i = 1:length(Recovery)
    
    stimulus_times_2 = [stimulus_times stimulus_times(end)+Recovery(i)];
 
    for t = 1:length(stimulus_times_2) %simulate calcium influx

        spike_start_index = round(stimulus_times_2(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
        spike_peak_index = round((stimulus_times_2(t)+mu)/delta_t) + 1;
        
        Ca_sim(i,spike_peak_index:end) = (Ca_sim(i,spike_peak_index) + Ca_residual - Ca_rest)*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay)+Ca_rest; 
        
    end
%     for t = 1:length(stimulus_times_2) %simulate calcium influx
% 
%         spike_start_index = round(stimulus_times_2(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
%         spike_peak_index = round((stimulus_times_2(t)+mu)/delta_t) + 1;
% 
%         Ca_sim(i,spike_start_index:end) = Ca_sim(i,spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one
%         
%     end
end
semilogy(ts,Ca_sim(round(ts/.1)+1));
save('CF50HzCa.mat','Ca_sim');
