% Felipe A. Mendez, INB, UNAM, Qro. México, 2020

% Continuation of neonatal hippocampal slice recordings
% Pre-processed with wavelet decomposition and reconstruction

% Peak detection and extraction: fGDP, bursts and spikes
% fGDPs must be accompanied by a spike burst to be valid

% files to use: clean_mat
% clean_mat(:,1) raw signal
% clean_mat(:,2) fGDP
% clean_mat(:,3) spikes
clear all, close all, clc

load('KO_P7_210920_tidy.mat');

v_raw = clean_mat(:,1);
fGDP = clean_mat(:,2);
spikes = clean_mat(:,3);

Fs= 10000;
% Burst detection; an approximation of the analog integrator 

window = 501; % samples = 50ms
smth_spks = smoothdata(abs(spikes), 'gaussian', window);
%%
% check it out
% figure(1)                   %pretty plot to show off
% plot(spikes), hold on
% plot(smth_spks, 'LineWidth', 1.5), hold on
% plot(fGDP, 'LineWidth', 2)
%% Peak detection
perd= .5; % minimum distance in seconds
% find burst peaks
brst_th = 1.4826*median(smth_spks)+10*mad(smth_spks);
[brst_pks, brst_locs, brst_w] = findpeaks(smth_spks, Fs,...
    'MinPeakHeight', brst_th, 'MinPeakDistance', perd);
% find GDP peaks
fGDP_th = 1.4826*median(-fGDP)+10*mad(-fGDP);
[gdp_pks, gdp_locs, gdp_w] = findpeaks(-fGDP, Fs, ...
    'MinPeakHeight', fGDP_th, 'MinPeakDistance', perd, 'MaxPeakWidth', 0.4);
%%
% Testify!!
% figure(2)
% findpeaks(-fGDP, Fs, 'MinPeakHeight', fGDP_th, 'MinPeakDistance', perd)
% hold on
% findpeaks(smth_spks, Fs, 'MinPeakHeight', brst_th, 'MinPeakDistance', perd)
%% Determine Coincidences
% Not very efficient but does the trick
m_idx = 1;
for gp = 1:length(gdp_locs)
    wk_gdp = gdp_locs(gp);
    for sb = 1:length(brst_locs)
        wk_brst = brst_locs(sb);
        jitter = abs(wk_gdp-wk_brst);
        if jitter <= .1
            match(m_idx) = gp;
            m_idx = m_idx+1;
            break
        end
    end
end

%% GDP Intervals
% Clean variables pks, locs, w

gdp_pks = gdp_pks(match); % absolutes values of the valley
gdp_locs = gdp_locs(match); % still in seconds
gdp_w = gdp_w(match)*1000; % transform to ms

% calculate intervals
gdp_int = NaN(length(gdp_locs),1);
for jj = 2:length(gdp_locs)
int= gdp_locs(jj)-gdp_locs(jj-1);
gdp_int(jj) = int;
end
%%
% figure(3)
% subplot(131)
% histogram(gdp_int, 'BinMethod', 'fd')
% title('GDP Interval')
% ylabel('Count')
% xlabel('Seconds')
% subplot(132)
% histogram(gdp_pks, 'BinMethod', 'fd')
% title('GDP Amplitude')
% ylabel('Count')
% xlabel('\muV')
% subplot(133)
% histogram(gdp_w, 'BinMethod', 'fd')
% title('GDP half-width')
% ylabel('Count')
% xlabel('Milliseconds')
%%
clearvars jj m_idx gp sb wk_brst wk_gdp int jitter
%% Spikes extraction
ref_pd = 0.002;
spks_th = 1.4826*median(spikes)+5*mad(spikes);
[spks_pks, spks_locs] = findpeaks(-spikes, Fs, 'MinPeakHeight', spks_th,...
    'MinPeakDistance', ref_pd, 'MinPeakWidth', 0.001);
%% Spike Interval 

% calculate intervals
spks_int = NaN(length(spks_locs),1);
for jj = 2:length(spks_locs)
int= spks_locs(jj)-spks_locs(jj-1);
spks_int(jj) = int;
end

clearvars jj int
%% Find and validate the positive thrusts in fGDP
gdp_ps_th = 1.4826*median(fGDP)+2*mad(fGDP);
[pos_pks, pos_locs] = findpeaks(fGDP, Fs, 'MinPeakHeight', gdp_ps_th, 'MinPeakDistance', perd);

m_idx = 1;
for gp = 1:length(gdp_locs)
    wk_gdp = gdp_locs(gp);
    for sb = 1:length(pos_locs)
        wk_brst = pos_locs(sb);
        jitter = abs(wk_brst-wk_gdp);
        if jitter <= .5
            pos_match(m_idx) = sb;
            m_idx = m_idx+1;
            break
        end
    end
end
pos_locs = pos_locs(pos_match);
clearvars m_indx gp sb wk_brst jitter
%% Find the spikes inside fGDPs
intra_gdp_isi = [];
for mm = 1:length(gdp_locs)
    spikes_add = spks_locs(spks_locs>(gdp_locs(mm)-gdp_w(mm)/1000) & spks_locs<= pos_locs(mm));
    isis = [];
    for ll = 2:length(spikes_add) % And get isi
        isi = spikes_add(ll)-spikes_add(ll-1);
        isis = horzcat(isis, isi);
    end
    intra_gdp_isi = horzcat(intra_gdp_isi, isis);
end
%% Get the frequency of all non GDPs spikes
inter_gdp_isi = [];
for mm = 1:length(gdp_locs)
    if mm == 1
        spikes_add = spks_locs(spks_locs<(gdp_locs(mm)-gdp_w(mm)/1000));
    elseif mm == length(gdp_locs)
        spikes_add = spks_locs(spks_locs> pos_locs(mm));
    else
        spikes_add = spks_locs(spks_locs> pos_locs(mm) & spks_locs<(gdp_locs(mm+1)-gdp_w(mm+1)/1000));
    end
    isis = [];
    for ll = 2:length(spikes_add) % And get isi
        isi = spikes_add(ll)-spikes_add(ll-1);
        isis = horzcat(isis, isi);
    end
    inter_gdp_isi = horzcat(inter_gdp_isi, isis);
end
%% Metrics Histogram

% figure(4)
% subplot(131)
% histogram(spks_int, 'BinMethod', 'fd')
% title('All Spikes Interval')
% xlabel('Seconds')
% subplot(132)
% histogram(inter_gdp_isi, 'BinMethod', 'fd')
% title('Inter-fGDPs Spike Interval')
% xlabel('Seconds')
% subplot(133)
% histogram(intra_gdp_isi, 'BinMethod', 'fd')
% title('Intra-fGDPs Spike Interval')
% xlabel('Seconds')
%% Extract Centrality Metrics for all Variables of Interest
% For fGDPs: Amplitude, half-width and interval
mean_GDP_amp = mean(gdp_pks);
mean_GDP_width = mean(gdp_w);
mean_GDP_int = mean(gdp_int, 'omitnan');
% For MUAs
mean_allspikes_fq = 1/mean(spks_int, 'omitnan');
mean_interspikes_fq = 1/mean(inter_gdp_isi);
mean_intraspikes_fq = 1/mean(intra_gdp_isi);
% Error estimate
spike_split_accuracy = 100*(((length(inter_gdp_isi)+length(intra_gdp_isi)))/length(spks_locs));
gdp_neg_pos_accuracy = 100*(length(pos_locs)/length(gdp_locs));
% And tabulate
metrics_tab = table(mean_GDP_amp, mean_GDP_width, mean_GDP_int...
    ,mean_allspikes_fq, mean_interspikes_fq, mean_intraspikes_fq, spike_split_accuracy, gdp_neg_pos_accuracy)
%%
spk_rate = firingrate(spks_locs, 0:0.1:spks_locs(end), 0.2);
gdp_rate = firingrate(gdp_locs, 0:0.1:gdp_locs(end), 0.2);

[xcor, lags] = xcorr(spk_rate(1:length(gdp_rate)), gdp_rate, 'coeff');

[peak_cor, lag_idx] = max(abs(xcor));
peak_cor
lag = lags/10;

figure(5)
plot(lag, xcor)
cor_seg = xcor(lag_idx-10:lag_idx+10);
lag_seg = lags(lag_idx-10:lag_idx+10)/10;
