
%% Explore the file
clear all, close all
file = '/Volumes/KINGSTON/WT_P7_111220.abf';
preview_rec = abfload(file);
v_sig = preview_rec(:,1);
plot(v_sig)
length(v_sig)/600000 % Recording duration in minutes
%%
% Wavelet decomposition to extract peaks
%   THIS IS JUST FOR CLEANING THE RECORDINGS AND GET THE WAVELET 
%   DECOMPOSITION AND RECONSTRUCTION
% Since artifacts are heterogeneous, cleaning will have to be semi-manual
% load abf 
clear all, close all
file = '/Volumes/KINGSTON/WT_P7_111220.abf';

Fs = 10000; % sample rate
start = [334, 1801]; %first 30min
stop = [1801, 3601]; %last 30min

clean_rec = [];
fGDP = [];
spikes = [];
for half = 1:2

rec = abfload(file, 'start', start(half), 'stop', stop(half));
% get the first hr
disp('Recording Loaded');
v_sig = rec(:,1);% channel 1 contains V signal
v_sig = v_sig-median(v_sig); % correct offset
disp('Zeroed');

disp('Working in Wavelet decomposition');
dwt = modwt(v_sig, 'sym3', 10); %perform a discrete wavelet transform
disp('Done')
disp('Reconstruction');
recons= modwtmra(dwt, 'sym3'); % reconstruct the signal by wv coefficients
disp('Done');

clean_rec = vertcat(clean_rec, v_sig);
fGDP = horzcat(fGDP, sum(recons(7:11, :))); % low components fGDPs
spikes = horzcat(spikes, sum(recons(3:6, :))); % de-noised spikes

 
% figure()
% subplot(3,1,1)
% plot(v_sig)                 %plot original
% subplot(3,1,2)
% plot(sum(recons(3:8, :)))   %plot denoised high frequency components (spikes)
% subplot(3,1,3)
% plot(sum(recons(6:11, :)))  %plot low frequency components (fGDPs)
% shg

end
clean_mat = [clean_rec fGDP' spikes'];

figure()                   %pretty plot to show off
plot(clean_mat(:,1)), hold on
plot(clean_mat(:,3), 'LineWidth', 1.5), hold on
plot(clean_mat(:,2), 'LineWidth', 2)
%%
% Run if only 30min recodring is needed or useful
clear all, close all
file = '/Volumes/Seagate Backup Plus Drive/GDPs/WT_P9_300620_1.abf';
Fs = 10000; % sample rate
start = 1; %first 30min
stop = 1801;

rec = abfload(file, 'start', start, 'stop', stop);
% get the first hr
disp('Recording Loaded');
v_sig = rec(:,1);% channel 1 contains V signal
v_sig = v_sig-median(v_sig); % correct offset
disp('Zeroed');

disp('Working in Wavelet decomposition');
dwt = modwt(v_sig, 'sym3', 10); %perform a discrete wavelet transform
disp('Done')
disp('Reconstruction');
recons= modwtmra(dwt, 'sym3'); % reconstruct the signal by wv coefficients
disp('Done');

clean_rec = v_sig;
fGDP = sum(recons(7:11, :)); % low components fGDPs
spikes = sum(recons(3:6, :)); 

clean_mat = [clean_rec fGDP' spikes'];

figure()                   %pretty plot to show off
plot(clean_mat(:,1)), hold on
plot(clean_mat(:,3), 'LineWidth', 1.5), hold on
plot(clean_mat(:,2), 'LineWidth', 2)