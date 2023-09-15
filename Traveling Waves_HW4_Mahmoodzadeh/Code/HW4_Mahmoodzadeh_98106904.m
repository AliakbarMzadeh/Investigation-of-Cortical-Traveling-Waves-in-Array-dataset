%% part.1 - load 
clc; clear; close all;
ArrayData = load('ArrayData.mat');
CleanTrials = load('CleanTrials.mat');

%% -  trials
clc; close all;

LFPsignals = zeros(length(ArrayData.chan),size(ArrayData.chan(1).lfp,1),...
    length(CleanTrials.Intersect_Clean_Trials));

for i = 1:size(LFPsignals,1)
    LFPsignals(i,:,:) = ArrayData.chan(i).lfp(:,CleanTrials.Intersect_Clean_Trials);
end

%% part.1.a - removing pink 
clc; close all;




Fs = 1/(ArrayData.Time(2)-ArrayData.Time(1));
DominantFreqs = zeros(1,size(LFPsignals,1));
powAllelecs = zeros(size(LFPsignals,1),size(LFPsignals,3),320);
powAllelecsPinkFree = zeros(size(LFPsignals,1),size(LFPsignals,3),320);
pinkNoiseAll = zeros(size(LFPsignals,1),size(LFPsignals,3),320);

for i = 1:size(LFPsignals,1)
    for j = 1:size(LFPsignals,3)
        targetTrial = j;
        targetElectrode = i;
        targetSig = (LFPsignals(targetElectrode,:,j));

        % signal fft
        [fSignal, powSignal] = calFFT(targetSig,Fs);
        powAllelecs(i,j,:) = powSignal;

        % pink noise 
        a = -1;
        b = mean(log10(powSignal) - a*log10((fSignal)));
        pinkNoiseAll(i,j,:) = a*log10(fSignal) + b;
        
        % remove pink noise
        dePinkedSignal = (log10(powSignal)-...
            reshape(pinkNoiseAll(i,j,:),[1 320]));
        
        pinkNoiseAll(i,j,:) = 10.^pinkNoiseAll(i,j,:);
        powAllelecsPinkFree(i,j,:) = 10.^dePinkedSignal;
    end
    
    %%% plot raw signal and pink noise - plot 48 figures! be careful.
    plotRawSigPinkNoise(log10(fSignal),...
    reshape(log10(mean(pinkNoiseAll(i,:,:),2)),[1 320]),...
    reshape(log10(mean(powAllelecs(i,:,:),2)),[1 320]),...
    reshape(log10(mean(powAllelecsPinkFree(i,:,:),2)),[1 320]),...
    targetElectrode);

    % find the dominant frequency
    finded = (find(mean(powAllelecsPinkFree(i,:,:),2) == ...
        max(mean(powAllelecsPinkFree(i,:,:),2))));
    DominantFreqs(i) = fSignal(finded);
    
end

% Add this code after your for loop
figure;

% Plot all raw data and their average
subplot(1, 2, 1);
hold on;
for i = 1:size(LFPsignals,1)
    plot(10.^fSignal, 10.^squeeze(log10(mean(powAllelecs(i,:,:),2))), 'LineWidth', 1);
end
plot(10.^fSignal, 10.^squeeze(log10(mean(mean(powAllelecs,1),2))), 'k', 'LineWidth', 2); % Raw data average in bold
hold off;
title('All Raw Data and their Average', 'interpreter', 'latex');
xlabel('Frequency(Hz)', 'interpreter', 'latex');
ylabel('Power', 'interpreter', 'latex');
grid on; grid minor;

% Plot all pink noise-free data and their average
subplot(1, 2, 2);
hold on;
for i = 1:size(LFPsignals,1)
    plot(10.^fSignal, 10.^squeeze(log10(mean(powAllelecsPinkFree(i,:,:),2))), 'LineWidth', 1);
end
plot(10.^fSignal, 10.^squeeze(log10(mean(mean(powAllelecsPinkFree,1),2))), 'k', 'LineWidth', 2); % Pink noise-free average in bold
hold off;
title('All Pink Noise-Free Data and their Average', 'interpreter', 'latex');
xlabel('Frequency(Hz)', 'interpreter', 'latex');
ylabel('Power', 'interpreter', 'latex');
grid on; grid minor;




%%
clc
clear
close all
load("ArrayData.mat")
load("CleanTrials.mat")
fs = 200;

num_channels = numel(chan);
% removing bad trials
for ch_no = 1:num_channels
    chan(ch_no).lfp = chan(ch_no).lfp(:, Intersect_Clean_Trials);
end
num_trials = size(chan(1).lfp, 2);


Ps = 0;
for ch_no = 1:num_channels
    lfp_data = chan(ch_no).lfp;
    lfp_data = zscore(lfp_data);
    for trial_no = 1:num_trials
        trial_data = lfp_data(:, trial_no);
        m = length(trial_data);
        n = pow2(nextpow2(m));
        Y = fft(trial_data, n);
        Y = fftshift(Y);
        Ps = Ps+abs(Y);
    end
end
normalize_constant = 10*log10((num_channels*num_trials)^2/n);

f = (-n/2:n/2-1)*(fs/n);
Ps_plot = 10*log10(Ps.^2/n);
pink_noise = 1./f(n/2+2:end);
[~,~,spectrum_regressed] = regress(Ps_plot(n/2+2:end), pink_noise');
pink_spectrum = Ps_plot(n/2+2:end) - spectrum_regressed;

% part b - clustering electrodes
clc
close all

figure
hold on
dominant_freq_mat = ChannelPosition*nan;
normalize_constant = 10*log10((num_trials)^2);
for ch_no = 1:num_channels
    lfp_data = chan(ch_no).lfp;
    Ps = 0;
    for trial_no = 1:num_trials
        trial_data = lfp_data(:, trial_no);
        trial_data = zscore(trial_data);
        m = length(trial_data);
        n = pow2(nextpow2(m));
        Y = fft(trial_data, n);
        Y = fftshift(Y);
        Ps = Ps+abs(Y);
    end
    f = (-n/2:n/2-1)*(fs/n);
    Ps = 10*log10(Ps.^2/n);
    Ps_plot = removePinkNoise(Ps, f, n, 1);
    plot(f(n/2+2:end), Ps_plot(n/2+2:end)-normalize_constant)
    [row, ~] = find(Ps_plot(n/2+2:end) == max(Ps_plot(n/2+2:end)));
    f_tmp = f(n/2+2:end);
    dominant_freq = f_tmp(row);
    [row, col] = find(ChannelPosition==ch_no);
    dominant_freq_mat(row, col) = dominant_freq;
end





figure
plt = imagesc(dominant_freq_mat);
set(plt,'AlphaData', ~isnan(dominant_freq_mat))

% Add black color for NaN values in the colormap
custom_colormap = [0, 0, 0;linspace(1, 0, 10)',linspace(0, 1, 10)', linspace(0, 1, 10)'];
colormap(custom_colormap)

c_bar = colorbar;
% Set the color axis
caxis([0, 10])  % Including the black color for NaN values
title('Dominant Frequencies','interpreter','latex')
ylabel(c_bar,'Frequency (Hz)','interpreter','latex')

% Display frequency values in each cell
textStrings = num2str(dominant_freq_mat(:), '%0.2f');
textStrings = strtrim(cellstr(textStrings));  
[x, y] = meshgrid(1:size(dominant_freq_mat,2), 1:size(dominant_freq_mat,1));  
hStrings = text(x(:), y(:), textStrings(:), ...
                'HorizontalAlignment', 'center', 'Color', 'k');

% Set y-axis labels
yticks(1:5)
yticklabels({'1', '2', '3', '4', '5'})

            
% Create bar graph for the frequency ranges
freq_ranges = [0, 10, 10, 12, 12.2, 12.5, 12.5, 12.7, 12.7, Inf];
freq_labels = {'<10', '10-12', '12.2-12.5', '12.5-12.7', '>12.7'};
counts = zeros(1, length(freq_labels));

for i = 1:2:length(freq_ranges)
    range_idx = ceil(i / 2);
    counts(range_idx) = sum(dominant_freq_mat(:) > freq_ranges(i) & dominant_freq_mat(:) <= freq_ranges(i+1));
end

figure
bar(counts)
xticklabels(freq_labels)
xlabel('Frequency Ranges (Hz)','interpreter','latex')
ylabel('Number of Electrodes','interpreter','latex')
title('Electrode/Pixel Distribution in Dominant Frequency Plot','interpreter','latex')
grid on; grid minor;
if save_figures
    set(gcf,'PaperPositionMode','auto')
    print("Report/images/b_2",'-dpng','-r0')
    
    
  
end





close all


% Create a custom colormap
n_colors = 256;
key_colors = [0, 0, 1; 0, 1, 0; 1, 1, 0; 1, 0, 0];
x_key = linspace(1, n_colors, size(key_colors, 1));
x_query = 1:n_colors;
custom_colormap = interp1(x_key, key_colors, x_query, 'linear');



% Welch Method
window_length = 80;
num_overlap_samples = 60;
pxx_mean = 0;
for ch_no = 1:num_channels
    lfp_data = chan(ch_no).lfp;
    lfp_data = zscore(lfp_data);
    for trial_no = 1:num_trials
        trial_data = lfp_data(:, trial_no);
        tmp = buffer(trial_data, window_length, num_overlap_samples);
        [pxx,f] = pwelch(tmp, 40, 20, fs, fs);
        pxx_mean = pxx_mean + pxx;
    end
end

pxx_clean = [];
for t = 1:size(pxx_mean, 2)
    Ps = pxx_mean(:, t);
    n = length(Ps);
    Ps_plot = removePinkNoise(Ps, f', n, 2);
    pxx_clean(:, t) = Ps_plot;
end

figure
imagesc(linspace(-1.2, 2, size(tmp, 2)),f, pxx_mean/(num_channels*num_trials))
c_bar = colorbar;
set(gca,'YDir','normal')
title("Power Spectrum  Welch",'interpreter','latex')
xlabel('Time (s)','interpreter','latex')
ylabel('Frequency (Hz)','interpreter','latex')
ylabel(c_bar,'Power','interpreter','latex')
ylim([0, 40])
xline(0,'k','LineWidth',5);
% Set colormap to custom colormap
colormap(gca, custom_colormap)

% Update colorbar tick labels
c_ticks = linspace(-0.05, 1, 6); % 6 ticks from -0.05 to 1
c_tick_labels = round(linspace(1, 70, 6)); % Corresponding tick labels from 1 to 70
c_bar.Ticks = c_ticks;
c_bar.TickLabels = c_tick_labels;



figure
imagesc(linspace(-1.2, 2, size(tmp, 2)),f, pxx_clean/(num_channels*num_trials))
c_bar = colorbar;
set(gca,'YDir','normal')
title("Power Spectrum Welch - Without Noise",'interpreter','latex')
xlabel('Time (s)','interpreter','latex')
ylabel('Frequency (Hz)','interpreter','latex')
ylabel(c_bar,'Power','interpreter','latex')
ylim([0, 40])
xline(0,'k','LineWidth',5);
% Set colormap to custom colormap
colormap(gca, custom_colormap)

% Update colorbar tick labels
c_bar.Ticks = c_ticks;
c_bar.TickLabels = c_tick_labels;


%% part c - WELCH
clc
close all






% Welch Method
window_length = 80;
num_overlap_samples = 60;
pxx_mean = 0;
for ch_no = 1:num_channels
    lfp_data = chan(ch_no).lfp;
    lfp_data = zscore(lfp_data);
    for trial_no = 1:num_trials
        trial_data = lfp_data(:, trial_no);
        tmp = buffer(trial_data, window_length, num_overlap_samples);
        [pxx,f] = pwelch(tmp, 40, 20, fs, fs);
        pxx_mean = pxx_mean + pxx;
    end
end

pxx_clean = [];
for t = 1:size(pxx_mean, 2)
    Ps = pxx_mean(:, t);
    n = length(Ps);
    Ps_plot = removePinkNoise(Ps, f', n, 2);
    pxx_clean(:, t) = Ps_plot;
end

figure
imagesc(linspace(-1.2, 2, size(tmp, 2)),f, pxx_mean/(num_channels*num_trials))
c_bar = colorbar;
set(gca,'YDir','normal')
title("Power Spectrum over Time - Welch")
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylabel(c_bar,'Power')
caxis([-0.05, 0.1])
ylim([0, 40])



figure
imagesc(linspace(-1.2, 2, size(tmp, 2)),f, pxx_clean/(num_channels*num_trials))
c_bar = colorbar;
set(gca,'YDir','normal')
title("Power Spectrum over Time - Welch - Denoised (no Pink noise)")
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylabel(c_bar,'Power')
caxis([-0.05, 0.1])
ylim([0, 40])



%% part 2.a - bandpass around dominant frequency
clc; close all;



% design a second order butterworth
fc = 12.5;
w_filter = 1.5;
Fs = 1/(ArrayData.Time(2)-ArrayData.Time(1));
order = 2;
[b,a] = butter(order, [fc-w_filter, fc+w_filter]/Fs*2,'bandpass');

% as you can see, the phase of the designed filter is linear in our
% frequency range of interest ([12 13])
figure;
freqz(b,a,1000,Fs);
grid on; 

%% Filter




%% Load Data_ 
clc
clearvars

load('ArrayData.mat');
load('CleanTrials.mat');

Fs = 200;
data_PSD = [];

for i=1:length(chan)
    x = chan(i).lfp;
    x = x(:,Intersect_Clean_Trials);
    for j=1:size(Intersect_Clean_Trials,1)

        N = length(x(:,j));
        xdft = fft(x(:,j));
        xdft = xdft(1:floor(N/2)+1);
        psdx = (1/(Fs*N)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:Fs/length(x(:,j)):Fs/2;
        data_PSD(i,j,:) = psdx;
    end
end
data_mean_PSD = squeeze(mean(data_PSD,2));


dominant_freq = ChannelPosition;
Coef = [];
for i=1:size(data_mean_PSD,1)
    i

    x = log10(freq);
    y = log10(data_mean_PSD(i,:));
    c = polyfit(x(2:end),y(2:end),1);
    Coef(i,:) = c;

    y = y - min(y(2:end));
    z = log10(data_mean_PSD(i,:))-polyval(c,log10(freq));
    [ma, mai] = max(z-min(z(2:end))); 

    dominant_freq(ChannelPosition == i) = freq(mai);

end
save('Dominant_freqs.mat','dominant_freq','Coef');


%%%%%%




clc
clearvars

load('ArrayData.mat');
load('CleanTrials.mat');
load('Dominant_freqs.mat');

fd = mean(mean([dominant_freq(1,2:end) dominant_freq(2,2:end) dominant_freq(3,:) dominant_freq(4,:) dominant_freq(5,2:end)]));

Fs = 200;
[b,a] = butter(2,[fd-0.5 fd+0.5]./(Fs/2));
figure;
set(gcf,'Color',[1 1 1]);

freqz(b,a,[],Fs)
title('Filter Specification');
export_fig('Filter.png','-r600');

num_windows = 12;
channels_per_window = 4;

for w = 1:num_windows
    fi = figure;
    set(gcf,'Color',[1 1 1]);
    
    for p = 1:channels_per_window
        i = (w-1) * channels_per_window + p;
        
        if i > length(chan)
            break;
        end
        
        x = chan(i).lfp;
        x = x(:,Intersect_Clean_Trials);
        x = x(:,10);
        y = filter(b,a,x);
        
        subplot(channels_per_window, 1, p);
        plot(Time, x, 'r', 'LineWidth', 2);
        hold on
        plot(Time, y, 'k', 'LineWidth', 2);
        legend({'RAW Data' 'Filtered(Butter2th) Data'},'interpreter','latex')
        xlim([Time(1),Time(end)]);
        xlabel('Time (s)','interpreter','latex');
        ylabel(['Ch' num2str(i) ' Mag (V)'],'interpreter','latex');
        grid on; grid minor;
    end
    export_fig(['result/Filter_W' num2str(w)],'-r600');
    close(fi)
end
close all


% Phase
clc
clearvars

load('ArrayData.mat');
load('CleanTrials.mat');
load('Dominant_freqs.mat');

fd = mean(mean([dominant_freq(1,2:end) dominant_freq(2,2:end) dominant_freq(3,:) dominant_freq(4,:) dominant_freq(5,2:end)]));

Fs = 200;
[b,a] = butter(2,[fd-0.5 fd+0.5]./(Fs/2));
sigphase_data = [];
for i=1:length(chan)
    i
    x = chan(i).lfp;
    x = x(:,Intersect_Clean_Trials);
    sigphase = [];
    for j=1:size(x,2)
        xs = x(:,j);
        y = filter(b,a,xs);
        y = hilbert(y);
        sigphase = [sigphase unwrap(angle(y))];
    end
    sigphase_data(i,:,:) = sigphase;
end



cmap = [linspace(0, 1, 128)', linspace(0, 1, 128)', linspace(1, 1, 128)';
        linspace(1, 1, 128)', linspace(1, 0, 128)', linspace(1, 0, 128)'];

for i=1:size(sigphase_data,3)
for i=[10]
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
        fi = figure;
        set(gcf,'Color',[1 1 1]);
        set(gca,'FontName','arial','FontSize',10);
        heatmap(phase_map);
        colormap(gca, cmap);
        caxis([-1, 1]);
        export_fig(['result/S' num2str(j) '.png'],'-r600');
        close(fi);
    end
end
close all

writerObj = VideoWriter('TravellingWave.avi');
writerObj.FrameRate = 4;
open(writerObj);
for j=1:size(sigphase_data,2)
    j
    img =  imread(['result/S' num2str(j) '.png']);
    frame = im2frame(img);
    writeVideo(writerObj, frame);
end
close(writerObj);

%% PGD
%% PGD
pgd_data = [];
for i=[123]
    i
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end

        [GPx, GPy] = gradient(phase_map);
        GPx = GPx  ./ (400*10^(-6));
        GPy = GPy  ./ (400*10^(-6));
        
        pgd_data(j,1,:,:) = GPx;
        pgd_data(j,2,:,:) = GPy;
    end
end

pgd_1 = [];
pgd_2 = [];
pgd_3 = [];
for i=1:641
    x = squeeze(pgd_data(i,:,:,:));
    y = sqrt(x(1,:,:).^2+x(2,:,:).^2);
    pgd_1(i,:,:) = squeeze(y);
end
pgd_1 = reshape(pgd_1, [size(pgd_1,1) 50]);
for i=1:size(pgd_1,1)
    pgd_2 = [pgd_2 squeeze(mean(pgd_1(i,~isnan(pgd_1(i,:))),2))];
end
pgd_t = reshape(pgd_data, [size(pgd_data,1) 2 50]);
for i=1:size(pgd_1,1)
    pgd_3(i,1,:) = pgd_t(i,1,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
    pgd_3(i,2,:) = pgd_t(i,2,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
end
pgd_3 = mean(pgd_3,3);
pgd_3 = sqrt(pgd_3(:,1).^2+pgd_3(:,2).^2);

PGD = pgd_3 ./ pgd_2';

i = 123



histogram(PGD(Time>=0),10);
hold on
histogram(PGD(Time<0),10);
title(['PGD - Trial ' num2str(i)],'interpreter','latex');
xlabel('PGD','interpreter','latex');
ylabel('Number','interpreter','latex');
legend({'After Cue' 'Before Cue'},'interpreter','latex');
export_fig(['PGD_Trial' num2str(50) '2.png'],'-r100');




function plotRawSigPinkNoise(f, pinkNoise, rawSignal, dePinkedSignal, targetElectrode)
    % plot raw data
    figure;
    subplot(1, 3, 1);
    hold on;
    plot(10.^f, 10.^rawSignal, 'LineWidth', 2);
    hold off;
    title("Raw Data - Electrode = " + targetElectrode, 'interpreter', 'latex');
    xlabel('Frequency(Hz)', 'interpreter', 'latex');
    ylabel('Power', 'interpreter', 'latex');
    grid on; grid minor;

    % plot pink noise-free data
    subplot(1, 3, 2);
    hold on;
    plot(10.^f, 10.^dePinkedSignal, 'LineWidth', 2);
    hold off;
    title("Pink Noise Free - Electrode = " + targetElectrode, 'interpreter', 'latex');
    xlabel('Frequency(Hz)', 'interpreter', 'latex');
    ylabel('Power', 'interpreter', 'latex');
    grid on; grid minor;
% plot raw data, pink noise, and pink noise-free signals
subplot(1, 3, 3);
hold on;
plot(10.^f, 10.^rawSignal, 'k', 'LineWidth', 2); % Raw data in black
plot(10.^f, 10.^pinkNoise, 'Color', [1 0.6 0.6], 'LineWidth', 2); % Pink noise in pink
plot(10.^f, 10.^dePinkedSignal, 'b', 'LineWidth', 2); % Pink noise-free in blue
hold off;
title("Combined - Electrode = " + targetElectrode, 'interpreter', 'latex');
xlabel('Frequency(Hz)', 'interpreter', 'latex');
ylabel('Power', 'interpreter', 'latex');
legend('Raw Signal', 'Pink Noise', 'Pink Noise Free');
grid on; grid minor;
end