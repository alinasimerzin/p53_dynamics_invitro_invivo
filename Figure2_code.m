% Finding peaks in Sparse and Dense MCF7 cultures. 

%Find peaks
GoodCells = [];
NumPeaks = [];
myfig = 0;
counter = 0;

All(size(traces_YFP,1)) = struct();

for i=1:size(traces_YFP,1)
    sample = smooth(traces_YFP(i,:),6);
  
    
    [pks, locs, W, P] = findpeaks(sample,timepoints,'MinPeakDistance',2.1,'MinPeakWidth',1.3,'MinPeakProminence',8);
    if ~isempty(pks) 
        numpeaks = length(pks);
        temp = zeros(numpeaks,1);
        temp(P ./ pks > 0.4) = 1; %get rid of peaks with low prominence - for example - double peak with low trough 
        
        pks = pks(temp == 1);
        locs = locs(temp==1);
        W = W(temp==1);
        P = P(temp==1);

        
        if sum(temp)  > 0
            GoodCells = [GoodCells i];
        end
        
        if sum(temp) == 0
            set(gca,'Color','r')
        elseif sum(temp) == 1 
            set(gca,'Color',[1 0.8 0.8]) %pink
        elseif sum(temp) == 2
            set(gca,'Color',[0.8 0.8 1]) %purple
        elseif sum(temp) == 3
            set(gca, 'Color',[0.5 1 1]) %cyan
        elseif sum(temp) > 3
            set(gca,'Color','g')
        end
    end
    
    All(i).pks = pks';
    All(i).locs = locs;
    All(i).W = W;
    All(i).P = P';  
end

%% plot peaks data
Numpeaks = [];
for i = 1:length(All)
    Numpeaks = [Numpeaks length(All(i).pks)];
end

% Define Sparse and Dense based on the positions in the movie.
NumPulsesCtrl = Numpeaks(frame>10 & frame<21);
NumPulsesDense = [Numpeaks(frame>30 & frame<41),  Numpeaks(frame>50 & frame<61)];


figure(1);
histogram(NumPulsesCtrl,[-0.5 0.5  1.5  2.5  3.5  4.5 5.5],'facecolor',[0.2 1 0.2],'facealpha',.5,'edgecolor','k'); ylim ([0 200])
% axis tight
%title('Control')
ylabel('# cells')
xlabel('# peaks')
set (gca, 'Fontsize', 12);

figure(2);
histogram(NumPulsesDense,[-0.5 0.5  1.5  2.5  3.5  4.5 5.5],'facecolor',[0.2 1 0.2],'facealpha',.5,'edgecolor','k'); ylim ([0 200])
% axis tight
%title('Control')
ylabel('# cells')
xlabel('# peaks')
set (gca, 'Fontsize', 12);

%% peak amplitudes
Amplitute1 = [];
for i = 1:length(All)
    temp = All(i).pks;
    if ~isempty(temp)
    Amplitute1 = [Amplitute1 temp(1)];
    else
        Amplitute1 = [Amplitute1 NaN];
    end
end

AmplSparse = Amplitute1(frame>10 & frame<21);
AmplDense = [Amplitute1(frame>30 & frame<41), Amplitute1(frame>50 & frame<61)];
AmplSparse = AmplSparse(~isnan(AmplSparse));
AmplDense = AmplDense(~isnan(AmplDense));
AvAmplS = mean(AmplSparse)
AvAmplD = mean(AmplDense)

data1 = {AmplSparse,AmplDense};
figure; 
handles=plotSpread(data1, 'showMM',5); set(handles{1},'color','k'); box on;
set(gca, 'FontSize', 14, 'XTickMode','manual', 'XTick',[1,2],'XTickLabel',{'Sparse','Dense'});
ylabel('Amplitude of 1st peak (a.u)');

%U-test
[pAmp,hAmp,statsAmp] = ranksum(AmplSparse,AmplDense); %pval:  1.214435599899223e-14
%t-test
[pvalAmp, pvalAmp] = ttest2(AmplSparse,AmplDense); %pval:  3.476866803900628e-20

figure;
 plot_histogram_shaded(AmplSparse, 'alpha', 0.3, 'bins', 10, 'color', [0 0 0],...
   'edges', [-10, 0, 200, 400,600, 800, 1000, 1200, 1400, 1600, 1800],   'normalization', 'probability'); xlim ([0 1000])
hold on;
plot_histogram_shaded(AmplDense, 'alpha', 0.3, 'bins', 10, 'color', [1 0 0],...
 'edges', [-10, 0, 200, 400,600, 800, 1000, 1200, 1400, 1600, 1800],   'normalization', 'probability');
set (gca, 'Fontsize', 14, 'box', 'on'); ylim([0 1]);
xlabel('Time of 1st peak (hr)');
ylabel('Probability');



%% peak 1 timings
Time1 = [];
for i = 1:length(All)
    temp = All(i).locs;
    if ~isempty(temp)
    Time1 = [Time1 temp(1)];
    else
        Time1 = [Time1 NaN];
    end
end

TimeSparse = Time1(frame>10 & frame<15);
TimeSparse1 = TimeSparse(TimeSparse<7);
TimeDense = [Time1(frame>30 & frame<41), Time1(frame>50)];
TimeSparse = TimeSparse (~isnan(TimeSparse));
TimeDense = TimeDense (~isnan(TimeDense));


data2 = {TimeSparse,TimeDense};
figure; 
handles=plotSpread(data2, 'showMM',5); set(handles{1},'color','k'); box on;
set(gca, 'FontSize', 14, 'XTickMode','manual', 'XTick',[1,2],'XTickLabel',{'Sparse','Dense'});
ylabel('Time of 1st peak (hr)');

%U-test
[pTime,hTime,statsTime] = ranksum(TimeSparse,TimeDense); %pval:   2.048996104855806e-37
%t-test
[pvalTime, pvalTime] = ttest2(TimeSparse,TimeDense); %pval:  9.556518184126213e-29


figure;
 plot_histogram_shaded(TimeSparse, 'alpha', 0.3, 'bins', 10, 'color', [0 0 0],...
   'edges', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],   'normalization', 'probability'); xlim ([0 10])
hold on;
plot_histogram_shaded(TimeDense, 'alpha', 0.3, 'bins', 10, 'color', [1 0 0],...
 'edges', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],   'normalization', 'probability');
set (gca, 'Fontsize', 14, 'box', 'on'); ylim([0 1]);
xlabel('Time of 1st peak (hr)');
ylabel('Probability');

%% peak 1 widths (half max prominence)
Width1 = [];
for i = 1:length(All)
    temp = All(i).W;
    if ~isempty(temp)
    Width1 = [Width1 temp(1)];
    else
        Width1 = [Width1 NaN];
    end
end

WidthCtrl = Width1(frame>10 & frame<21);
WidthDense = [Width1(frame>30 & frame<41), Width1(frame>50)];


WidthCtrl = WidthCtrl(~ isnan(WidthCtrl));
WidthDense = WidthDense(~ isnan(WidthDense));

%U-test
[pWidth,hWidth,statsWidth] = ranksum(WidthCtrl,WidthDense); %pval:   2.470297985439836e-36
%t-test
[pvalWidth, pvalWidth] = ttest2(WidthCtrl,WidthDense); %pval:  6.034677008044913e-37



data3 = {WidthCtrl,WidthDense};
figure; 
handles=plotSpread(data3, 'showMM',5); set(handles{1},'color','k'); box on;
set(gca, 'FontSize', 14, 'XTickMode','manual', 'XTick',[1,2],'XTickLabel',{'Sparse','Dense'});
ylabel('Duration of 1st peak (hr)'); ylim ([0 25]);


figure;
plot_histogram_shaded(WidthCtrl, 'alpha', 0.3, 'bins', 10, 'color', [0 0 0],...
   'edges', [0, 1, 2, 3, 4, 5, 6,7,8,9,10],   'normalization', 'probability'); xlim ([0 10])
hold on;
plot_histogram_shaded(WidthDense, 'alpha', 0.3, 'bins', 10, 'color', [1 0 0],...
 'edges', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],   'normalization', 'probability');
set (gca, 'Fontsize', 14, 'box', 'on'); ylim([0 1]);
xlabel('Duration of 1st peak (hr)');
ylabel('Probability');

%% Pulse duration based on autocorrelations 

mean_p53 = mean(traces_YFP,2);
[~,p53_idx] = sortrows(mean_p53);

% Sparse 
traces_Sparse = traces_YFP(group_number == 1,:);
p53_AC = autocorrelationMatrix(traces_Sparse(:,1:50),1:30); 
%1:50 - to see only 2 first peaks (there is a limit on the time in which to perform auto-correlation. 1:20 frames - 5 hours
p53_mean_AC = nanmean(p53_AC,2);
[p53_sorted_mean_AC, p53_sorted_idx_AC] = sortrows(p53_mean_AC);
p53_sorted_AC = p53_AC(p53_sorted_idx_AC,:);
frame_sort53 = frame(p53_idx);
frame_sorted_AC = frame_sort53(p53_sorted_idx_AC);

figure; 
plotnfill_auto_quantiles_exclude(1:30,p53_sorted_AC(:,:),0.25,'black', 0)
pbaspect([1 1 1]);
xlabel('Time (hr)');
ylabel('Autocorrelation');
set(gca,'XTickMode','manual');
set(gca,'XTick',[0,4,8,12,16, 20, 24,28]);
set(gca,'XtickLabels',[0,1,2,3,4,5,6,7]);
set (gca, 'Fontsize', 14);box on;

% Set the figure properties
fig = gcf; % Get current figure handle
fig.PaperOrientation = 'landscape'; % Set orientation to landscape
fig.PaperUnits = 'normalized'; % Set paper units to normalized
fig.PaperPosition = [0 0 1 1]; % Set paper position to full figure size

% Save the figure as a PDF file
saveas(fig, '2023-12-10_Sparse_autocorr.pdf', 'pdf')

% Dense cultures 
traces_Dense = traces_YFP(group_number == 3,:);
p53_AC2 = autocorrelationMatrix(traces_Dense(:,1:50),1:30); 
%1:50 - to see only 2 first peaks (there is a limit on the time in which to perform auto-correlation. 1:20 frames - 5 hours
p53_mean_AC2 = nanmean(p53_AC2,2);
[p53_sorted_mean_AC2, p53_sorted_idx_AC2] = sortrows(p53_mean_AC2);
p53_sorted_AC2 = p53_AC2(p53_sorted_idx_AC2,:);

figure; 
plotnfill_auto_quantiles_exclude(1:30,p53_sorted_AC2(:,:),0.25,'red', 0)
pbaspect([1 1 1]);
xlabel('Time (hr)');
ylabel('Autocorrelation');
set(gca,'XTickMode','manual');
set(gca,'XTick',[0,4,8,12,16,20,24,28]);
set(gca,'XtickLabels',[0,1,2,3,4,5,6,7]);
set (gca, 'Fontsize', 14);box on;

% Set the figure properties
fig = gcf; % Get current figure handle
fig.PaperOrientation = 'landscape'; % Set orientation to landscape
fig.PaperUnits = 'normalized'; % Set paper units to normalized
fig.PaperPosition = [0 0 1 1]; % Set paper position to full figure size

% Save the figure as a PDF file
saveas(fig, '2023-12-10_Dense_autocorr.pdf', 'pdf')

%% plotting autocorrelaion pattern (duration)

mean_p53 = mean(traces_YFP,2);
[~,p53_idx] = sortrows(mean_p53);

% All 
p53_AC3 = autocorrelationMatrix(traces_YFP(:,1:50),1:30); 
%1:50 - to see only 2 first peaks (there is a limit on the time in which to perform auto-correlation. 1:20 frames - 5 hours
p53_mean_AC3 = nanmean(p53_AC3,2);
[p53_sorted_mean_AC3, p53_sorted_idx_AC3] = sortrows(p53_mean_AC3);
p53_sorted_AC3 = p53_AC3(p53_sorted_idx_AC3,:);
frame_sort53 = frame(p53_idx);
frame_sorted_AC3 = frame_sort53(p53_sorted_idx_AC3);


num_cells = size(p53_AC3, 1);  % Number of cells
time_points = size(p53_AC3, 2);  % Number of time points

% Initialize a vector to store the durations
pattern_durations = zeros(num_cells, 1);

% Loop through each cell
for i = 1:num_cells
    % Extract autocorrelation trace for the current cell
    autocorr_trace = p53_AC3(i, :);
    
    % Find the first peak after time 0
    [~, peak_index] = findpeaks(autocorr_trace, 'MinPeakDistance',0.7,'MinPeakWidth',0.7,'MinPeakProminence',0.1); % Adjust the threshold as needed
    
    % Check if a peak is found
    if ~isempty(peak_index)
        % Extract the time of the first peak
        first_peak_time = peak_index(1);
        
        % Calculate the duration by subtracting time 0
        pattern_duration = first_peak_time;
        
        % Store the duration in the vector
        pattern_durations(i) = pattern_duration;
    else
        % Handle the case where no peak is found (pattern duration is undefined)
        pattern_durations(i) = NaN;
    end
end



Sparse_pattern = pattern_durations(frame_sorted_AC3>10 & frame_sorted_AC3<21,:,:);
Sparse_valid_durations = Sparse_pattern(~isnan(Sparse_pattern));
Sparse_valid_durations = Sparse_valid_durations./4; % hours 

Dense_pattern = pattern_durations(frame_sorted_AC3>30 & frame_sorted_AC3<41,:,:);
Dense_valid_durations = Dense_pattern(~isnan(Dense_pattern));
Dense_valid_durations = Dense_valid_durations./4; % hours 

figure;
plot_histogram_shaded(Dense_valid_durations, 'alpha', 0.3, 'bins', 8, 'color', [1 0 0],...
   'edges', [-1, 3, 4, 5, 6, 7, 8],   'normalization', 'probability'); xlim ([0 10]);ylim ([0 1]);
hold on;
plot_histogram_shaded(Sparse_valid_durations, 'alpha', 0.3, 'bins', 8, 'color', [0 0 0],...
   'edges', [-1, 3, 4, 5, 6, 7, 8],   'normalization', 'probability'); xlim ([0 10]);ylim ([0 1]);
xlabel('Duration of 1st peak (hr)');
ylabel('Probability');
set (gca, 'Fontsize', 14);box on;
set(gca,'XTick',[0,1,2,3,4,5]);
set(gca,'XtickLabels',[0,1,2,3,4,5]);
set(gca,'YTick',[0,0.2,0.4,0.6,0.8,1]);
set(gca,'YtickLabels',[0,0.2,0.4,0.6,0.8,1]);
% Set the figure properties
fig = gcf; % Get current figure handle
fig.PaperOrientation = 'landscape'; % Set orientation to landscape
fig.PaperUnits = 'normalized'; % Set paper units to normalized
fig.PaperPosition = [0 0 1 1]; % Set paper position to full figure size

% Save the figure as a PDF file
saveas(fig, 'Autocorr_Durations_distr.pdf', 'pdf')


% Calculate mean and standard deviation for each dataset
mean_Sparse_autocorr = mean(Sparse_valid_durations);
std_Sparse_autocorr = std(Sparse_valid_durations);

mean_Dense_autocorr = mean(Dense_valid_durations);
std_Dense_autocorr = std(Dense_valid_durations);

% Create a bar graph
figure;
bar([mean_Dense_autocorr, mean_Sparse_autocorr], 'BarWidth', 0.6, 'FaceColor', [0.2, 0.2, 0.2]);
% Add error bars using the standard deviation
hold on;
errorbar([1, 2], [mean_Sparse_autocorr, mean_Dense_autocorr], [std_Sparse_autocorr, std_Dense_autocorr], 'k.', 'LineWidth', 1.5);
ylabel('Duration (h)');
xticks([1, 2]);
xticklabels({'Sparse', 'Dense'});
set (gca, 'Fontsize', 14);box on;

% Set the figure properties
fig = gcf; % Get current figure handle
fig.PaperOrientation = 'landscape'; % Set orientation to landscape
fig.PaperUnits = 'normalized'; % Set paper units to normalized
fig.PaperPosition = [0 0 1 1]; % Set paper position to full figure size

% Save the figure as a PDF file
saveas(fig, 'Autocorr_durations.pdf', 'pdf')


data5 = {Sparse_valid_durations,Dense_valid_durations};
figure; 
handles=plotSpread(data5, 'showMM',5); set(handles{1},'color','k'); box on;
set(gca, 'FontSize', 14, 'XTickMode','manual', 'XTick',[1,2],'XTickLabel',{'Sparse','Dense'});
ylabel('Duration of 1st peak (hr)'); ylim ([0 10]);



