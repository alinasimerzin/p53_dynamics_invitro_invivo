%% Find peaks in EMT6 cultured cells:
GoodCells = [];
NumPeaks = [];
myfig = 0;
counter = 0;
% tempYFP=traces_YFP (:, 1:44);
% temptimepoints = timepoints (1:44);
All(size(traces_YFP,1)) = struct();

for i=1:size(traces_YFP,1)
    sample = smooth(traces_YFP(i,:),0.08, 'loess');
    
    [pks, locs, W, P] = findpeaks(sample,timepoints,'MinPeakDistance',2.5,'MinPeakWidth',1.3,'MinPeakProminence',20);
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
% 
%         if W(1) < 4.25 
%             set(gca,'Color','g')
%         elseif W(1) >= 4.25 
%             set (gca,'Color',[0.7 0.7 0.7])
    end
 
    
    All(i).pks = pks';
    All(i).locs = locs;
    All(i).W = W;
    All(i).P = P';

    
   
end

%% plot: number of peaks
frame = cell2mat(annotation(:,2));

Numpeaks = [];
for i = 1:length(All)
    Numpeaks = [Numpeaks length(All(i).pks)];
end


SparseP = Numpeaks(frame>40 & frame<48);
DenseP = Numpeaks(frame>30 & frame<41);


fig = figure;
subplot(1,2,1)
histogram(SparseP,[0:8], 'facecolor',[0 0 0],'facealpha',.5,'edgecolor','k');
set (gca, 'Fontsize', 14); title ('Sparse + Toposar'); ylim([0 200]); ylabel('Frequency');
subplot(1,2,2)
histogram(DenseP,[0:8],'facecolor',[0 0 0],'facealpha',.5,'edgecolor','k'); 
title ('Dense + Toposar'); set(gca, 'Fontsize', 14); ylim([0 200])

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
SparseT = Time1(frame>40 & frame<51);
DenseT = Time1(frame>30 & frame<41);

figure
subplot(1,2,1)
histogram(SparseT,[0.5:1:20],'facecolor','k','facealpha',.5,'edgecolor','k')
title ('Sparse + Toposar'); set (gca, 'Fontsize', 14); ylim([0 100]);
ylabel('Frequency')
subplot(1,2,2)
histogram(DenseT,[0.5:1:20],'facecolor','k','facealpha',.5,'edgecolor','k')
title('Dense + Toposar'); set (gca, 'Fontsize', 14); ylim([0 100]);

data1 = {SparseT,DenseT};
figure; 
handles=plotSpread(data1, 'showMM',5); set(handles{1},'color','k'); box on;
set(gca, 'FontSize', 14, 'XTickMode','manual', 'XTick',[1,2],'XTickLabel',{'Sparse','Dense'});
ylabel('Time of 1st peak (hr)'); ylim([0 8])


n=[];
label={'Sparse','Dense'};
lengths = [length(SparseT) length(DenseT)];
for i = 1:2
    n = [n ; repmat(label(i),lengths(i),1)];
end
figure
boxplot(Time1',n)
ylabel('Time of 1st pulse')



label={'Control','Dense Old'};
lengths = [length(ControlT) length(DenseOldT)];
for i = 1:2
    n = [n ; repmat(label(i),lengths(i),1)];
end

figure
boxplot(Time1',n)
ylabel('Time of 1st pulse')

[~, ~,stats]=anova1(Time1',n); multcompare(stats)
set(gcf, 'Position', [144   942   496   287])

% Boxplot 

x1 = Time1(frame>50);
x2 = Time1(frame>30 & frame<41);

x = [x1'; x2']; 


% Estimating the length of the groups 
lx1 = length (x1); 
lx2 = length (x2);


% Generating x-axis labels 
g1 = repmat({'Sparse'},lx1,1);
g2 = repmat({'Dense'},lx2,1);

g = [g1; g2];


% figure of the plots 
% On each box, the central mark indicates the median, and the bottom and top edges of the box indicate the 25th and 75th percentiles, respectively. 
% The whiskers extend to the most extreme data points not considered outliers, and the outliers are plotted individually using the '+' symbol.
figure
boxplot(x,g, 'BoxStyle','outline','Widths',0.5, 'Colors','k'); set(gca, 'FontSize', 14);

% change Median color and appearance
h = findobj(gca,'tag','Median');
set(h,'linestyle','-');
set(h,'Color',[1 0 0])
ylabel('Time of 1st peak')

figure
boxplot(z,f, 'BoxStyle','outline','Widths',0.5, 'Colors','k'); set(gca, 'FontSize', 14);

% change Median color and appearance
h = findobj(gca,'tag','Median');
set(h,'linestyle','-');
set(h,'Color',[1 0 0])
ylabel('Time of 1st peak')

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


SparseW = Width1(frame>40 & frame<51);
DenseW = Width1(frame>30 & frame<41);

figure
subplot(1,2,1)
histogram(SparseW,[0:1:5],'facecolor','k','facealpha',.5,'edgecolor','k')
title('Sparse + xray');  set (gca, 'Fontsize', 14); ylim([0 200]);ylabel('Frequency')
subplot(1,2,2)
histogram(DenseW,[0:1:5],'facecolor','k','facealpha',.5,'edgecolor','k')
title('Dense + xray');  set (gca, 'Fontsize', 14); ylim([0 200]);


data2 = {SparseW,DenseW};
figure; 
handles=plotSpread(data2, 'showMM',5); set(handles{1},'color','k'); box on;
set(gca, 'FontSize', 14, 'XTickMode','manual', 'XTick',[1,2],'XTickLabel',{'Sparse','Dense'});
ylabel('Width of 1st pulse (hr)'); ylim([0 6]);


n=[];
label={'Sparse','Dense'};
lengths = [length(SparseW) length(DenseW)];
for i = 1:2
    n = [n ; repmat(label(i),lengths(i),1)];
end
figure
boxplot(Width1',n)
ylabel('Width of 1st pulse')

[~, ~,stats]=anova1(Width1',n); multcompare(stats)
set(gcf, 'Position', [144   942   496   287])

% Boxplot for sparse and dense HMFL
y1 = Width1(frame>50);
y2 = Width1(frame>30 & frame<41);

y = [y1'; y2']; 



% Estimating the length of the groups 
ly1 = length (y1); 
ly2 = length (y2);



% Generating x-axis labels 
e1 = repmat({'Sparse'},ly1,1);
e2 = repmat({'Dense'},ly2,1);
e = [e1; e2];


% figure of the plots 
% On each box, the central mark indicates the median, and the bottom and top edges of the box indicate the 25th and 75th percentiles, respectively. 
% The whiskers extend to the most extreme data points not considered outliers, and the outliers are plotted individually using the '+' symbol.

boxplot(y,e, 'BoxStyle','outline','Widths',0.5, 'Colors','k'); set(gca, 'FontSize', 14);

% change Median color and appearance
h = findobj(gca,'tag','Median');
set(h,'linestyle','-');
set(h,'Color',[1 0 0])
ylabel('Duration of 1st pulse')


boxplot(c,k, 'BoxStyle','outline','Widths',0.5, 'Colors','k'); set(gca, 'FontSize', 14);

% change Median color and appearance
h = findobj(gca,'tag','Median');
set(h,'linestyle','-');
set(h,'Color',[1 0 0])
ylabel('Duration of 1st pulse')


[t,p]=ttest2(y1',y2', 'VarType', 'unequal')

%% Compare in vitro and in vivo data 

% Create a cell array to store Time of the peak values
Time2_topo = cell(273, 1); % Data from findpeaks All.loc

extracted_values = cell(273, 1); 

% Loop through each cell of Time2_topo
for i = 1:273
    % Check if the cell is empty or has fewer than 2 elements
    if ~isempty(Time2_topo{i}) && numel(Time2_topo{i}) >= 2
        % Extract the second value from the current row
        extracted_values{i} = Time2_topo{i}(2);
    else
        % Handle the case where there is no second value
        extracted_values{i} = NaN;  % You can use a placeholder value like NaN
    end
end

% Logical indexing to find cells where Numpeaks is equal to 2
indices_numpeaks_2 = (NumpeaksSpTopo == 2);

% Extract Time1 and Time2 values only where Numpeaks is equal to 2
NumpeaksSpTopo_only2 = extracted_values(indices_numpeaks_2);
Time1_numpeaks_2 = TimeSpTopo(indices_numpeaks_2);
Time2_numpeaks_2 = cell2mat(NumpeaksSpTopo_only2)'; 

% Calculate the difference Time2 - Time1 for cells with Numpeaks = 2
time_difference_numpeaks_2 = Time2_numpeaks_2 - Time1_numpeaks_2;

time_difference_only2pkscells_invitro = time_difference_numpeaks_2;
time_difference_only2pkscells_invivo = [];



%% plot bars 
% Calculate mean and standard deviation for each dataset
invitro_data = mean(time_difference_only2pkscells_invitro);
std_invitro_data = std(time_difference_only2pkscells_invitro);

invivo_data = mean(time_difference_only2pkscells_invivo);
std_invivo_data = std(time_difference_only2pkscells_invivo);



% Create a bar graph
figure;
bar([invitro_data, invivo_data], 'BarWidth', 0.6, 'FaceColor', [0.2, 0.2, 0.2]);
% Add error bars using the standard deviation
hold on;
errorbar([1, 2], [invitro_data, invivo_data], [std_invitro_data, std_invivo_data], 'k.', 'LineWidth', 0.5);
ylabel('{\Delta}T(peak2-peak1)');
xticks([1, 2]);
xticklabels({'in vitro', 'in vivo'});
set (gca, 'Fontsize', 14);box on;

% Set the figure properties
fig = gcf; % Get current figure handle
fig.PaperOrientation = 'landscape'; % Set orientation to landscape
fig.PaperUnits = 'normalized'; % Set paper units to normalized
fig.PaperPosition = [0 0 1 1]; % Set paper position to full figure size

% Save the figure as a PDF file
saveas(fig, 'time difference between pulses.pdf', 'pdf')

%% Extract durations of first and second 

Width2_topo = cell(273, 1);  % paste from All3.W

Width_extracted_values = cell(273, 1);

% Loop through each cell of Time2
for i = 1:273
    % Check if the cell is empty or has fewer than 2 elements
    if ~isempty(Width2_topo{i}) && numel(Width2_topo{i}) >= 2
        % Extract the second value from the current row
        Width_extracted_values{i} = Width2_topo{i}(2);
    else
        % Handle the case where there is no second value
        Width_extracted_values{i} = NaN;  % You can use a placeholder value like NaN
    end
end

% Logical indexing to find cells where Numpeaks is equal to 2
indices_numpeaks_2 = (NumpeaksSpTopo == 2);

% Extract Time1 and Time2 values only where Numpeaks is equal to 2
WidthSpTopo_only2_1peak = Width_extracted_values(indices_numpeaks_2);
WidthSpTopo_only2_2ndpeak = 
Time1_numpeaks_2 = TimeSpTopo(indices_numpeaks_2);
Time2_numpeaks_2 = cell2mat(NumpeaksSpTopo_only2)'; 

