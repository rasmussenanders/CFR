%% 
% CFR anal

clear; clc; close all;
folder = '/Users/med-anr/Google Drive/Documents (GDrive)/Projects/Current/CFR shapes/Analysis';
figures = '/Users/med-anr/Google Drive/Documents (GDrive)/Projects/Current/CFR shapes/Analysis/Figures';
functions = '/Users/med-anr/Google Drive/Documents (GDrive)/Lab/MATLAB/ARscripts';
addpath(genpath(functions));
cd(folder); addpath(genpath(folder)); % Goto folder
load("CFR_6.mat") % Load the data

% Graphics
set(0, 'DefaultFigureRenderer', 'painters'); % Save figures in vector format

% Define the colormap
lu = [
    [0, 0, 0] / 255;     % Black
    [156, 97, 20] / 255;   % Brown
    [185, 211, 220] / 255; % Light blue
    [223, 196, 199] / 255; % Pink
    [0, 0, 128] / 255;     % Blue
    [173, 202, 184] / 255; % Turquoise
    [214, 210, 196] / 255; % Beige  
    [191, 184, 175] / 255  % Dark gray
    [191, 184, 175] / 255  % Dark gray
];

% Alternative colormap (without black)
lu2 = [
    [0, 0, 128] / 255; % Blue
    [156, 97, 20] / 255; % Brown
    [223, 196, 199] / 255; % Pink
    [185, 211, 220] / 255; % Light blue
    [173, 202, 184] / 255; % Turkos
    [214, 210, 196] / 255; % Beige
    [191, 184, 175] / 255; % Dark grey
    ]; 

% 
Cells = categorical({ ...
    'e785pc1', ...
    'e785PC1', ...
    'e785pc2.2', ...
    'e785pc3', ...
    'e793pc6' ...
    'e806pc2', ...
    'e806pc6', ...
    'e822pc1', ...
    'e822pc3.1', ...
    'e822pc3.2', ...
    'e822pc4', ...
    'e822pc5', ...
    'e822pc6', ...
    'e822pc7', ...
    'e822pc9', ...
    'e825pc4', ...
    'e825pc5', ...
    'e825pc6', ...
    'e830pc1', ...
    'e830pc2', ...
    'e835pc3', ...
    'e835pc4', ...
    'e835pc6', ...
    'e839pc5', ...
    'e839pc8', ...
    'e844pc2', ...
    'e1545pc1', ...
    'e1545', ...
    'e1551pc1', ...
    'e1551pc2', ...
    'e1551pc3', ...
    });

%Cells = unique(EPSPdata.Cell); % Find all unique cells in the entire dataset  


%% Color graph of all EPSPs

AllEPSPs = nan(numel(Cells),1000); % Create variable

% Loop through all cells and get the number of EPSPs in that cell
for j = 1 : numel(Cells)
    temp = EPSPdata.Epsp(EPSPdata.Cell == Cells(j)); % Find all observations in the first cell
    AllEPSPs(j,1:numel(temp)) = temp;
end

AllEPSPs(:,50:1000) = []; % Remove columns with only nans

% Create figure (including zeroes)
figure()
colormap(lu) % Set the colormap to the specified colors
imagesc(AllEPSPs, 'AlphaData', ~isnan(AllEPSPs)) % Plot the data
set(gca, 'TickDir', 'out'); box off; set(gca,'color','none');
set(gca,'YDir','normal') 
set(gca, 'TickDir', 'out')
xticks(5:5:50);
ylabel('Cell'); xlabel('CFR responses')
colorbar
str = 'All CFRs';
title(str)
cd(figures); saveas(gcf, str, 'pdf'); cd(folder)

% Create a new variable with zeroes removed
AllEPSPs(AllEPSPs == 0) = NaN;  % Replace zeroes with NaNs
AllEPSPsNoZero = nan(numel(Cells),49); 

% Loop through cells and gather the relevant data
for j = 1 : numel(Cells)
    temp = AllEPSPs(j,:); temp = temp(~isnan(temp));
    AllEPSPsNoZero(j,1:numel(temp)) = temp;
end

% Update the colormap
lu = [
    [156, 97, 20] / 255;   % Brown
    [185, 211, 220] / 255; % Light blue
    [223, 196, 199] / 255; % Pink
    [0, 0, 128] / 255;     % Blue
    [173, 202, 184] / 255; % Turquoise
    [214, 210, 196] / 255; % Beige  
    [191, 184, 175] / 255  % Dark gray
    [191, 184, 175] / 255  % Dark gray
];

figure() % Create figure
colormap(lu)  % Set the colormap to the specified colors
imagesc(AllEPSPsNoZero, 'AlphaData', ~isnan(AllEPSPsNoZero))
set(gca, 'TickDir', 'out'); box off; set(gca,'color','none');
set(gca,'YDir','normal') 
set(gca, 'TickDir', 'out')
xticks(5:5:50);
ylabel('Cell'); xlabel('CFR responses')
colorbar
str = 'All CFRs, zeroes excluded';
title(str)
cd(figures); saveas(gcf, str, 'pdf'); cd(folder)



%% Barchart showing the number of EPSPs in the selected cells  

Cells = unique(EPSPdata.Cell); % Find all unique cells in the entire dataset   
AllEPSPs = nan(100, numel(Cells)); % Create variable 
SpontEPSPs = nan(100, numel(Cells)); % Create second variable

% Go through all cells and find the number of EPSPs in different conditions

for j = 1:numel(Cells) 
    cellData = Cells(j);
    temp = EPSPdata.Epsp(EPSPdata.Cell == cellData);
    AllEPSPs(1:numel(temp), j) = temp;

    temp = EPSPdata.Epsp(EPSPdata.Cell == cellData & EPSPdata.Stimtyp=='Spon'); % Find only spontaneous
    SpontEPSPs(1:numel(temp), j) = temp;
end

% Remove zeroes
AllEPSPs(AllEPSPs==0) = nan;
SpontEPSPs(AllEPSPs==0) = nan;

% Compute mean EPSPs and standard deviations
meanEPSPs = nanmean(AllEPSPs);
stdEPSPs = nanstd(AllEPSPs);
meanSpontEPSPs = nanmean(SpontEPSPs);
stdSpontEPSPs = nanstd(SpontEPSPs);

% Sort the Cells based on mean values
[sortedMeanEPSPs, sortedIndices] = sort(meanEPSPs);
sortedCells = Cells(sortedIndices);
[sortedMeanSpontEPSPs, sortedIndices] = sort(meanSpontEPSPs);
sortedSpontCells = Cells(sortedIndices);

% Plotting bar chart with error bars
figure;
hold on;
bar(sortedMeanEPSPs, 'FaceColor', 'b', 'EdgeColor', 'k');
errorbar(1:numel(Cells), sortedMeanEPSPs, stdEPSPs(sortedIndices), 'k', 'LineStyle', 'none');

ylabel('EPSP');
title('EPSP Bar Chart (Sorted by Mean EPSPs)');
set(gca, 'TickDir', 'out'); box off; set(gca,'color','none');
str = ['EPSP distribution in 31 Purkinje cells']; title(str); 
saveas(gcf, str, 'pdf');

% Do the same for spontaneous
figure;
hold on;
bar(sortedMeanSpontEPSPs, 'FaceColor', 'b', 'EdgeColor', 'k');
errorbar(1:numel(Cells), sortedMeanSpontEPSPs, stdSpontEPSPs(sortedIndices), 'k', 'LineStyle', 'none');
xlim([0 20])
ylabel('EPSP');
title('Spontaneous EPSP Bar Chart (Sorted by Mean EPSPs)');
set(gca, 'TickDir', 'out'); box off; set(gca,'color','none');
str = ['EPSP distribution in 19 Purkinje cells']; title(str); 
saveas(gcf, str, 'pdf');

%% Test difference in the number of spontaneous EPSPs using one way anova


SpontTest = SpontEPSPs(:,sum(isnan(SpontEPSPs))<100)
SpontCell = SpontTest; 
for j =1: size(SpontCell,2);
    SpontCell(:,j) = j;
end

AnovaData = [SpontTest(:) SpontCell(:)];
AnovaData(isnan(AnovaData(:,1)),:) = [];
% [p,tbl,stats] = anova1(AnovaData(:,1),AnovaData(:,2));
% multcompare(stats);


% Perform one-way ANOVA
[p, tbl, stats] = anova1(AnovaData(:, 1), AnovaData(:, 2));

% Compute effect size (eta-squared)
SSGroups = tbl{2, 2};
SSTotal = tbl{4, 2};
etaSquared = SSGroups / SSTotal;

% Display results and effect size
disp(['F = ', num2str(etaSquared)]);
disp(['Effect size (eta-squared): ', num2str(tbl{2,5}), ' p = ', num2str(tbl{2,6})]);

% Perform multiple comparisons (optional)
multcompare(stats);


%% Compare Spontaneous and elicited CFRs

SponElicited = nan(31,2); % define variable

% Find mean number of EPSPs for 
for j = 1 : numel(Cells)
    cell = Cells(j);
    Spontemp = EPSPdata.Epsp(EPSPdata.Cell == cell & EPSPdata.Stimtyp == 'Spon')
    Eltemp = EPSPdata.Epsp(EPSPdata.Cell == cell & EPSPdata.Stimtyp == '08' & EPSPdata.Epsp ~=0)
    SponElicited(j,1) = mean(Spontemp);
    SponElicited(j,2) = mean(Eltemp);
end

sample = find(~isnan(SponElicited(:,1)) & ~isnan(SponElicited(:,2)));
SponElicited = SponElicited(sample,:)

x = SponElicited(:,1);
y = SponElicited(:,2);

% Calculate means and standard deviations
mean_x = mean(x);
mean_y = mean(y);
std_x = std(x);
std_y = std(y);

% Perform paired t-test
[h, p, ci, stats] = ttest(x, y, 'Alpha', 0.05, 'Tail', 'both');

% Calculate effect size (Cohen's d)
mean_diff = mean(x - y);
pooled_std = sqrt((std_x^2 + std_y^2) / 2);
effect_size = mean_diff / pooled_std;

% Display results
disp('Paired t-test:');
disp(['   t-value: ', num2str(stats.tstat)]);
disp(['   p-value: ', num2str(p)]);
disp(['Mean (x): ', num2str(mean_x)]);
disp(['Mean (y): ', num2str(mean_y)]);
disp(['Standard Deviation (x): ', num2str(std_x)]);
disp(['Standard Deviation (y): ', num2str(std_y)]);
disp(['Effect size (Cohen''s d): ', num2str(effect_size)]);

%% Linear mixed effects model to examine the effect of NO stimulation on EPSPs

% Linear mixed effects model with zeroes
Sample = EPSPs_Stim;
lme = fitlme(Sample, 'Epsp ~ StimstyrkaNO + (1|Cell) ') % Linear mixed effects model with NO stimulation strenght as fixed variable and cell as random variable

% Linear mixed effects model without zeroes
Sample = EPSPs_Stim(EPSPs_Stim.Epsp~=0,:);
lme = fitlme(Sample, 'Epsp ~ StimstyrkaNO + (1|Cell) ') % Linear mixed effects model with NO stimulation strenght as fixed variable and cell as random variable


%% Compare EPSPs between periorbital stimulation with and without NO

lu = [
    [0, 0, 128] / 255; % Blue
    [156, 97, 20] / 255; % Brown
    [223, 196, 199] / 255; % Pink
    [185, 211, 220] / 255; % Light blue
    [173, 202, 184] / 255; % Turkos
    [214, 210, 196] / 255; % Beige
    [191, 184, 175] / 255; % Dark grey
    ]; 

Title = 'Periorbital ESPSs with and without NO stimualtion zeroes included'; y = [-0.25 0.5];
Cells = unique(EPSPs_Stim.Cell);

NOstim = nan(numel(Cells),2); % Create variable to add values to, first column = periorbital only; second column; periorbital + NO

for j = 1:numel(Cells)
    NOstim(j,1) = mean(EPSPs_Stim.Epsp(EPSPs_Stim.Cell == Cells(j) & EPSPs_Stim.Stimtyp == "08" ));
    NOstim(j,2) = mean(EPSPs_Stim.Epsp(EPSPs_Stim.Cell == Cells(j) & EPSPs_Stim.Stimtyp == "14+08" ));
end

% Plot data using raincloud
A = NOstim(:,1); Name_A = 'Periorbital';
B = NOstim(:,2); Name_B = 'NO + Periorbital';

% Do a paired t-test
[h,p,CI,stats] = ttest(A,B)

str = ['A (', num2str(sum(~isnan(A))), '), Mean = ', num2str(nanmean(A)), ' ± ', num2str(nanstd(A))]; disp(str)
str = ['B (', num2str(sum(~isnan(B))), '), Mean = ', num2str(nanmean(B)), ' ± ', num2str(nanstd(B))]; disp(str)

% Raincloud figure 
figure()
dat = A; color = lu(1,:); h1 = raincloud_plot(dat, 'box_on', 1, 'color', color, 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 1);
dat = B; color = lu(2,:); h2 = raincloud_plot(dat, 'box_on', 1, 'color', color, 'alpha', 0.5, 'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);

ylim(y)
set(gca, 'TickDir', 'out'); box off; set(gca,'color','none'); ylabel('Boxplot and Probability density'); xlabel('mean EPSPs')

legend([h1{1} h2{1}], { ...
    [Name_A,' ', num2str(nanmean(A)), ' ± ', num2str(nanstd(A))], ... 
    [Name_B,' ', num2str(nanmean(B)), ' ± ', num2str(nanstd(B))], ...
    }, "Location", "northeast", "Box", "off");

str = [Title,' p = ', num2str(p)]; title(str); 
cd(figures); saveas(gcf, Title, 'pdf'); cd(folder)


