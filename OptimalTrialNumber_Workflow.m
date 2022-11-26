%   OptimalTrialNumber_Workflow.m
%   2022/11/25
%   Updated 2022/06/27
%
%    If this code is used in a publication, please cite the manuscript:
%    "A Quantitative Method to Determine Optimal Number of CCEP Stimulation Trials"
%    by R Colbur, H Huang, NM Gregg, BH Brinkmann, BN Lundstrom, GA Worrell,
%    G Ojeda Valencia, D Hermes, and KJ Miller.
%
%    Copyright (C) 2022  Ryan Colburn
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
%

%% Load Test Dataset
clearvars

load([name of data file]);

fig_path = [put pathname to save figures to];

subject = 1; % select subject number

setCur = '6.0 mA'; % set stim current to include in analysis

setCur2 = '4.0 mA'; % set alternate current to include in analysis

tt = dataCcepsAllSubs(1).tt'; % get timepoints

%% Subject Sorting and Raw CCEP Visualization
data_struct = dataCcepsAllSubs(subject); % get struct for sub

V = data_struct.V; % get data for sub
V_zeroed = cell(size(V));
sig = cell(size(V));

events = data_struct.events; % get events file for sub
stim_sites = cell(size(V)); % stim sites for sub
stim_trials = cell(size(V)); % stim sites counts for sub

% Get processed signal, list of stimulation sites, information about stimulation trials, and plot & save raw CCEPs for all recording electrodes
for jj = 1:length(events)
    stim_sites{jj} = events{jj}.electrical_stimulation_site;
    unique_cur = unique(events{jj}.electrical_stimulation_current); % find unique current types
    
    if length(unique_cur) > 1                                    % if there is more than one current used in stimulation
        string = strcmp(events{jj}.electrical_stimulation_current, setCur);
        
        if sum(string) >= length(string)/2 % use primary current specified in setCur variable if it makes up half or more of the total trials
            stim_sites{jj} = stim_sites{jj}(string);
            V{jj} = V{jj}(string,:);
            tot = length(string);
            used = sum(string);
            disp(['More than One Current Present in ' data_struct.sub ' ' data_struct.chs{jj} ', Using Only ' setCur ' Trials (' num2str(used) ' Out of ' num2str(tot) ' Total Trials)'])
        else % use the secondary current specified in setCur2 variable if the primary current makes up less than half of the total trials
            string = strcmp(events{jj}.electrical_stimulation_current, setCur2);
            stim_sites{jj} = stim_sites{jj}(string);
            V{jj} = V{jj}(string,:);
            tot = length(string);
            used = sum(string);
            disp(['More than One Current Present in ' data_struct.sub ' ' data_struct.chs{jj} ', Using Only ' setCur2 ' Trials (' num2str(used) ' Out of ' num2str(tot) ' Total Trials)'])
        end

    end

    stim_site_names = unique(stim_sites{jj},'stable'); % find unique stim sites
    stim_trials{jj} = cell(length(stim_site_names),2);
    resp_avg = zeros(length(V{jj}),length(stim_site_names));

    figure; hold on;
    for ii = 1:length(stim_site_names)
        string = strcmp(stim_site_names{ii}, stim_sites{jj}); % find indices for each stim site
        stim_trials{jj}{ii,1} = stim_site_names{ii};
        stim_trials{jj}{ii,2} = sum(string);
        resp_avg(:,ii)  = mean(V{jj}(string,:)); % find average response for each stim site
    end

    plotTrials(tt, resp_avg, 200, stim_site_names);
    xlabel('Time (s)')
    ylabel('Response (Î¼V)')
    title(['Average ' data_struct.sub ' ' data_struct.chs{jj} ' CCEPs for Each Stim Site'])
    hold off

    saveas(gcf, [fig_path data_struct.sub '_' data_struct.chs{jj} '_RawCCEPs.fig'])

    V_zeroed{jj} = (V{jj} - median(V{jj}(:,tt > -0.5 & tt < -0.05),2))'; % baseline correction

    sig{jj} = V_zeroed{jj}(tt > 0.02 & tt < 1,:); % get post-stim signal (20 to 1000 ms)

end

%% Run Correlation Confidence Interval Function
sim = 1000; % specify number of simulations to use for bootstrapping

avg = cell(size(sig));
confInt = cell(size(sig));
r = cell(size(sig));
stim_trials_fun = cell(size(sig));
nosignif = cell(size(sig));
signif = cell(size(sig));
percent_signif = cell(size(sig));

 
% run function for all recording sites using 3-10 input trials and plot & 
% save graph of percentage of significant trials vs # of input trials for 
% each recording site
for ii = 1:length(stim_trials)
    maxTrials = 10;

    avg{ii} = cell(maxTrials-2,1);
    confInt{ii} = cell(maxTrials-2,1);
    r{ii} = cell(maxTrials-2,1);
    stim_trials_fun{ii} = cell(maxTrials-2,1);

    delTri = cell(maxTrials-2,1);

    nosignif{ii} = zeros(maxTrials-2,1);
    signif{ii} = zeros(maxTrials-2,1);
    percent_signif{ii} = zeros(maxTrials-2,1);

    for numTrials = 3:maxTrials % number of trials to use
       
        % delete stimulation sites that have less than number of trials specified in numTrials, not usable for this method
        low_counts = cell2mat(stim_trials{ii}(:,2)) < numTrials; % find indices for low trials
        low_triNames = stim_trials{ii}(low_counts,1); % find names for stim sites with low trials

        if numTrials == 3
            delTri{1} = low_triNames;
        else

            unique_counts = cell2mat(stim_trials{ii}(:,2)) == numTrials-1;
            unique_triNames = stim_trials{ii}(unique_counts,1); % find names for stim sites with number of trials one less than numTrials
            delTri{numTrials-2} = unique_triNames;
        end

        indexes = zeros(length(stim_sites{ii}),length(low_triNames)); % make indexes matrix for indexing for each stim site with low trials
        for jj = 1:length(low_triNames)
            indexes(:,jj) = strcmp(stim_sites{ii},low_triNames{jj}); % find index for particular stim site with low trials
        end

        index = logical(sum(indexes,2)); % make one index matrix with "1" for all stim sites with low trials
        stim_sites{ii}(index) = []; % delete stim sites
        sig{ii}(:,index) = []; % delete signals from those stim sites
        
        % Run function to get confidence intervals
        [avg{ii}{numTrials-2}, confInt{ii}{numTrials-2}, r{ii}{numTrials-2}, stim_trials_fun{ii}{numTrials-2}] = estimateCorConf(sig{ii}, stim_sites{ii}, sim, numTrials);
        
        % Find significant trials from confidence intervals
        signif_log = confInt{ii}{numTrials-2}(:,1) > 0 & confInt{ii}{numTrials-2}(:,2) > 0;
        nosignif_log = signif_log == 0;
        nosignif{ii}(numTrials-2) = sum(nosignif_log(:,1));
        signif{ii}(numTrials-2) = sum(signif_log(:,1));
        percent_signif{ii}(numTrials-2) = sum(signif_log(:,1))/length(signif_log(:,1));

    end

    figure; hold on
    plot(3:maxTrials, 100*percent_signif{ii})
    text(3:maxTrials, 100*percent_signif{ii}-15, delTri,'FontSize',7)
    xlim([3, maxTrials])
    ylim([0, 100])
    xlabel('Number of Trials')
    ylabel('Percentage of Significant Stim Sites')
    title([data_struct.sub ' ' data_struct.chs{ii} ' Percent Significant Trials'])
    hold off
    
    saveas(gcf, [fig_path data_struct.sub '_' data_struct.chs{ii} '_PerSig.jpg'])
 
end