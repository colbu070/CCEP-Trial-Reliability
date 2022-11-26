%% This function calculates the correlation of input time-series signal and provides confidence intervals to indicate significance against r <= 0.
%
%  
%       [avg, confInt, r, stim_trials] = estimateCorConf(signal, trials, sim, numTrials);
%       signal =          mxn num, signal length for trials by total number of trials
%       trials =          nx1 cell, stimulation sites for each trial in signal
%                         Example may look like [LA1-LA2 LA1-LA2 LA1-LA2 LA2-LA3 LA2-LA3 LA3-LA4 ...]
%       sim =             positive integer value, number of simulations to run in bootstrap
%       numTrials =       positive integer value > 2, number of trials to correlate for each stimulation site
%
%   Returns:
%       avg =             kx1 double, mean Pearson's correlation value for k stimulation sites
%       confInt =         kx2 double, confidence interval for mean for k stimulation sites
%       r =               kxt double, pairwise correlation values for all combinations t for the number of trials specified in numTrials for k stimulation sites
%       stim_trials =     kx2 cell, names of k stimulation sites in first column, number of trials for each site in second
%
%   Based on noise ceiling definition from:
%   Allen, E.J., St-Yves, G., Wu, Y. et al. A massive 7T fMRI dataset to bridge cognitive neuroscience and artificial intelligence. Nat Neurosci 25, 116–126 (2022). https://doi.org/10.1038/s41593-021-00962-x
%
%   estimateCorConf.m
%   2022/11/25
%   Updated 2022/06/27
%
%    If this code is used in a publication, please cite the manuscript:
%    "A Quantitative Method to Determine Optimal Number of CCEP Stimulation Trials"
%    by R Colburn, H Huang, NM Gregg, BH Brinkmann, BN Lundstrom, GA Worrell,
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

function [avg, confInt, r, stim_trials] = estimateCorConf(signal, trials, sim, numTrials)

    stim_site_names = unique(trials, 'stable'); % get unique stimulation sites
    
    % create cell array with stim sites in first column and number of trials for each stim site in second column 
    stim_trials = cell(length(stim_site_names),2);

    for ii = 1:length(stim_site_names)
        string = strcmp(stim_site_names{ii}, trials); % find indices for each stim site
        stim_trials{ii,1} = stim_site_names{ii};
        stim_trials{ii,2} = sum(string);
    end


    % Get r values and put in "r" matrix
     numPairs = nchoosek(numTrials,2); % find n choose 2 to get number of all possible trial pairs
     r = zeros(length(stim_site_names),numPairs); 

    for ii = 1:length(stim_site_names) % for 1 to the number of stim sites
        string = strcmp(stim_site_names{ii},trials); % find specific stim site index
        sig_1stim = signal(:,string); % index for that stim site
        randomTri = randsample(size(sig_1stim,2),numTrials); % randomly choose "numTrials" trials out of total number available
        sig_1stim = sig_1stim(:,randomTri);
        comb = nchoosek(1:numTrials,2); % get indices for all trial pairs
        cor = zeros(numPairs,1); % create empty nx1 double matrix where n is the number of trial pairs
        
        for jj = 1:numPairs
            covar_mat = cov(sig_1stim(:,comb(jj,1)),sig_1stim(:,comb(jj,2))); % find correlations between paired trials
            covar = covar_mat(1,2);
            denominator = sqrt((std(sig_1stim(:,comb(jj,1)))^2)*(std(sig_1stim(:,comb(jj,2)))^2)); % get correlation coefficient
            cor_coef = covar/denominator;
            cor(jj) = cor_coef;
        end
        r(ii,:) = cor;
    end
    
    % Find mean of r values
    avg = zeros(size(r,1),1);
    
    for ii = 1:size(r,1)
        avg(ii) = mean(r(ii,:));
    end
    

    % Use boostrap function to get right tailed 95% confidence intervals
    confInt = zeros(size(r,1),2);

    for ii = 1:size(confInt,1)
        bootstat = bootstrp(sim, @mean, r(ii,:));
        confInt(ii,:) = [prctile(bootstat, 0), prctile(bootstat(:, end:-1:1), 95)];
    end
    
end
