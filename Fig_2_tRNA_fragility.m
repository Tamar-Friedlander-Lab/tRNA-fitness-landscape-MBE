% Calculate fragility of the tRNA data
% We weigh the different neighbors as follows: Divide each
% genotype's 1-neighbors into 3 subsets: closer to WT (-1), same distance
% to WT (0), further from WT (+1). Calculate for each genotype the total
% number of genotypes from each subset. Then weigh the actual measured ones
% accordingly.
% No random drawing of 1-neighbors. Simply use all available data. 

% Written by Tamar Friedlander
% 15/4/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% download the data file "GSE111508_FitnessDataMultiEnv.txt.gz" from:
% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111508,
% save, unzip and then save
% as 'mat' file. 
% in the following use 'all_trna_data.mat' data file.
%----------------------------------------------------------------------

clearvars; 
load('all_trna_data.mat');

min_N1_Num = 5;
[~, N1_5_ind] = find(N1_num > min_N1_Num - 1); % genotypes with at least 5 neighbors. 
[~, N1_10_ind] = find(N1_num > 9); % genotypes with at least 5 neighbors. 
OnlyDel_Weighted = zeros(1,length(N1_5_ind));


% calculate correct fragility measures for the non-random choice of closer
% to WT 1-neighbors. 
max_total = 207; 
Closer_Num = zeros(size(N1_5_ind)); 
Eq_Num = zeros(size(N1_5_ind)); 
Further_Num = zeros(size(N1_5_ind)); 

for ii = 1:length(N1_5_ind)
    N1 = N1_indices{N1_5_ind(ii)} + 1; % nearest neighbors
    max_Closer = Num(N1_5_ind(ii)); 
    max_eq = 2 * Num(N1_5_ind(ii)); 
    max_further = max_total - max_eq - max_Closer;
    N1_WT_dist = [];
    for mm=1:length(N1)
        N1_WT_dist = [N1_WT_dist, Num(N1(mm))];
    end
    ClosertoWT_N1_ind = N1_WT_dist < Num(ii);
    SameDistWT_N1_ind = N1_WT_dist == Num(ii);
    FurtherfromWT_ind = N1_WT_dist > Num(ii);
    Closer_Num(ii) = sum(ClosertoWT_N1_ind); 
    Eq_Num(ii) = sum(SameDistWT_N1_ind); 
    Further_Num(ii) = sum(FurtherfromWT_ind); 
    Wc = double(max_Closer) / max_total;
    We = double(max_eq) / max_total;
    Wf = double(max_further) / max_total;
    if Wc+We+Wf ~=1
        error('weights don''t sum to 1.');
    end
    % take care of empty sets - mean of an empty set is NaN!!! specifically
    % the WT gets NaN
     % the full expressions - only work if all 3 sets are non-empty
     % In this version I adjust the weights if one or two sets are empty. 
    if sum (FurtherfromWT_ind > 0) 
        if sum(SameDistWT_N1_ind > 0) & sum(ClosertoWT_N1_ind > 0) % all 3 subsets are non-empty
            OnlyDel_Weighted(ii) = mean(max(0, Fit30(N1_5_ind(ii)) - Fit30(N1(ClosertoWT_N1_ind))))*Wc +...
            mean(max(0, Fit30(N1_5_ind(ii)) - Fit30(N1(SameDistWT_N1_ind))))*We +...
            mean(max(0, Fit30(N1_5_ind(ii)) - Fit30(N1(FurtherfromWT_ind))))*Wf;
        elseif sum(SameDistWT_N1_ind == 0) & sum(ClosertoWT_N1_ind > 0)
             OnlyDel_Weighted(ii) = mean(max(0, Fit30(N1_5_ind(ii)) - Fit30(N1(ClosertoWT_N1_ind))))*(Wc+We) +...
                                    mean(max(0, Fit30(N1_5_ind(ii)) - Fit30(N1(FurtherfromWT_ind))))*Wf;
            elseif sum(SameDistWT_N1_ind > 0) & sum(ClosertoWT_N1_ind == 0)
             OnlyDel_Weighted(ii) = ...
            mean(max(0, Fit30(N1_5_ind(ii)) - Fit30(N1(SameDistWT_N1_ind))))*(We+Wc) +...
            mean(max(0, Fit30(N1_5_ind(ii)) - Fit30(N1(FurtherfromWT_ind))))*Wf;
            
        elseif sum(SameDistWT_N1_ind == 0) & sum(ClosertoWT_N1_ind == 0) % only further set is non-empty
            OnlyDel_Weighted(ii) = mean(max(0, Fit30(N1_5_ind(ii)) - Fit30(N1(FurtherfromWT_ind))))*Wf;

        end
    end

       
        
end

% Plot only genotypes for which the number of single mutant that are
% further from the WT is larger of equal to the number of closer and
% equi-distant ones. 
figure; figure_set;
plot(Fit30(N1_5_ind(Further_Num>=Eq_Num+Closer_Num & Num(N1_5_ind)==1)), OnlyDel_Weighted(Further_Num>=Eq_Num+Closer_Num & Num(N1_5_ind)==1), '.b','MarkerSize',12); hold on; 
plot(Fit30(N1_5_ind(Further_Num>=Eq_Num+Closer_Num & Num(N1_5_ind)==2)), OnlyDel_Weighted(Further_Num>=Eq_Num+Closer_Num & Num(N1_5_ind)==2), '.','Color',[0.3 0.75 0.93],'MarkerSize',12); hold on; 
plot(Fit30(N1_5_ind(Further_Num>=Eq_Num+Closer_Num & Num(N1_5_ind)>=3)), OnlyDel_Weighted(Further_Num>=Eq_Num+Closer_Num & Num(N1_5_ind)>=3), '.','Color',[0.65 0.65 0.65], 'MarkerSize',12); hold on; 
plot(Fit30(1), OnlyDel_Weighted(1),'k.','MarkerSize',36); grid on;
xlabel('Fitness'); ylabel('Mutational fragility');
legend('N1', 'N2','N3 and further', 'WT','Location', 'NorthWest');


% Plot only genotypes that have at least 3 single mutant that are further from the WT
figure; figure_set;
plot(Fit30(N1_5_ind(Further_Num>=3 & Num(N1_5_ind)==1)), OnlyDel_Weighted(Further_Num>=3 & Num(N1_5_ind)==1), '.b','MarkerSize',12); hold on; 
plot(Fit30(N1_5_ind(Further_Num>=3 & Num(N1_5_ind)==2)), OnlyDel_Weighted(Further_Num>=3 & Num(N1_5_ind)==2), '.','Color',[0.3 0.75 0.93],'MarkerSize',12); hold on; 
plot(Fit30(N1_5_ind(Further_Num>=3 & Num(N1_5_ind)>=3)), OnlyDel_Weighted(Further_Num>=3 & Num(N1_5_ind)>=3), '.','Color',[0.65 0.65 0.65], 'MarkerSize',12); hold on; 
plot(Fit30(1), OnlyDel_Weighted(1),'k.','MarkerSize',36); grid on;
xlabel('Fitness'); ylabel('Mutational fragility');
legend('N1', 'N2','N3 and further', 'WT','Location', 'NorthWest');



fitness_N1 = Fit30(N1_5_ind(Further_Num>=3 & Num(N1_5_ind)==1)); fragility_N1 = OnlyDel_Weighted(Further_Num>=3 & Num(N1_5_ind)==1);
fitness_N2 = Fit30(N1_5_ind(Further_Num>=3 & Num(N1_5_ind)==2)); fragility_N2 = OnlyDel_Weighted(Further_Num>=3 & Num(N1_5_ind)==2);
fitness_N3 = Fit30(N1_5_ind(Further_Num>=3 & Num(N1_5_ind)==3)); fragility_N3 = OnlyDel_Weighted(Further_Num>=3 & Num(N1_5_ind)==3);
fitness_WT = Fit30(1); fragility_WT = OnlyDel_Weighted(1); 

save('fragility_fig_data_v7','-v7','fitness_N1','fitness_N2','fitness_N3','fragility_N1','fragility_N2','fragility_N3','fitness_WT','fragility_WT');




