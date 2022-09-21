% download the data from:
% https://github.com/lichuan199010/tRNAMultiEnv/blob/master/Reads%20data%20per%20replicate.xlsx
% https://github.com/lichuan199010/tRNAMultiEnv, file name 'Reads data per replicate.xlsx'
% once this file is saved, read it using the following python code and save
% as 'mat' file:
% import scipy.io as sio
% import os
% Reads_data = pd.read_excel('Reads_data.xlsx')
% sio.savemat(('Reads_data.mat'),  {name: col.values for name, col in Reads_data.items()})

% in the following use 'Reads_data.mat' data file.
%--------------------------------------------------------------------------------------


% variable names:
% Num: number of mutations
% Pos: Position of mutation, note that it counts from the variable region
% Nuc: Nucleotide mutated to
% T0All: starting pool
% T23: 23C
% T30: 30C
% D1: DMSO
% T37: 37C
% C1~C5: replicate 1 to 5.


load('Reads_data.mat');

% calculate fitness following Chuan Li's code
zero_vec = zeros(size(T0All));
Data = [Num',zero_vec',zero_vec',T0All',T23C1',T23C2',T23C3',T23C4',T23C5',T30C1',T30C2',T30C3',T30C4',T30C5',D1',D2',D3',T37C1',T37C2',T37C3'];
size(Data);
% number of generations in 24 hrs at each condition
G23 = 13.666; 
G30 = 11.305; 
GDM = 11.894; 
G37 = 11.506; 

% If we look at the porportion of t0 among all the cells with the correct length of insertion t0
t0All = double(sum(T0All));
t0Wt = double(max(T0All)); 
t1All23 = double(sum(sum(Data(:,5:9))));
t1Wt23 = double(sum(Data(1,5:9)));
ngenadd23 = log2(t1Wt23/t0Wt*t0All/t1All23);

t1All30 = double(sum(sum(Data(:,10:14))));
t1Wt30 = double(sum(Data(1,10:14)));
ngenadd30 = log2(t1Wt30/t0Wt*t0All/t1All30);

t1AllDM = double(sum(sum(Data(:,15:17))));
t1WtDM = double(sum(Data(1,15:17)));
ngenaddDM = log2(t1WtDM/t0Wt*t0All/t1AllDM);

t1All37 = double(sum(sum(Data(:,18:20))));
t1Wt37 = double(sum(Data(1,18:20)));
ngenadd37 = log2(t1Wt37/t0Wt*t0All/t1All37);

% add in these generations after correcting for the fraction of WT cells
ngen23 = G23+ngenadd23;
ngen30 = G30+ngenadd30;
ngenDM = GDM+ngenaddDM;
ngen37 = G37+ngenadd37;

% remove genoytpes that don't have at least 100 reads in all experiments. 
MinReadNum = 50; 
Data50 = prod(Data(:,5:20)>MinReadNum, 2); 


RatioEach  = double(Data(:,5:20))./double(T0All)';
wtRatio = double(Data(1,5:20))/double(T0All(1)); 
RatioEach = RatioEach ./ wtRatio; 

% Calculate Fitness 
% Calculate fitness for all
ngenAll = [repmat(ngen23, 1,5), repmat(ngen30, 1, 5),repmat(ngenDM, 1,3), repmat(ngen37,1,3)];
FitnessEach = double(RatioEach).^(1./ngenAll);
FitnessEach(FitnessEach < 0.5) = 0.5;

% average fitness over biological repeats
Mean23C = mean(FitnessEach(:,1:5),2); 
Mean30C = mean(FitnessEach(:,6:10),2); 
MeanDM  = mean(FitnessEach(:,11:13),2); 
Mean37C = mean(FitnessEach(:,14:16),2); 
% Standardize the raw measurements by subtracting the mean for each genotype
Fit23_standard = FitnessEach(:,1:5) - Mean23C;
Fit30_standard = FitnessEach(:,6:10) - Mean30C;
FitDM_standard = FitnessEach(:,11:13) - MeanDM;
Fit37_standard = FitnessEach(:,14:16) - Mean37C;

sd23C = std(FitnessEach(:,1:5),0,2); 
sd30C = std(FitnessEach(:,6:10),0,2); 
sdDM  = std(FitnessEach(:,11:13),0,2); 
sd37C = std(FitnessEach(:,14:16),0,2); 
% standard error of the mean = std divided by sqrt of # experiments.
ExpNum23C = 5; ExpNum30C = 5; ExpNumDM = 3; ExpNum37C = 3;
sem23C = sd23C/sqrt(ExpNum23C); 
sem30C = sd30C/sqrt(ExpNum30C); 
semDM  = sdDM/sqrt(ExpNumDM); 
sem37C = sd37C/sqrt(ExpNum37C); 


figure;
subplot(2,2,1); semilogy(Mean23C, sd23C./Mean23C,'.'); title('23 C');
subplot(2,2,2); semilogy(Mean30C, sd30C./Mean30C,'.'); title('30 C');
subplot(2,2,3); semilogy(MeanDM, sdDM./MeanDM,'.'); title('DMSO');
subplot(2,2,4); semilogy(Mean37C, sd37C./Mean37C,'.'); title('37 C');

% filter out measurements with too few reads
Ind = 1:32576;
Ind50 = Ind(Data50>0); % take only genotypes for which all measurements in all conditions had at least 50 reads - 2110 genotypes.
IndT0100 = Ind(T0All>100); %23200  - this is the filtering used in the original papers. 
Ind_0_5 = Ind(logical(prod(FitnessEach(:,1:16)>0.5,2)' ));
Ind_filt = Ind50; 
Ind23 = (prod(FitnessEach(:,1:5)>0.5, 2)).* (T0All' > 100);
Ind30 = (prod(FitnessEach(:,6:10)>0.5, 2)).* (T0All' > 100);
IndDM = (prod(FitnessEach(:,11:13)>0.5, 2)).* (T0All' > 100);
Ind37 = (prod(FitnessEach(:,14:16)>0.5, 2)).* (T0All' > 100);

Err23 = Fit23_standard(logical(repmat(Ind23,1,5))); 
Err30 = Fit30_standard(logical(repmat(Ind30,1,5))); 
ErrDM = FitDM_standard(logical(repmat(IndDM,1,3))); 
Err37 = Fit37_standard(logical(repmat(Ind37,1,3))); 


figure; 
subplot(2,2,1); histogram(Fit23_standard(logical(repmat(Ind23,1,5)))); grid on; title('23 C');
subplot(2,2,2); histogram(Fit30_standard(logical(repmat(Ind30,1,5)))); grid on; title('30 C');
subplot(2,2,3); histogram(FitDM_standard(logical(repmat(IndDM,1,3)))); grid on; title('DMSO');
subplot(2,2,4); histogram(Fit37_standard(logical(repmat(Ind37,1,3)))); grid on; title('37 C');

% CDF of deviation of SE
figure; 
subplot(2,2,1); histogram(Fit23_standard(logical(repmat(Ind23,1,5))) / sqrt(5), 'Normalization', 'cdf'); grid on; title('23 C'); hold on;
subplot(2,2,2); histogram(Fit30_standard(logical(repmat(Ind30,1,5))) / sqrt(5), 'Normalization', 'cdf'); grid on; title('30 C');
subplot(2,2,3); histogram(FitDM_standard(logical(repmat(IndDM,1,3))) / sqrt(3), 'Normalization', 'cdf'); grid on; title('DMSO');
subplot(2,2,4); histogram(Fit37_standard(logical(repmat(Ind37,1,3)))/ sqrt(3), 'Normalization', 'cdf'); grid on; title('37 C');

figure; 
figure_set_with_title;
histogram(Fit23_standard(logical(repmat(Ind23,1,5))) / sqrt(5), 'Normalization', 'pdf'); grid on; title('23 C'); ylabel('fitness error pdf'); xlim([-0.15 0.15]); hold on;
[mu23,sigma23] = normfit(Fit23_standard(logical(repmat(Ind23,1,5))) / sqrt(5));
plot(-0.15:0.01:0.15, normpdf(-0.15:0.01:0.15, mu23, sigma23), 'r', 'linewidth',2);


figure; 
figure_set_with_title;
histogram(Fit30_standard(logical(repmat(Ind30,1,5))) / sqrt(5), 'Normalization', 'pdf'); grid on; title('30 C'); ylabel('fitness error pdf');  xlim([-0.15 0.15]); hold on;
[mu30,sigma30] = normfit(Fit30_standard(logical(repmat(Ind30,1,5))) / sqrt(5));
plot(-0.15:0.01:0.15, normpdf(-0.15:0.01:0.15, mu30, sigma30), 'r', 'linewidth',2);


figure; 
figure_set_with_title;
histogram(FitDM_standard(logical(repmat(IndDM,1,3))) / sqrt(3), 'Normalization', 'pdf'); grid on; title('DMSO'); ylabel('fitness error pdf');  xlim([-0.15 0.15]); hold on;
[muDM,sigmaDM] = normfit(FitDM_standard(logical(repmat(IndDM,1,3))) / sqrt(3));
plot(-0.15:0.01:0.15, normpdf(-0.15:0.01:0.15, muDM, sigmaDM), 'r', 'linewidth',2);


figure; 
figure_set_with_title;
histogram(Fit37_standard(logical(repmat(Ind37,1,3)))/ sqrt(3), 'Normalization', 'pdf'); grid on; title('37 C'); ylabel('fitness error pdf'); xlim([-0.15 0.15]); hold on;
[mu37,sigma37] = normfit(Fit37_standard(logical(repmat(Ind37,1,3))) / sqrt(3));
plot(-0.15:0.01:0.15, normpdf(-0.15:0.01:0.15, mu37, sigma37), 'r', 'linewidth',2);


% fitness values of 1- and 2-neighbors
[~,J1]= find(Num == 1);
[~,J2]= find(Num == 2);
figure; 
subplot(2,2,1); histogram(Fit23_standard(logical(repmat(Ind23(2:208),1,5))) / sqrt(5), 'Normalization', 'pdf'); grid on; title('23 C, N1');
subplot(2,2,2); histogram(Fit30_standard(logical(repmat(Ind30(2:208),1,5))) / sqrt(5), 'Normalization', 'pdf'); grid on; title('30 C, N1');
subplot(2,2,3); histogram(FitDM_standard(logical(repmat(IndDM(2:208),1,3))) / sqrt(3), 'Normalization', 'pdf'); grid on; title('DMSO, N1');
subplot(2,2,4); histogram(Fit37_standard(logical(repmat(Ind37(2:208),1,3)))/ sqrt(3), 'Normalization', 'pdf'); grid on; title('37 C, N1');

figure; 
subplot(2,2,1); histogram(Fit23_standard(logical(repmat(Ind23(J2),1,5))) / sqrt(5), 'Normalization', 'pdf'); grid on; title('23 C, N2');
subplot(2,2,2); histogram(Fit30_standard(logical(repmat(Ind30(J2),1,5))) / sqrt(5), 'Normalization', 'pdf'); grid on; title('30 C, N2');
subplot(2,2,3); histogram(FitDM_standard(logical(repmat(IndDM(J2),1,3))) / sqrt(3), 'Normalization', 'pdf'); grid on; title('DMSO, N2');
subplot(2,2,4); histogram(Fit37_standard(logical(repmat(Ind37(J2),1,3)))/ sqrt(3), 'Normalization', 'pdf'); grid on; title('37 C, N2');


std23_filt = std(FitnessEach(Ind_filt,1:5),0,2);
std30_filt = std(FitnessEach(Ind_filt,6:10),0,2);
stdDM_filt = std(FitnessEach(Ind_filt,11:13),0,2);
std37_filt = std(FitnessEach(Ind_filt,14:16),0,2);
figure; 
subplot(2,2,1); plot(mean(FitnessEach(Ind_filt,1:5),2), std23_filt,'.'); grid on; xlabel('mean over repeats'); ylabel('std between repeats'); title('23C');
subplot(2,2,2); plot(mean(FitnessEach(Ind_filt,6:10),2), std30_filt,'.'); grid on; xlabel('mean over repeats'); ylabel('std between repeats'); title('30C');
subplot(2,2,3); plot(mean(FitnessEach(Ind_filt,11:13),2), stdDM_filt,'.'); grid on; xlabel('mean over repeats'); ylabel('std between repeats'); title('DMSO');
subplot(2,2,4); plot(mean(FitnessEach(Ind_filt,14:16),2), std37_filt,'.'); grid on; xlabel('mean over repeats'); ylabel('std between repeats'); title('37C');

figure; 
subplot(2,2,1); plot(mean(FitnessEach(Ind_filt,1:5),2), log(FitnessEach(Ind_filt,1:5)) - log(mean(FitnessEach(Ind_filt,1:5),2)),'.'); grid on; xlabel('mean over repeats'); ylabel('log(x) - log(\mu)'); title('23C');
subplot(2,2,2); plot(mean(FitnessEach(Ind_filt,6:10),2), log(FitnessEach(Ind_filt,6:10)) - log(mean(FitnessEach(Ind_filt,6:10),2)),'.'); grid on; xlabel('mean over repeats'); ylabel('log(x) - log(\mu)'); title('30C');
subplot(2,2,3); plot(mean(FitnessEach(Ind_filt,11:13),2), log(FitnessEach(Ind_filt,11:13)) - log(mean(FitnessEach(Ind_filt,11:13),2)),'.'); grid on; xlabel('mean over repeats'); ylabel('log(x) - log(\mu)'); title('DMSO');
subplot(2,2,4); plot(mean(FitnessEach(Ind_filt,14:16),2), log(FitnessEach(Ind_filt,14:16)) - log(mean(FitnessEach(Ind_filt,14:16),2)),'.'); grid on; xlabel('mean over repeats'); ylabel('log(x) - log(\mu)'); title('37C');

%%
% Fig. S2
% Fitness variation between biological repeats under each condition.

figure; figure_set; semilogy(Mean23C(Ind_filt), sem23C(Ind_filt)./Mean23C(Ind_filt),'.');  grid on; ylim([0.0001 0.1]); xlabel('Mean over biological repeats, 23C'); ylabel('sem/mean');
figure; figure_set; semilogy(Mean30C(Ind_filt), sem30C(Ind_filt)./Mean30C(Ind_filt),'.');  grid on; ylim([0.0001 0.1]); xlabel('Mean over biological repeats, 30C'); ylabel('sem/mean');
figure; figure_set; semilogy(MeanDM(Ind_filt), semDM(Ind_filt)./MeanDM(Ind_filt),'.');  grid on; ylim([0.0001 0.1]); xlabel('Mean over biological repeats, DMSO'); ylabel('sem/mean');
figure; figure_set; semilogy(Mean37C(Ind_filt), sem37C(Ind_filt)./Mean37C(Ind_filt),'.');  grid on; ylim([0.0001 0.1]); xlabel('Mean over biological repeats, 37C'); ylabel('sem/mean');

%%
% create the empirical distribution of errors in the geometric mean
% Fig. S3
IndAll = Ind23.*Ind30.*IndDM.*Ind37;
GeoMeanFit = (Mean23C.*Mean30C.*MeanDM.*Mean37C).^0.25;
[IndAll_filt,~] = find(IndAll>0);

r23 = mean(Err23(randi(length(Err23), length(IndAll_filt), 5)),2); 
r30 = mean(Err30(randi(length(Err30), length(IndAll_filt), 5)),2); 
rDM = mean(ErrDM(randi(length(ErrDM), length(IndAll_filt), 3)),2); 
r37 = mean(Err37(randi(length(Err37), length(IndAll_filt), 3)),2); 
GeoMeanFit_Err = ((Mean23C(IndAll_filt) + r23).*(Mean30C(IndAll_filt) + r30).*(MeanDM(IndAll_filt) + rDM).*(Mean37C(IndAll_filt) + r37)).^0.25 - GeoMeanFit(IndAll_filt); 
figure; figure_set_with_title;
histogram(GeoMeanFit_Err,'Normalization','pdf'); grid on; title('error in geometric mean fitness'); ylabel('pdf'); xlim([-0.06 0.06]); hold on;
[muGeoMean,sigmaGeoMean] = normfit(GeoMeanFit_Err)
plot(-0.06:0.001:0.06, normpdf(-0.06:0.001:0.06, muGeoMean,sigmaGeoMean), 'r', 'linewidth',2);
figure; figure_set_with_title; 
histogram(GeoMeanFit_Err,'Normalization','cdf'); grid on; title('error in geometric mean fitness'); ylabel('cdf'); xlim([-0.06 0.06]);
hold on;
Q = quantile(GeoMeanFit_Err, [1-0.95^(1/3), 1-0.95^(1/2), 0.05, 0.1, 0.2]);
line([Q(1) Q(1)], [0 1]); 
line([Q(2) Q(2)], [0 1]);
line([Q(3) Q(3)], [0 1]);
legend('CDF','0.95 confidence 3-step','0.95 confidence 2-step','0.95 confidence 1-step')
% calculate thresholds for the trajectory count: assume the whole
% trajectory should have confidence of 0.95, then for a 2-step or 3-step
% trajectory take the 2- or 3- root of 0.95. 



%%
% Fig. S1
% estimate read count noise assuming Binomial distribution of the measured number of reads. 
draw_num = 20; % 50 is too much. 20 still worked. For large number of draws It's possible to do it with a loop drawing one a time - takes the same time. 
gene_ind = 5383;
p30_1 = double(T30C1(gene_ind))/double(sum(T30C1));  
p30_2 = double(T30C2(gene_ind))/double(sum(T30C2));
p30_3 = double(T30C3(gene_ind))/double(sum(T30C3));
p30_4 = double(T30C4(gene_ind))/double(sum(T30C4));
p30_5 = double(T30C5(gene_ind))/double(sum(T30C5));
p0 = double(T0All(gene_ind))/double(sum(T0All));
sample30_1 = binornd(sum(T30C1),p30_1,1,draw_num); 
sample30_2 = binornd(sum(T30C2),p30_2,1,draw_num);
sample30_3 = binornd(sum(T30C3),p30_3,1,draw_num);
sample30_4 = binornd(sum(T30C4),p30_4,1,draw_num);
sample30_5 = binornd(sum(T30C5),p30_5,1,draw_num);
sample_t0 =  binornd(double(sum(T0All)),p0,1,draw_num);

fitness_dist_1 = (kron(sample30_1,1./sample_t0)/wtRatio(6)).^(1/ngen30);
fitness_dist_2 = (kron(sample30_2,1./sample_t0)/wtRatio(7)).^(1/ngen30);
fitness_dist_3 = (kron(sample30_3,1./sample_t0)/wtRatio(8)).^(1/ngen30);
fitness_dist_4 = (kron(sample30_4,1./sample_t0)/wtRatio(9)).^(1/ngen30);
fitness_dist_5 = (kron(sample30_5,1./sample_t0)/wtRatio(10)).^(1/ngen30);

figure; figure_set;
subplot(2,3,1); histogram(fitness_dist_1); xlim([1.5 1.65]); grid on;
subplot(2,3,2); histogram(fitness_dist_2); xlim([1.5 1.65]); grid on;
subplot(2,3,3); histogram(fitness_dist_3); xlim([1.5 1.65]); grid on;
subplot(2,3,4); histogram(fitness_dist_4); xlim([1.5 1.65]); grid on;
subplot(2,3,5); histogram(fitness_dist_5); xlim([1.5 1.65]); grid on;
figure; figure_set; 
histogram(fitness_dist_1); hold on;
histogram(fitness_dist_2); hold on;
histogram(fitness_dist_3); hold on;
histogram(fitness_dist_4); hold on;
histogram(fitness_dist_5); hold on;
xlabel('fitness at 30C, genotype #5383');
ylabel('Counts');

