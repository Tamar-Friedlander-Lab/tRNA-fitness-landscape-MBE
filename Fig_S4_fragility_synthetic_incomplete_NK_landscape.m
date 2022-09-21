% define a correlated FLS using the NK model where each genotype is specified by a binary vector
% of length L.
% fitness values are randomly drawn from a uniform distribution.
% We then calculate fragility values using full single-mutant data.

% later, we dilute the landscape in a similar manner to the experimental data and then re-calculate the fragility measure on the diluted landscape.
% Try multiple different dilutions and repeat the fragility calculation to verify that it is unbiased. 

% Written by Tamar Friedlander, 16/06/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
save_FLAG = 0; % to save the data to files change this to 1.
N = 14; % vecotr length
K_vals = 6; %[6, 10, 14]; % epistasis (ruggedness) parameter of the NK model. The higher is K the more rugged (epistatic)the FLS is
Genotypes_dec = 0:2^N - 1;
for ii=1:length(Genotypes_dec)
    Genotypes(ii,:) = dec2binvec(Genotypes_dec(ii), N);
end
 

for k_ind = 1:length(K_vals)
    K = K_vals(k_ind);
    % create all binary k-tuples and define their fitness values
    K_fitness = rand(N,2^K)/ N ; % make the fitness in the range [0,1]
    Fitness = zeros(size(Genotypes_dec));
    for jj = 1:length(Genotypes_dec)
        for nn=1:N % position in the vector
            if nn + K - 1 <= N
                K_tuple = binvec2dec(Genotypes(jj, nn:nn+K-1));
            elseif nn + K - 1 > N
                K_tuple = binvec2dec([Genotypes(jj, nn:N), Genotypes(jj, 1 : mod(nn+K-1,N))]);
            end;
        Fitness(jj) = Fitness(jj) + K_fitness(nn, K_tuple + 1);
        end
    end;
    %scale fitness values to be between [0 1]
    max_F = max(Fitness);
    min_F = min(Fitness);
    Fitness_scaled = (Fitness - min_F) / (max_F - min_F);
    Fitness = Fitness_scaled; 


    % fragility calculation - full data
    Fragility = zeros(size(Genotypes,1),1);
    Neighbors_1 = zeros(N,N);
    Fitness_N1 = zeros(size(Genotypes,1),N);
 
    % calculate fragility (with full data) for each genotype
    N1_dec_all = zeros(size(Genotypes,1), N); % list all 1-neighbors of all genotypes
    for kk=1:size(Genotypes,1) 
        Neighbors_1 = repmat(Genotypes(kk,:), N, 1);
        for nn = 1 : N
            % for every genotype nn, list all its single mutants ('Neighbors_N1')
            % and their fitness values 'Fitness_N1'. 
            Neighbors_1(nn,nn) = ~Neighbors_1(nn,nn); % invert 1 bit to get the nearest neighbor
            N1_dec_all(kk,nn) = binvec2dec(Neighbors_1(nn,:)) + 1; % for every genotype list its 1-mutants (decimal format). Each row stands for the mutants of one genotype. 
            Fitness_N1(kk,nn) = Fitness(binvec2dec(Neighbors_1(nn,:)) + 1); % list of fitness values of all 1-mutants. Every row is the mutants of a single genotype. 
            % we add one because the binary number starts from zero
        end
        Fragility(kk) = mean(max(0, Fitness(kk) - Fitness_N1(kk,:)));
     end

    figure; figure_set_with_title;
    plot(Fitness, Fragility, '.'); xlim([0 1]); hold on;
    xlabel('Fitness');
    ylabel('Mutational fragility');
    title(['K = ',num2str(K-1)]);
    hold on;
    grid on;


    % sample the landscape around a focal genotype
    [~,fittest] = max(Fitness);
    % create also a Focal genotype which is 'wild-type'-like: maximal fitness
    % amongst the low fragility ones. 
    [ii,jj] = find(Fragility==0); % find the fittest amongst the least fragile genotypes, similar to the WT. 
    [mm,ind] = max(Fitness(ii));
    WT_like = ii(ind);
    Focals = [fittest, WT_like];
    plot(Fitness(fittest), Fragility(fittest),'.', 'MarkerSize',20, 'Color','r'); hold on;
    plot(Fitness(WT_like), Fragility(WT_like),'.', 'MarkerSize',20, 'Color','k');
    

    for ff = 1:length(Focals)
        Focal_genotype = Genotypes(Focals(ff),:);
        Distances_from_Focal = sum(repmat(Focal_genotype, 2^N, 1) ~= Genotypes, 2); % distances from the focal of all genotypes (also those not included in the incomplete FL) 
        %Dilution = ones(size(Genotypes));
        p = 0.14;
        Diluted_Num = 1000;
        Fragility_incomplete = cell(1,2^N);
        min_1n_num = 4; % minimal number of 1-neighbors required to calculate fragility
        repeats = 100;
        for mm = 1:repeats  % try different random samplings of the landscape
            rr = rand(Diluted_Num, N) < p;
            mutants = (unique(rr,'rows')); % which posiitons in the focal to mutate
            Num = hist(sum(unique(rr,'rows'),2),0:N); % how many mutations each has, and then count how many mutations with 0, 1, 2, ...N mutations are there.
            Incomplete_FL_bin = xor(repmat(Focal_genotype,sum(Num),1), mutants); % replicate the focal genotype. xor actually applies the mutation scheme of 'mutants' on the focal genotype.
            % the result is a list of genotypes in the incomplete FL - binary representation, each row is a genotype.
            %Distances_from_Focal = sum(repmat(Focal_genotype,sum(Num),1) ~= Incomplete_FL_bin, 2); % distance of each genotype in the incomplete Lanscape from the focal
            Incomplete_FL_dec = zeros(1,size(mutants,1)); 
            for ii = 1:size(mutants,1)
                Incomplete_FL_dec(ii) = binvec2dec(mutants(ii,:)) + 1;
                % list of genotypes in the incomplete FL - decimal representation, each number is a genotype.
            end

        % calculate fragility with incomplete data
            for kk = 1:size(Genotypes, 1)
                if ismember(kk, Incomplete_FL_dec) % genotype kk is included in the incomplete dataset
                    known_N1 = intersect(N1_dec_all(kk,:), Incomplete_FL_dec); % check which of its neighbors are included
                    D = Distances_from_Focal(kk); % distance of genotypes kk from the focal
                    N1_dist = Distances_from_Focal(known_N1);
                    % every genotype has N nearest neighbors, of which D are closer
                    % to the focal at distance D-1 (invert one of the bits separating this genotype
                    % from the focal) and N-D are further away at distance D+1. 
                    % Now distinguish closer and further 1-neighbors and weigh
                    % accordingly to calculate fragility. 
                    Wc = D / N; % relative weight of closer 1-neighbors
                    Wf = (N-D) / N;  % relative weight of further 1-neighbors
                    ClosertoFocal_N1_ind = N1_dist < D;
                    FurtherfromFocal_N1_ind = N1_dist > D;
                    if sum(FurtherfromFocal_N1_ind > 0) & sum(ClosertoFocal_N1_ind > 0) & length(known_N1) > min_1n_num % Both subsets are non-empty and theere is a total minimal number of nearest neighbors
                        Frag = mean(max(0, Fitness(kk) - Fitness(known_N1(ClosertoFocal_N1_ind)) )) * Wc +...
                                mean(max(0, Fitness(kk) - Fitness(known_N1(FurtherfromFocal_N1_ind))))*Wf;
                        Fragility_incomplete{kk} = [Fragility_incomplete{kk}, Frag];
                    elseif D == 0 % it is the focal genotype and hence the 'closer' set is necessarily empty.
                        Frag = mean(max(0, Fitness(kk) - Fitness(known_N1(FurtherfromFocal_N1_ind))))*Wf;
                        Fragility_incomplete{kk} = [Fragility_incomplete{kk}, Frag];
                    end
                end % if ismember(kk,...
            end   % for kk=1:size(Genotypes...
        end % for mm

        [nrows,ncols] = cellfun(@size,Fragility_incomplete);
        figure; figure_set_with_title;
        errorbar(Fragility(ncols>9), cellfun(@mean, Fragility_incomplete(ncols>9)), cellfun(@std, Fragility_incomplete(ncols>9)), 's','MarkerSize',6,...
            'MarkerEdgeColor','red','MarkerFaceColor','red'); hold on; 
        R = corrcoef(Fragility(ncols>9), cellfun(@mean, Fragility_incomplete(ncols>9))); 
        V = axis; 
        text(0.7*(V(2)-V(1)) + V(1), V(3) + 0.2*(V(4) - V(3)), ['\rho = ',num2str(R(1,2))],'FontSize',14); 
        H = refline(1,0);
        set(H, 'LineStyle','--','Color','k','Linewidth',2);
        xlabel('Fragility (full data)');
        ylabel('Fragility (partial data)');
        grid on;
        if ff == 1
            title(['sample around fittest, K = ',num2str(K - 1)]);
        elseif ff==2
            title(['sample around WT-like, K = ',num2str(K - 1)]);
        end

    end % ff = 1:length(Focals)

    if save_FLAG
        fname = ['fragility_incomp_NK_',num2str(K-1),'_v7'];
        save(fname, '-v7','Fitness', 'Fragility','K', 'Fragility_incomplete');
    end
    % to read in python it is needed to save in older Matlab version.

end % for k_ind



%%
% how many different genotypes are sampled?
mutant_num = zeros(1,repeats); 
Neigh_Num = zeros(repeats, N+1);
 for mm = 1:repeats  % try different random samplings of the landscape
            rr = rand(Diluted_Num, N) < p;
            mutants = (unique(rr,'rows')); % which positions in the focal to mutate
            mutant_num(mm) = size(mutants,1);
            Neigh_Num(mm,:) = hist(sum(unique(rr,'rows'),2),0:N);
 end

mean(mutant_num)
std(mutant_num)
mean(Neigh_Num)
std(Neigh_Num)

