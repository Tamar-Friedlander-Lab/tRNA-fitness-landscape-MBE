% define a correlated FLS using the NK model where each genotype is specified by a binary vector
% of length L.
% fitness values are randomly drawn from a uniform distribution.
% We then calculate fragility values and plot fragility vs. fitness and
% fitness value distributions. 

% Written by Tamar Friedlander
% April 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;

N = 14; % vecotr length
K_vals = [6  14]; % epistasis (ruggedness) parameter of the NK model. The higher is K the more rugged (epistatic)the FLS is
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
        for nn=1:N
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


    figure; figure_set_with_title;
    hist(Fitness,50);
    xlim([0 1]);
    xlabel('fitness values');
    title(['K = ',num2str(K-1)]);
    hold on;

    OnlyDel = zeros(size(Genotypes,1),1);
    Neighbors_1 = zeros(N,N);
    Fitness_N1 = zeros(1,N);

    % calculate steepness at each genotype
    for kk=1:size(Genotypes,1) 
        Neighbors_1 = repmat(Genotypes(kk,:), N, 1);
        for nn = 1 : N
            Neighbors_1(nn,nn) = ~Neighbors_1(nn,nn); % invert 1 bit to get the nearest neighbor
            Fitness_N1(nn) = Fitness(binvec2dec(Neighbors_1(nn,:)) + 1); % locate the neighbor's Fitness value
            % we add one because the binary number starts from zero
        end
        OnlyDel(kk) = mean(max(0, Fitness(kk) - Fitness_N1));
    end

    figure; figure_set_with_title;
    plot(Fitness, OnlyDel, '.'); xlim([0 1]); hold on;
    xlabel('Fitness');
    ylabel('Mutational fragility');
    title(['K = ',num2str(K-1)]);
    hold on;
    grid on;


    fname = ['fragility_NK_',num2str(K-1),'_v7'];
    save(fname, '-v7','Fitness', 'OnlyDel','K');
    % to read in python it is needed to save in older Matlab version.

end % for k_ind
