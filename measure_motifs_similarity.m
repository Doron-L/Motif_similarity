function [D,Dscaled,pval,qval,motif1_num,motif2_num] = measure_motifs_similarity(PWM_motif_set1,PWM_motif_set2,varargin)

% In this function we estimate the similarity and its significance (P-value) between two sets of motifs. The 
% similarity is taken to be one minus the Shannon-Jensen distance. More specifically we follow the procedure 
% presented in Itzkovitz et al. 2006, "Coding limits on the number of transcription factors".

% Input:
% -----
% PWM_motif1: is a cell array of Position Weight Matrix (PWM) of the first motif set. Its length is the number of 
% motifs,k. Each cell is an n X m matrix, where the number of rows is the number of different options each element 
% can have (for example, for nucleotides there are 4 options (A, T, C, and G)), and the number of columns is the 
% number of the motif's base number. In each column, we have the probability for each element kind to be in this 
% position. The sum of each column must give 1.
% PWM_motif2: is a cell array of Position Weight Matrix (PWM) of the second
% motif set. Its structure is exactly as the motif1 structure, though its length may be different.

% Output:
% D: the similarity k1 X k2 matrix, where k1 is for the different motifs in set 1 and k2 is for the different motifs 
% in set 2.
% Dscaled: is the scaled similarity between motifs1 and motifs2 by the minimum of their self similarity. 
% pval and qval: are the significance of the values in the D matrix in P-values and q-values, respectively.
% motif1_num,motif2_num: the indexes of the motifs in the first and second sets that where found significantly 
% similar to each other.  

% Example:
% PWM_motif_set1 =  {
% [0.1672  0.0000  0.5615  0.1728  0.2313  0.0000  0.0000  0.0000
% 0.8328  0.0964  0.3550  0.0000  0.0000  0.1510  1.0000  0.0847
% 0.0000  0.5523  0.0835  0.8272  0.0755  0.1512  0.0000  0.9153
% 0.0000  0.3513  0.0000  0.0000  0.6933  0.6978  0.0000  0.0000],
% [0.0000  0.0000  0.7418  0.6715  0.0000  0.0000  0.0828  0.0000
% 1.0000  0.0000  0.0000  0.0000  0.3364  0.2821  0.3586  0.0000
% 0.0000  1.0000  0.2582  0.0967  0.5339  0.0000  0.5586  1.0000
% 0.0000  0.0000  0.0000  0.2318  0.1297  0.7179  0.0000  0.0000]
% };
% 
% PWM_motif_set2 =  {
% [0.1672  0.0000  0.5615  0.1728  0.2313  0.0000  0.0000  0.0000
% 0.8328  0.0964  0.3550  0.0000  0.0000  0.1510  1.0000  0.0847
% 0.0000  0.5523  0.0835  0.8272  0.0755  0.1512  0.0000  0.9153
% 0.0000  0.3513  0.0000  0.0000  0.6933  0.6978  0.0000  0.0000],
% [0.0000  0.0000  0.7418  0.6715  0.0000  0.0000  0.0828  1.0000
% 0.0000  0.5000  0.0000  0.0000  0.3364  0.2821  0.3586  0.0000
% 0.5000  0.0000  0.2582  0.0967  0.5339  0.0000  0.5586  0.0000
% 0.5000  0.5000  0.0000  0.2318  0.1297  0.7179  0.0000  0.0000]
% };
%
% [~,~,~,qval] = measure_motifs_similarity(PWM_motif1,PWM_motif2)

%% Setting default values for the input
% The default input is set for:
% Nsim, Sim_flag, Dscaled_thresh, qvalue_thresh
numvarargs = length(varargin);

% set defaults for optional inputs
optargs = {1e4 1 0.7 0.1};

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[Nsim, Sim_flag, Dscaled_thresh, qvalue_thresh] = optargs{:};
%%

% Finding the number of motifs in each set.
N_PWM_motif_set1 = length(PWM_motif_set1);
N_PWM_motif_set2 = length(PWM_motif_set2);

% Estimating the self similarity of motif1.
Dself_PWM_motif_set1 = zeros(1,N_PWM_motif_set1);
for ind_motif_set1 = 1:N_PWM_motif_set1
    Dself_PWM_motif1(ind_motif_set1) = compare_two_PSSMs(PWM_motif_set1{ind_motif_set1},PWM_motif_set1{ind_motif_set1});
end

% Estimating the self similarity of motif2.
Dself_PWM_motif_set2 = zeros(1,N_PWM_motif_set2);
for ind_motif_set2 = 1:N_PWM_motif_set2
    Dself_PWM_motif_set2(ind_motif_set2) = compare_two_PSSMs(PWM_motif_set2{ind_motif_set2},PWM_motif_set2{ind_motif_set2});
end

D_sim = zeros(1,Nsim);
for ind_motif_set2 = 1:N_PWM_motif_set2
    for ind_motif_set1 = 1:N_PWM_motif_set1
        
        D(ind_motif_set2,ind_motif_set1) = compare_two_PSSMs(PWM_motif_set1{ind_motif_set1},PWM_motif_set2{ind_motif_set2});
        clear D_sim Dscaled_sim
        for ind_sim = 1:Nsim
            PWM_motif_set1_sim = make_PWM_realization(PWM_motif_set1{ind_motif_set1},Sim_flag);
            PWM_motif_set2_sim = make_PWM_realization(PWM_motif_set2{ind_motif_set2},Sim_flag);
            D_sim(ind_sim) = compare_two_PSSMs(PWM_motif_set1_sim,PWM_motif_set2_sim);
            
        end
        % Estimating the Z score
        Z=(D(ind_motif_set2,ind_motif_set1)-mean(D_sim))/std(D_sim);
        % Estimating the P-vlaue
        pval(ind_motif_set2,ind_motif_set1) = 1-normcdf(Z);

    end
end

% Scaling the similarity between motifs1 and motifs2 by the minimum of their self similarity. 
Dscaled = D./min(repmat(Dself_PWM_motif_set2',1,size(D,2)),repmat(Dself_PWM_motif_set1,size(D,1),1));

% Estimating the q-value.
[~,qval] = mafdr(pval(:));

[motif2_num,motif1_num] = ind2sub(size(pval),find(qval < qvalue_thresh & Dscaled(:) > Dscaled_thresh));

qval = reshape(qval,size(pval));
% for ind_RBP_num = 1:length(RBP_num)
%     seqlogo(PWM_RBP{RBP_num(ind_RBP_num)})
% end


