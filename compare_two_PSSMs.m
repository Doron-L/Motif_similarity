function D = compare_two_PSSMs(PWM1,PWM2,minimum_position_overlap)

% This function estimates the similarity between two motifs. It based on 
% “Measurement of sequence similarity” in Itzkovitz et al., 2006.

% PWM is the position weight matrix, it has 4 row for the 4 nucleotid and n
% columns, where n is for the motif length (in nucleotids number).

N1 = size(PWM1,2); N2 = size(PWM2,2);

if nargin == 2, minimum_position_overlap = 5; end

Nshifts = N1+N2-1-(minimum_position_overlap-1)*2;
ind1_end = N1;
ind2_start = 1; ind2_end = minimum_position_overlap;
D_s = zeros(1,Nshifts);
for s = 1:Nshifts % s for shift
    ind1_start = N1-minimum_position_overlap-(s-1)+1;
    if ind1_start < 1, ind1_start = 1; end
    if N2-(minimum_position_overlap+(s-1)) < 0, ind1_end = ind1_end-1; end
    ind1 = ind1_start:ind1_end;
    
    if N1-(minimum_position_overlap+(s-1)) < 0, ind2_start = ind2_start+1; end
    if ind2_end > N2, ind2_end = N2; end
    ind2 = ind2_start:ind2_end;
    ind2_end = ind2_end + 1;
        
    PWM1_s = PWM1(:,ind1); PWM2_s = PWM2(:,ind2);
    H1_s = calc_H(PWM1_s); H2_s = calc_H(PWM2_s);
    I1_s = 2 - H1_s; I2_s = 2 - H2_s;
    
    d = 1 - (calc_H((PWM1_s+PWM2_s)/2) - (H1_s+H2_s)/2);
    
    D_s(s) = max(sum(I1_s.*I2_s.*d));
end

D = max(D_s);

function H = calc_H(PWM)

% This function calculates the entropy (H) of PWM.

% H = -sum(PWM.*log2(PWM));
H_state = PWM.*log2(PWM);
H_state(PWM == 0) = 0; % As PWM goes to 0 faster than log2(PWM) goes to Inf.
H = -sum(H_state);
