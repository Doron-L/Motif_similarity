function PWM = make_PWM_realization(PWM,flag)

% For checking the significance of binding to a specific motif, we randomly exchange 
% the A-T and C-G positions (to preserve the GC content) in each column of the motif's PWM (position
% weight matrix). We can also shuffle between the PWM columns, if it is not important to preserve the 
% information profile.
% flag = 1 and 2 are for don't (default) and do shuffle. 

% We assume the PWM rows are ordered in the following order: A C G T/U.

if nargin == 1, flag = 1; end

% flag = 1
% For each column, we have two random decisions: 1) to make an A-T exchange
% and 2) to make C-G excahnge.
N = size(PWM,2);
exchange = randi([0 1], 2, N); % 0 and 1 means don't and do exchange.

for i = 1:N
   ind = [1 2 3 4];
   if exchange(1,i) == 1, ind(1) = 4; ind(4) = 1; end; if exchange(2,i) == 1, ind(2) = 3; ind(3) = 2; end 
   PWM(:,i) = PWM(ind,i);
end

% flag = 2, add shuffling
if flag == 2, PWM = PWM(:,randperm(N)); end


