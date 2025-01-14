function C = cluster_certainty(W, Ci)
%CLUSTER_CERTAINTY Computes the probability of cluster affiliation using softmax
%   based on both within-module and without-module degrees.
%
%   C = cluster_certainty(W, Ci);
%
%   Computes a certainty score (probability) for cluster affiliation based on the 
%   within-module degree and without-module degree using softmax normalisation.
%
%   Inputs:
%       W       - connection matrix (binary/weighted, undirected)
%       Ci      - community affiliation vector
%
%   Output:
%       C       - probability of community affiliation for each node.
%

W=weight_conversion(W,"normalize");
n = length(W);  % Number of nodes
C = zeros(n, 1); % Initialize probability array

% Iterate through each community
for i = 1:max(Ci)
    % Compute within-community degree (degree within the community)
    K_within = sum(W(Ci == i, Ci == i), 2);  % Sum of weights within the community
    
    % Compute without-community degree (degree outside the community)
    K_without = sum(W(Ci == i, Ci ~= i), 2); % Sum of weights outside the community
    
    % Compute total degree (within + without) for each node
    K_total = K_within + K_without;
    
    % Apply softmax over within-community degree and without-community degree
    % Softmax function: exp(x_i) / sum(exp(x)) for each node's degree
    K_exp = exp([K_within, K_without]);  % Exponentiate both within and without degrees
    softmax_vals = K_exp ./ sum(K_exp, 2);  % Softmax over the row
    
    % The first column of softmax_vals corresponds to within-community degree probability
    C(Ci == i) = softmax_vals(:, 1);  % Probability of affiliation to community i
end

end

