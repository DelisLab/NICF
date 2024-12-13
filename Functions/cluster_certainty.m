function [C] = cluster_certainty(Vs, cluster_labels)

%%%Patient-level cluster affiliation certainties
    
    %%This approach determines the stability of cluster affiliations for
    %%each patient and converts the stability scores into cluster
    %%affiliation probabilities using softmax normalisation.

    %Input:
            %Vs = The Co-Membership matrix.
            %cluster_labels = The cluster membership vector.

    %Output:
            %C = A [No. of patients x No. of clusters] matrix of cluster
                 %affiliation probabilities

% Initialize the stability vector with zeros
C = zeros(size(Vs, 1), length(unique(cluster_labels)));
%Vs=weight_conversion(Vs,"normalize");
% Iterate over all pairs of nodes in the matrix
for i = 1:size(Vs, 1)
    for j = 1:size(Vs, 1)
        if i ~= j
            if cluster_labels(i) == cluster_labels(j)
                % Add the co-membership weight if they belong to the same cluster
                C(i, cluster_labels(i)) = C(i, cluster_labels(i)) + Vs(i, j);
                C(j, cluster_labels(j)) = C(j, cluster_labels(j)) + Vs(j, i);
                % C(i, ~cluster_labels(i)) = C(i, ~cluster_labels(i)) - Vs(i, j);
                % C(j, ~cluster_labels(j)) = C(j, ~cluster_labels(j)) - Vs(j, i);
            else
                % % Subtract the co-membership weight if they do not belong to the same cluster
                C(i, ~cluster_labels(i)) = C(i, ~cluster_labels(i)) + Vs(i, j);
                C(j, ~cluster_labels(j)) = C(j, ~cluster_labels(j)) + Vs(j, i);
                % C(i, cluster_labels(i)) = C(i, cluster_labels(i)) - Vs(i, j);
                % C(j, cluster_labels(j)) = C(j, cluster_labels(j)) - Vs(j, i);
            end
        end
    end
end

% Normalize the stability values by the number of nodes to avoid large values
C = C ./ (size(Vs, 2));

% % Apply softmax normalization to convert the stability scores to probabilities
C_exp = exp(C - max(C,[],2));% Subtract row-wise max values for numerical stability
C = C_exp ./ sum(C_exp,2);   % Normalize to get probabilities
end
