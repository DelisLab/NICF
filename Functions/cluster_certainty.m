function [C] = cluster_certainty(Vs, cluster_labels)
   
    %%This approach determines the participant-wise cluster affiliation certainty by applying a SoftMax normalisation to the local and global components of the modular centrality.
    %Input:
            %Vs = The Co-Membership matrix.
            %cluster_labels = The cluster membership vector.
    %Output:
            %C = A [No. of participants x [Intra-community centrality , Inter-Community centrality] matrix of cluster affiliation probabilities

%References 
% Ghalmane Z, Cherifi C, Cherifi H, Hassouni ME. Centrality in complex networks with overlapping community structure. Scientific reports. 2019 Jul 12;9(1):10133.

% Initialize the stability vector with zeros
C_intra = nan(size(Vs, 1), size(Vs,2));
C_inter = nan(size(Vs, 1), size(Vs,2));
% Iterate over all pairs of nodes in the matrix
for i = 1:size(Vs, 1)
    for j = 1:size(Vs, 1)
        if i ~= j
            if cluster_labels(i) == cluster_labels(j)
                % Add the co-membership weight if they belong to the same cluster
                C_intra(i,j) = Vs(i, j);
                C_intra(j,i)= Vs(j, i);
            else
                % % Subtract the co-membership weight if they do not belong to the same cluster
                C_inter(i,j) = Vs(i, j);
                C_inter(j,i)= Vs(j, i);
            end
        end
    end
end


C=[mean(C_intra,2,'omitnan'),mean(C_inter,2,'omitnan')];
% % Apply softmax normalization to convert the stability scores to probabilities
C_exp = exp(C);
C = C_exp ./ sum(C_exp,2);   % Normalize to get probabilities
end

