function result_tree = recursive_Divisive_Clustering(X, max_depth, current_depth)
    % Recursive Divisive Clustering with the Louvain algorithm
    %
    % Arguments:
    % X: 2D cell array of cell arrays, where each subcell represents data for a graph layer.
    % max_depth: Maximum recursion depth.
    % current_depth: Current recursion depth (default is 1).
    %
    % Returns:
    % result_tree: A structure representing the tree of community detection results.

    % Default values for optional arguments
    if nargin < 2
        max_depth = 3;
    end
    if nargin < 3
        current_depth = 0;
    end

    % Stop recursion if the current depth exceeds the max_depth
    if current_depth >= max_depth
        result_tree = [];
        return;
    end

    % Apply the Divisive Louvain algorithm to the 2D cell array X
    [Opt_rank, M, aggregated_vs] = Divisive_Louvain(X);

    % If Opt_rank is 1, stop recursion and return empty (this branch will not have any output)
    if Opt_rank == 1
        result_tree = [];
        return;
    end

    % Prepare the tree structure at the current depth
    result_tree.Depth = current_depth;
    result_tree.CommunityMembership = M;
    result_tree.AggregatedCoMembershipMatrix = aggregated_vs;

    % Create subgraphs based on community membership
    subgraphs = struct();
    N = length(M);
    unique_communities = unique(M);
    
    for i = 1:length(unique_communities)
        community_id = unique_communities(i);

        % Create a subgraph for each community
        community_nodes = find(M == community_id);
        combo_nodes = nchoosek(community_nodes, 2);
        
        if length(unique_communities) > 1 && length(community_nodes) >2% Only recurse on non-trivial subgraphs
            % Initialize a cell array to hold the subgraph data for each layer
            subgraph_data = cell(1, length(X));  % One entry per layer

            for layer = 1:length(X)
                % Get the current layer matrix
                layer_matrix = X{layer};
                A = nan(size(layer_matrix));  % Initialize adjacency matrix A

                % Fill A with edges between pairs of nodes in community_nodes
                for combo = 1:size(combo_nodes, 1)
                    A(combo_nodes(combo, 1), combo_nodes(combo, 2)) = layer_matrix(combo_nodes(combo, 1), combo_nodes(combo, 2));
                    A(combo_nodes(combo, 2), combo_nodes(combo, 1)) = layer_matrix(combo_nodes(combo, 2), combo_nodes(combo, 1));
                end

                % Remove rows and columns that are all zeros
                rows_to_delete = all(isnan(A), 2);
                cols_to_delete = all(isnan(A), 1);
                rows_and_cols_to_delete = rows_to_delete' & cols_to_delete;
                A(rows_and_cols_to_delete, :) = [];
                A(:, rows_and_cols_to_delete) = [];
                A(isnan(A))=0;
                % Store the subgraph data for the current layer
                subgraph_data{layer} = A;
            end

            % Recursive call for the subgraph
            subgraphs.(sprintf('Community_%d', community_id)) = ...
                recursive_Divisive_Clustering(subgraph_data, max_depth, current_depth + 1);

        end
    end

    % Only add subgraphs if there are any
    if ~isempty(fieldnames(subgraphs))
        result_tree.Subgraphs = subgraphs;
    else
        result_tree.Subgraphs = struct();
    end
end