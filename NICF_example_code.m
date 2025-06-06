%%Example code for implementing the NICF


acals_R=reshape(permute(reshape(h,[Rn,%No. of participants , No. of task parameters etc.]),[2,1,3]),[%No. of participants,Rn,No. of task parameters etc.]);

acals_S=reshape(permute(reshape(h,[Sn,%No. of participants , No. of task parameters etc.]),[2,1,3]),[%No. of participants,Sn,No. of task parameters etc.]);

acals_UYZ=reshape(permute(reshape(h,[UYZn,%No. of participants , No. of task parameters etc.]),[2,1,3]),[%No. of participants,UYZn,No. of task parameters etc.]);

acals=cat(2,acals_R,acals_S,acals_UYZ);

permutations=50;p_value=0.01;
[Sf,MI_Fy,MI_FF,CMI_FFy] = mRMR_SFBS(acals,motor_assessment,permutations,p_value,[]);

A={};sigma=2;
combos=[nchoosek(1:Sf,2);[1:length(Sf);1:length(Sf)]'];%inclusive of self-similiarity
for i=1:length(combos)
    a=RBF(acals(:,combos(i,1)),acals(:,combos(i,2)),sigma); %   RBF Kernel
    a=a - diag(diag(a)); %Remove Diagonal
    [threshold] = modified_percolation_analysis(a);a(a<threshold)=0; %Sparsify
    A=cat(1,A,a);
end

A_=cat(3,A{:});
save('A_.mat',"A_"); %Move to python script provided for Leidin community detection or below for Louvain algorithm.


[Opt_rank,M,Vs]=Divisive_Louvain(A);% Cluster structure at first level using Louvain algorithm.

[C] = cluster_certainty(Vs, M); %Cluster affiliation certainties

max_depth=8;current_depth=1;
result_tree = recursive_Divisive_Clustering(A, max_depth, current_depth);% Multiscale cluster structur using Louvain algorithm.


