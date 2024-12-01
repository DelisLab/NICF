%%Patient clustering of redundant and synergistic network recruitment
%%patterns at pre- and post-sessions.

Rn=.... %Number of redundant modules to extract
DataR=...%A [[No. of affected-side muscle interactions + No. of unaffected-side muscle interactions] x [No. of Participants x No. of Sessions]] 2d matrix
[w,h]=opnmf(DataR,Rn); %Projective non-negative matrix factorisation (PNMF)
acals_R=reshape(h,[Rn,%No. of participants , No. of Sessions]);

R_dist_pre={};R_dist_post={};
combos=[nchoosek(1:Rn,2)];
for i=1:length(combos)
    A=[acals_R(combos(i,1),:,1).*acals_R(combos(i,2),:,1)']; %Linear kernel
    A=A - diag(diag(A)); %Remove Diagonal
    [threshold] = modified_percolation_analysis(A);A(A<threshold)=0; %Sparsify
    R_dist_pre=cat(1,R_dist_pre,A);


    A=[acals_R(combos(i,1),:,2).*acals_R(combos(i,2),:,2)']; %Linear kernel
    A=A - diag(diag(A)); %Remove Diagonal
    [threshold] = modified_percolation_analysis(A);A(A<threshold)=0;  %Sparsify
    R_dist_post=cat(1,R_dist_post,A);
end

[Opt_rank,M,Vs]=Consensus(R_dist_pre);%...Divisive clustering of pre-session redundant K
[Opt_rank,M,Vs]=Consensus(R_dist_post);%...Divisive clustering of post-session redundant K

Sn=.... %Number of synergistic modules to extract
DataS=...%A [[No. of affected-side muscle interactions + No. of unaffected-side muscle interactions] x [No. of Participants x No. of Sessions]] 2d matrix
[w,h]=opnmf(DataS,Sn); %Projective non-negative matrix factorisation (PNMF)
acals_S=reshape(h,[Sn,%No. of participants , No. of Sessions]);

S_dist_pre={};S_dist_post={};
combos=[nchoosek(1:Sn,2)];
sigma=2;
for i=1:length(combos)
    A=RBF(acals_S(combos(i,1),:,1),acals_S(combos(i,2),:,1),sigma); %RBF kernel
    A=A - diag(diag(A)); %Remove Diagonal
    [threshold] = modified_percolation_analysis(A);A(A<threshold)=0; %Sparsify
    S_dist_pre=cat(1,S_dist_pre,A);


    A=RBF(acals_S(combos(i,1),:,2),acals_S(combos(i,2),:,2),sigma); %RBF kernel
    A=A - diag(diag(A)); %Remove Diagonal
    [threshold] = modified_percolation_analysis(A);A(A<threshold)=0;  %Sparsify
    S_dist_post=cat(1,S_dist_post,A);
end

[Opt_rank,M,Vs]=Consensus(S_dist_pre);%...Divisive clustering of pre-session synergistic K
[Opt_rank,M,Vs]=Consensus(S_dist_post);%...Divisive clustering of post-session synergistic K

