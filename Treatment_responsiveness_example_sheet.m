%%Patient clustering of redundant and synergistic network recruitment
%%patterns from pre- to post-session.

Rn=.... %Number of redundant modules to extract
DataR=...%A [[No. of affected-side muscle interactions + No. of unaffected-side muscle interactions] x [No. of Participants x No. of Sessions]] 2d matrix
[w,h]=opnmf(DataR,Rn); %Projective non-negative matrix factorisation (PNMF)
acals_R=reshape(h,[Rn,%No. of participants , No. of Sessions]);

R_dist_prepost={};
combos=[nchoosek(1:Rn,2)];
for i=1:Rn
    for ii=1:Rn
        A=[acals_R(combos(i,1),:,1).*acals_R(combos(i,2),:,2)']; %Linear kernel
        [threshold] = modified_percolation_analysis(A);A(A<threshold)=0; %Sparsify
        R_dist_prepost=cat(1,R_dist_prepost,A);
    end

end

[Opt_rank,M,Vs]=Divisive_Louvain(R_dist_prepost);%...Divisive clustering of pre-to-post session redundant K


Sn=.... %Number of synergistic modules to extract
DataS=...%A [[No. of affected-side muscle interactions + No. of unaffected-side muscle interactions] x [No. of Participants x No. of Sessions]] 2d matrix
[w,h]=opnmf(DataS,Sn); %Projective non-negative matrix factorisation (PNMF)
acals_S=reshape(h,[Sn,%No. of participants , No. of Sessions]);

S_dist_prepost={};
combos=[nchoosek(1:Sn,2)];
for i=1:Sn
    for ii=1:Sn
        A=[acals_S(combos(i,1),:,1).*acals_S(combos(i,2),:,2)']; %Linear kernel
        [threshold] = modified_percolation_analysis(A);A(A<threshold)=0; %Sparsify
        S_dist_prepost=cat(1,S_dist_prepost,A);
    end

end

[Opt_rank,M,Vs]=Divisive_Louvain(S_dist_prepost);%...Divisive clustering of pre-to-post session synergistic K
