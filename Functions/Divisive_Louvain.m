function [Opt_rank,M,Vs]=Divisive_Louvain(X)
%%%%%%%%%%%%%%%%%%%%
%Divisive clustering using the louvain algorithm for network community
%detection

%Input:
        %X: A cell array containing the multiplex network
                
%Output:
        %Opt_rank: The optimal number of clusters to extract
        %M: Graph partition vector
        %Vs: Co-Membership matrix

%%%%%%%%%%%%%%%%%%%%%


Vs=[];
combos_s=nchoosek(1:size(X{1},1),2);
for i=1:length(X)
    v=[];
    net=X{i};
    [M,Q]=community_louvain(net);
    for ii=1:length(combos_s)
        m1=M(combos_s(ii,1),:);
        m2=M(combos_s(ii,2),:);
        if isequal(m1,m2)
            v=[v;1];
        else
            v=[v;0];
        end
    end
    Vs=cat(3,Vs,squareform(v));
end
Vs=sum(Vs,3);


[M,Q]=community_louvain(Vs);
Opt_rank=max(M);

end