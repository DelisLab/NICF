function [Sf,MI_Fy,MI_FF,CMI_FFy] = mRMR_SFBS(F,y,permutations,p_value,discrete)

% mRMR_SFBS: Sequential forward-backwards selection using Maximum relevance
% minimum redundancy criterion.

% INPUTS:
% F            - Matrix of features (rows = samples, columns = features).
% y            - Target variable (vector).
% permutations - Number of permutations for statistical significance testing.
% p_value      - P-value threshold for significance.
% discrete     - Boolean flag: true if the target variable y is discrete.
%
% OUTPUTS:
% Sf           - Selected features (indices of relevant features).
% MI_Fy        - MI between features and target variable.
% MI_FF        - MI between pairs of features.
% CMI_FFy      - CMI between feature pairs and target variable.
%
%References:
% 1) Wollstadt P, Schmitt S, Wibral M. A rigorous information-theoretic definition of redundancy and relevancy in feature selection 
%    based on (partial) information decomposition. Journal of Machine Learning Research. 2023;24(131):1-44.
% 2) Ince RA, Giordano BL, Kayser C, Rousselet GA, Gross J, Schyns PG. A statistical framework for neuroimaging data analysis based on mutual information 
%    estimated via a gaussian copula. Human brain mapping. 2017 Mar;38(3):1541-73.
% 3) Ince RA, Mazzoni A, Bartels A, Logothetis NK, Panzeri S. A novel test to determine the significance of neural selectivity to single and 
%    multiple potentially correlated stimulus features. Journal of neuroscience methods. 2012 Sep 15;210(1):49-65.
% 4) Brown G, Pocock A, Zhao M, Luján M. Conditional likelihood maximisation: a unifying framework for mutual information feature selection. 
%    J Mach Learn Res. 2011;13:26-46.


combos=nchoosek(1:size(F,2),2);
F=copnorm(F);
if discrete
else
    y=copnorm(y);
end

if discrete
    MI_Fy=[];
    for i=1:size(F,2)
        mi=mi_model_gd(F(:,i), y, max(y)+1,true,true);
        mi_Ip=mi_model_gd(F(:,i), y, max(y)+1,false,true);
        perms=[];
        for iter=1:permutations
            perms=[perms;mi_model_gd(F(randperm(size(F,1)),i),y,max(y)+1,false,true)];
        end
        mu=mean(perms);
        sd=std(perms);
        zscore=norminv(p_value);
        perm=mu+(zscore*sd);

        if mi_Ip<perm || mi<0
            MI_Fy=[MI_Fy;[0,i]];
        else
            MI_Fy=[MI_Fy;[mi,i]];
        end
    end
else
    MI_Fy=[];
    for i=1:size(F,2)
        mi=mi_gg(y,F(:,i),true,true);
        mi_Ip=mi_gg(y,F(:,i),false,true);
        perms=[];
        for iter=1:permutations
            perms=[perms;mi_gg(y,F(randperm(size(F,1)),i),false,true)];
        end
        mu=mean(perms);
        sd=std(perms);
        zscore=norminv(p_value);
        perm=mu+(zscore*sd);
        if mi_Ip<perm || mi<0
            MI_Fy=[MI_Fy;[0,i]];
        else
            MI_Fy=[MI_Fy;[mi,i]];
        end
    end
end

MI_sorted=sortrows(MI_Fy,'descend');
if MI_sorted(1,1)==0
    Sf=NaN;MI_FF=NaN;CMI_FFy=NaN;
    fprintf('No relevant features found, discontinuing....')
    return
else
    Si=[MI_sorted(1,2)];
    MI_FF=[];CMI_FFy=[];
    for i=1:length(combos)
        mi=mi_gg(F(:,combos(i,1)),F(:,combos(i,2)),true,true);
        mi_Ip=mi_gg(F(:,combos(i,1)),F(:,combos(i,2)),false,true);

        if discrete
            [cmi I] = gccmi_ccd(F(:,combos(i,1)),F(:,combos(i,2)),y,max(y)+1);
            [cmi_Ip I] =gccmi_ccd(F(:,combos(i,1)),F(:,combos(i,2)),y,max(y)+1);
        else
            cmi=cmi_ggg(F(:,combos(i,1)),F(:,combos(i,2)),y,true,true);
            cmi_Ip=cmi_ggg(F(:,combos(i,1)),F(:,combos(i,2)),y,false,true);
        end
        perms_mi=[];perms_cmi=[];
        for iter=1:permutations
            if discrete
                perms_mi=[perms_mi;mi_gg(F(randperm(size(F,1)),combos(i,1)),F(:,combos(i,2)),false,true)];
                perms_cmi=[perms_cmi;gccmi_ccd(F(randperm(size(F,1)),combos(i,1)),F(randperm(size(F,1)),combos(i,2)),y,max(y)+1)];
            else
                perms_mi=[perms_mi;mi_gg(F(randperm(size(F,1)),combos(i,1)),F(:,combos(i,2)),false,true)];
                perms_cmi=[perms_cmi;cmi_ggg(F(randperm(size(F,1)),combos(i,1)),F(randperm(size(F,1)),combos(i,2)),y,false,true)];
            end
        end
        mu=mean(perms_mi);
        sd=std(perms);
        zscore=norminv(p_value);
        perm_mi=mu+(zscore*sd);
        if mi_Ip<perm_mi || mi<0
            MI_FF=[MI_FF;[0,combos(i,:)]];
        else
            MI_FF=[MI_FF;[mi,combos(i,:)]];
        end

        mu=mean(perms_cmi);
        sd=std(perms);
        zscore=norminv(p_value);
        perm_cmi=mu+(zscore*sd);
        if cmi_Ip<perm_cmi || cmi<0
            CMI_FFy=[CMI_FFy;[0,combos(i,:)]];
        else
            CMI_FFy=[CMI_FFy;[cmi,combos(i,:)]];
        end
    end


    Xi=1:size(F,2);Xi(Xi==MI_sorted(1,2))=[];


    while length(Xi)>0
        CMI_FS=[];
        for i=1:length(Xi)

            mi_Fi_y=MI_sorted(MI_sorted(:,2)==Xi(i),1);
            for ii=1:length(Si)
                mi_fi_fjs = sum(MI_FF(((MI_FF(:,2)==Xi(i)) | (MI_FF(:,3)==Xi(i))) & ((MI_FF(:,2)==Si(ii)) | (MI_FF(:,3)==Si(ii))), 1));
                cmi_fi_fjs_y = sum(CMI_FFy(((CMI_FFy(:,2)==Xi(i)) | (CMI_FFy(:,3)==Xi(i))) & ((CMI_FFy(:,2)==Si(ii)) | (CMI_FFy(:,3)==Si(ii))), 1));


                CMI_FS=[CMI_FS;[mi_Fi_y-(mi_fi_fjs+cmi_fi_fjs_y),Xi(i),Si(ii)]];
            end
        end
        % Select the feature with the maximum CMI
        [~, max_idx] = max(CMI_FS(:, 1));
        if CMI_FS(max_idx, 1) > 0
            Si = [Si; CMI_FS(max_idx, 2)];
            Xi(Xi == CMI_FS(max_idx, 2)) = []; % Remove selected feature from Xi
        else
            break; % Stop if no more significant features
        end
    end

    Sf = Si;
    while length(Sf) > 1
        disp(['Current length of Sf: ', num2str(length(Sf))]);  % Debugging: print the length of Sf
        feature_removed = false;  % Flag to track if a feature is removed

        for i = length(Si):-1:1
            % Check if the length of Sf is already 1, exit immediately if so
            if length(Sf) <= 1
                break;
            end

            mi_Fi_y = MI_sorted(MI_sorted(:, 2) == Si(i), 1);
            CMI_FS = [];

            for ii = length(Si):-1:1
                if ii ~= i
                    % Calculating the CMI and adding to the list
                    mi_fi_fjs = sum(MI_FF((MI_FF(:, 2) == Si(i) | MI_FF(:, 3) == Si(i)) & (MI_FF(:, 2) == Si(ii) | MI_FF(:, 3) == Si(ii)), 1));
                    cmi_fi_fjs_y = sum(CMI_FFy((CMI_FFy(:, 2) == Si(i) | CMI_FFy(:, 3) == Si(i)) & (CMI_FFy(:, 2) == Si(ii) | CMI_FFy(:, 3) == Si(ii)), 1));

                    % Add results to CMI_FS
                    CMI_FS = [CMI_FS; [mi_Fi_y - (mi_fi_fjs + cmi_fi_fjs_y), Si(i), Si(ii)]];
                end
            end

            % Select the feature with the minimum CMI
            [~, min_idx] = min(CMI_FS(:, 1));

            % Debugging: Check what is happening with CMI_FS
            disp(['CMI_FS(min_idx, 1): ', num2str(CMI_FS(min_idx, 1))]);
            disp(['CMI_FS(min_idx, 2): ', num2str(CMI_FS(min_idx, 2))]);

            % Remove the feature if CMI is zero or less
            if CMI_FS(min_idx, 1) <= 0
                Si(Si == CMI_FS(min_idx, 2)) = [];  % Remove selected feature from Si
                Sf(Sf == CMI_FS(min_idx, 2)) = [];  % Remove selected feature from Sf

                feature_removed = true;  % Set the flag to true
                disp(['Removed feature: ', num2str(CMI_FS(min_idx, 2))]);  % Debugging: print removed feature
            end
        end

        % Break the loop if no features were removed in this pass
        if ~feature_removed
            disp('No features removed; breaking the loop.');
            break;
        end

        % Check if the minimum CMI is greater than 0
        if min(CMI_FS(:, 1)) > 0
            disp('Breaking the loop as no more insignificant features to remove.');
            break;
        end
    end

    disp('Finished, final Sf:');
    disp(Sf);

    MI_Fy=sortrows(MI_Fy,'descend');
    MI_FF=sortrows(MI_FF,'descend');
    CMI_FFy=sortrows(CMI_FFy,'descend');
end





end

