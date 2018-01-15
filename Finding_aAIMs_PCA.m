% Finding_aAIMs_PCA

% This script infer aAIMs in ancient individuals divided to 22 groups (see Annot.txt). 
% This script works as follows:
% For each group i, we plot the PCA using the complete SNP set (CSS) dataset. Then using
% the most 100 informative SNPs, we examine the 2 PCA plots visually (the 100-best SNPs and CSS). 
% If they are similar, we keep these 100 and moved on to the next group, otherwise,
% we increase the number of inofrmative SNPs by 100 and examine the plots again. 
% We write down the number of aAIMs found each time.

%% Load data
input_dir = 'D:\My Documents\University\Elhaik Lab\Data\Populations\AncientDNA\aGPS_300\';
load([input_dir 'AllData.mat']);
disp('Data were loaded');

%% load population annotation data
list = load_list([input_dir 'Annot.txt'], 0, 2);

% check that list matches Annot
Annot_Id = {};
Annot_group = {};
for i=1:size(Annot,1)
    spaces_pos = regexp(Annot{i},'\s');
    curr_sample = Annot{i}(spaces_pos+1:end);
    [~, b] = ismember(curr_sample, list{1});
    Annot_Id{i} = cell2str(list{1}(b));
    Annot_group{i} = cell2str(list{2}(b));
end;

%% PCA for each group
% To Run the analysis: 
%   1) Set the number of group_index 1..[number of population groups]
%   starting with 1.
%   2) run this section. You will see 2 PCA plots. If they are similar,
%   exit and start from 1). If they are not similar, click something and
%   you will see another 2 plots.
%   3) all the aAIMs are recordd automatically in total_aAIMs and
%   total_aAIMs_indexes - print them when you are done.
%   4) A suggestion: keep a masked vector of PCA_reduced and record the
%   optimal number of aAIMs you find in each step. If you'll lose count,
%   you won't start again, then just run this section on this vector
%   without the pause command so all the AIMs would be stored in the
%   total_aAIMs variable.

group_index=1;

emptyness = [];

group_names = unique(Annot_group);
PCA_reduced = ones(1,numel(group_names))*100; %starting number of the top SNPs

% Save the rs# of the aAIMs
final_rs = [];
final_rs_index = [];

% Save the indexes of the aAIMs
total_aAIMs = [];
total_aAIMs_indexes = [];

for i=group_index:group_index
    curr_group = group_names(i);
    emptyness_factor = 0.05;
    disp([num2str(i) '. Now analyzing ' char(cell2str(curr_group)) '...']);
    final_group = [];
    
    curr_group_indexes = find(ismember(Annot_group, curr_group));
    final_group = final(curr_group_indexes,:);
    curr_annot = Annot_Id(curr_group_indexes);
    disp(['    Found #' num2str(numel(curr_annot)) ' samples']);
    color_palate = rand(numel(curr_annot), 3);
    
    % Remove missingness
    miss = sum(final_group==5)./size(final_group,1); %find missing positions
    final_clean = final_group(:,miss<=emptyness_factor);
    final_rs = snps(miss<=emptyness_factor);
    final_rs_index = find(miss<=emptyness_factor);
    emptyness(i) = emptyness_factor;
    while size(final_clean,2)<100
        emptyness_factor = emptyness_factor + 0.05;
        final_clean = final_group(:,miss<=emptyness_factor);
        final_rs = snps(miss<=emptyness_factor);
        final_rs_index = find(miss<=emptyness_factor);
        emptyness(i) = emptyness_factor;
    end;
    disp(['    Found #' num2str(size(final_clean,2)) ' positions with emptyness_factor = ' num2str(emptyness_factor)]);

    % Calculate PCA
    disp('Caluclating PCA for the full set...');
    [~, SCORE, va] = pca(double(final_clean));
    explain_var = round((va/sum(va))*100);

    if size(final_clean,1)<=3
        disp('Insufficient samples, skipping...');
        continue;
    end;
    
    %Find informative SNPs
    [scores] = PCAscores(double(final_clean),4); %computes the scores used to identify PCA-correlated SNPs for n statisticall significant PCs
    % informative_SNPs = find(scores>0.025);
    [Y,I] = sort(scores, 2, 'descend'); 

%     %Find the minimal nubmer of informative SNPs
    while 1
%        
        %Do not mask these to get the file printed
        informative_SNPs_indexes = I(1:PCA_reduced(i));
        total_aAIMs = unique([total_aAIMs final_rs(informative_SNPs_indexes)']);
        total_aAIMs_indexes = unique([total_aAIMs_indexes final_rs_index(informative_SNPs_indexes)]);

        disp('Caluclating PCA for the reduced set...');
        [~, SCORE_, va_] = pca(double(final_clean(:,informative_SNPs_indexes)));
        explain_var_ = round((va_/sum(va_))*100);

        figure
        subplot(1,2,1)
            gscatter(SCORE(:,2), SCORE(:,1), curr_annot', color_palate);
            legend off;
    %         plot(SCORE(:,2), SCORE(:,1), 'o');
            hold on;
            text(SCORE(:,2)+1, SCORE(:,1), curr_annot, 'FontSize', 9);
            title([num2str(i) '. Full results for: ' char(cell2str(curr_group))]);
            xlabel(['PC2=' num2str(explain_var(2))]);
            ylabel(['PC1=' num2str(explain_var(1))]);

        subplot(1,2,2)
            gscatter(SCORE_(:,2), SCORE_(:,1), curr_annot', color_palate);
            legend off;
    %         plot(SCORE_(:,2), SCORE_(:,1), 'o', 'MarkerFaceColor', color_palate);
            hold on;
            text(SCORE_(:,2)+1, SCORE_(:,1), curr_annot, 'FontSize', 9);
            title(['Reduced results with ' num2str(PCA_reduced(i)) ' markers']);
            xlabel(['PC2=' num2str(explain_var_(2))]);
            ylabel(['PC1=' num2str(explain_var_(1))]);

        pause
        close;
        PCA_reduced = PCA_reduced + 100;
    end; %while loop

end; %for loop
%%

% The final list of aAIMs is in total_aAIMs


