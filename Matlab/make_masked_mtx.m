function [mtx_vals,mask_vals,gen_sim_scores,imp_sim_scores] = make_masked_mtx(mtx_name,mask_name,n_gal,n_probe)
% Format: Output = make_masked_mtx('name of matrix','name of
% mask','n_gal','n_probe)
% ====
% Masks/labels matrices. cd to dir with the matrix and mask
% Number of gallery and probe items are obtained by opening .mask file (in
% benchmarks dir) and looking at the numerical values in line 4. The first
% value is the number of probe times. The second value is the number of
% gallery items.
% ====

%% Read matrix and mask
mtx_vals = read_mtx_mat(mtx_name, n_gal, n_probe);
mask_vals = read_mask_mat(mask_name,n_gal,n_probe);

%% Find similarity scores for genuine and imposter data
%genuine similarity scores
gen_sim_scores = mtx_vals(mask_vals==255);
%imposter similarity scores
imp_sim_scores = mtx_vals(mask_vals==127);

%% Plot data
% All
histogram(gen_sim_scores);
hold on
histogram(imp_sim_scores);

% %Trim weird offshoot number
figure
histogram(gen_sim_scores(gen_sim_scores > min(gen_sim_scores)));
hold on
histogram(imp_sim_scores(imp_sim_scores > min(imp_sim_scores)));

end
