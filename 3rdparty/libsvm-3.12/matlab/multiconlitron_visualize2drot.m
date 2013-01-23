function multiconlitron_visualize2drot(label_matrix, instance_matrix, vis_matrix, downsample_level)
%% svmtoy(label_matrix, instance_matrix, options, contour_level)
%% label_matrix: N by 1, has to be two-class
%% instance_matrix: N by 3
%% options: default '',
%%          see libsvm-mat-8 README, has to be a classification formulation.
%% contour_level: default [0 0], 
%%                change to [-1 0 1] for showing the +/- 1 margin.
%%
%% svmtoy shows the two-class classification boundary of the 2-D data
%% based on libsvm-mat-2.8
%%
%% Hsuan-Tien Lin, htlin at caltech.edu, 2006/04/07


if nargin <= 4
  downsample_level = 1;
end
% 
% N = size(label_matrix, 1);
% if N <= 0
%   fprintf(2, 'number of data should be positive\n');
%   return;
% end
% 
% if size(label_matrix, 2) ~= 1
%   fprintf(2, 'the label matrix should have only one column\n');
%   return;
% end
% 
% if size(instance_matrix, 1) ~= N
%   fprintf(2, ['the label and instance matrices should have the same ' ...
%               'number of rows\n']);
%   return;
% end
% 
% if size(instance_matrix, 2) ~= 3
%   fprintf(2, 'multiconlitron_visualize2drot only works for 3-D data\n');
%   return;
% end

X = vis_matrix(:, 1);
Y = vis_matrix(:, 2);
Z = vis_matrix(:, 3);
V = vis_matrix(:, 4);

N = round((size(X, 1))^(1/3));

bigY = X(1:N, 1);
bigX = Y(1:N:N*N, 1);
bigZ = Z(1:N*N:N*N*N, 1);
bigV = reshape(V, N, N, N);

clf;
hold on;

p = patch(isosurface(bigX, bigY, bigZ, bigV, 0));
isonormals(bigX, bigY, bigZ, bigV, p);
alpha(0.5);

% save obj model 
[f,v] = isosurface(bigX, bigY, bigZ, bigV, 0);
vertface2obj(v,f,'model_multiconlitron_2d_rot.obj');

set(p,'FaceColor','red','EdgeColor','none');

view(3); daspect([1 1 1]); axis tight
camlight;  camlight(-80,-10); lighting phong; 


% label_matrix = downsample(label_matrix, downsample_level);
% instance_matrix = downsample(instance_matrix, downsample_level);
% 
% ispos = (label_matrix == label_matrix(1));
% if label_matrix(1) == 1
%     pos = find(ispos);
%     neg = find(~ispos);
% else
%     pos = find(~ispos);
%     neg = find(ispos);
% end
% 
% scatter3(instance_matrix(pos, 2), instance_matrix(pos, 1), instance_matrix(pos, 3), 'ko');
% scatter3(instance_matrix(neg, 2), instance_matrix(neg, 1), instance_matrix(neg, 3), 'gx');

