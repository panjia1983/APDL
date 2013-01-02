function multiconlitron_visualize(label_matrix, instance_matrix, vis_matrix, contour_level, downsample_level)
%% svmtoy(label_matrix, instance_matrix, options, contour_level)
%% label_matrix: N by 1, has to be two-class
%% instance_matrix: N by 2
%% options: default '',
%%          see libsvm-mat-8 README, has to be a classification formulation.
%% contour_level: default [0 0], 
%%                change to [-1 0 1] for showing the +/- 1 margin.
%%
%% svmtoy shows the two-class classification boundary of the 2-D data
%% based on libsvm-mat-2.8
%%
%% Hsuan-Tien Lin, htlin at caltech.edu, 2006/04/07


if nargin <= 3
  contour_level = [0 0];
end

if nargin <= 4
  downsample_level = 1;
end

N = size(label_matrix, 1);
if N <= 0
  fprintf(2, 'number of data should be positive\n');
  return;
end

if size(label_matrix, 2) ~= 1
  fprintf(2, 'the label matrix should have only one column\n');
  return;
end

if size(instance_matrix, 1) ~= N
  fprintf(2, ['the label and instance matrices should have the same ' ...
              'number of rows\n']);
  return;
end

if size(instance_matrix, 2) ~= 2
  fprintf(2, 'multiconlitron_visualize only works for 2-D data\n');
  return;
end

X = vis_matrix(:, 1);
Y = vis_matrix(:, 2);
Z = vis_matrix(:, 3);

N = sqrt(size(vis_matrix, 1));

bigX = X(1:N, 1);
bigY = Y(1:N:size(Y, 1), 1);
bigZ = reshape(Z, N, N)';

clf;
hold on;


label_matrix = downsample(label_matrix, downsample_level);
instance_matrix = downsample(instance_matrix, downsample_level);

ispos = (label_matrix == label_matrix(1));
if label_matrix(1) == 1
    pos = find(ispos);
    neg = find(~ispos);
else
    pos = find(~ispos);
    neg = find(ispos);
end

plot(instance_matrix(pos, 1), instance_matrix(pos, 2), 'ko');
plot(instance_matrix(neg, 1), instance_matrix(neg, 2), 'gx');
contour(bigX, bigY, bigZ, contour_level);
colormap winter

