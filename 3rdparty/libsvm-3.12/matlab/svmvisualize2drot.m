function svmvisualize2drot(label_matrix, instance_matrix, model_name, scaler_name, downsample_level)
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

if nargin <= 3
    scaler_name = '';
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

if size(instance_matrix, 2) ~= 3
  fprintf(2, 'svmvisualize2drot only works for 3-D data\n');
  return;
end

mdl = loadmodel(model_name, 3);

nclass = mdl.nr_class;
svmtype = mdl.Parameters(1);

if nclass ~= 2 || svmtype >= 2
  fprintf(2, ['cannot plot the decision boundary for these ' ...
              'SVM problems\n']);
  return
end

minX = min(instance_matrix(:, 1));
maxX = max(instance_matrix(:, 1));
minY = min(instance_matrix(:, 2));
maxY = max(instance_matrix(:, 2));
minZ = min(instance_matrix(:, 3));
maxZ = max(instance_matrix(:, 3));

gridX = (maxX - minX) / 100;
gridY = (maxY - minY) / 100;
gridZ = (maxZ - minZ) / 100;

minX = minX - 10 * gridX;
maxX = maxX + 10 * gridX;
minY = minY - 10 * gridY;
maxY = maxY + 10 * gridY;
minZ = minZ - 10 * gridZ;
maxZ = maxZ + 10 * gridZ;

[bigX, bigY, bigZ] = meshgrid(minX:gridX:maxX, minY:gridY:maxY, minZ:gridZ:maxZ);

mdl.Parameters(1) = 3; % the trick to get the decision values
ntest=size(bigX, 1) * size(bigX, 2) * size(bigX, 3);
instance_test=[reshape(bigX, ntest, 1), reshape(bigY, ntest, 1), reshape(bigZ, ntest, 1)];

if ~strcmp(scaler_name, '')
    scaler_matrix = load(scaler_name);
    instance_test(:, 1) = instance_test(:, 1) .* scaler_matrix(1, 1) + scaler_matrix(1, 2);
    instance_test(:, 2) = instance_test(:, 2) .* scaler_matrix(2, 1) + scaler_matrix(2, 2);
    instance_test(:, 3) = instance_test(:, 3) .* scaler_matrix(3, 1) + scaler_matrix(3, 2);
end

label_test = zeros(size(instance_test, 1), 1);

[V, acc] = svmpredict(label_test, instance_test, mdl);

bigV = reshape(V, size(bigX, 1), size(bigX, 2), size(bigX, 3));

clf;
hold on;

p = patch(isosurface(bigX, bigY, bigZ, bigV, 0));
isonormals(bigX, bigY, bigZ, bigV, p);
alpha(0.5);

set(p,'FaceColor','red','EdgeColor','none');

view(3); daspect([1 1 1]); axis tight
camlight;  camlight(-80,-10); lighting phong; 


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

scatter3(instance_matrix(pos, 1), instance_matrix(pos, 2), instance_matrix(pos, 3), 'ko');
scatter3(instance_matrix(neg, 1), instance_matrix(neg, 2), instance_matrix(neg, 3), 'gx');

