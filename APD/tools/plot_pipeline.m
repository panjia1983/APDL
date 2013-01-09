function plot_pipeline(filename)
fid = fopen(filename, 'r');
A = {};
i = 1;
while ~feof(fid)
    line = fgets(fid);
    A{i} = sscanf(line, '%f');
    i = i + 1;
end

f = A{1}(3);
for i = 2 : size(A, 2)
    f1 = A{i}(3);
    for j = 3 : 3: size(A{i}, 1)
        A{i}(j) = A{i}(j) * (f / f1);
    end
end

a = 1;
b = 2; % 6
plot(A{a}(1:3:size(A{a},1)), A{a}(3:3:size(A{a},1)), 'r', A{b}(1:3:size(A{b},1)), A{b}(3:3:size(A{b},1)));

%plot(x, A(1, :), 'r', x, A(2, :), 'g', x, A(3, :), 'b', x, A(4, :), 'c', x, A(5, :), 'm', x, A(6, :), 'y');
%plot(x, A(1, :), 'r', x, A(6, :), 'y');

