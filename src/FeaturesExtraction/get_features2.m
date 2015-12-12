function get_features2(prefix, N)
addpath('jsonlab');
TOTAL_ITV = 12;
bn = 0;
l = zeros(N, 1);
s = zeros(N, 1);
v_ave = zeros(N, 1);
a_var = zeros(N, 1);
a_votes = zeros(N, TOTAL_ITV); 
for i = 1:N
    filename = strcat(prefix, sprintf('/%d.csv', i));
    fprintf(repmat('\b', 1, bn));
    msg = sprintf('Processing %s ...', filename);
    bn = numel(msg);
    fprintf(msg);
    x = importdata(filename);
    x = x.data;
    n = size(x, 1);
    l(i) = norm(x(n,:)-x(1,:));
    s(i) = 0;
    v_ave(i) = 0;
    a = zeros(n-2, 1);
    a_directed = zeros(n-2, 2);
    for j = 2:n
        s(i) = s(i) + norm(x(j,:)-x(j-1,:));
        v_ave(i) = v_ave(i) + norm(x(j,:)-x(j-1,:));
        if (j>2)
            a(j-2) = (norm(x(j,:)-x(j-1,:)) - norm(x(j-1,:)-x(j-2,:))) / 2;
            
            % normalize according to travel direction
            theta = atan2(x(j-1,2)-x(j-2,2), x(j-1,1)-x(j-2,1));
            a_directed(j-2, :) = (x(j,:)-2*x(j-1,:)+x(j-2,:)) / 2 * ...
                [cos(theta), -sin(theta); sin(theta), cos(theta)];
        end
    end
    a_var(i) = std(a);                                  % var() is likely
    for j = 1:n-2
        mag = norm(a_directed(j, :));
        deg = atan2d(a_directed(j, 2), a_directed(j, 1));
        if (deg<-1e-6)
            deg = deg+360;
        end
        itv = mod(floor(deg/(360/TOTAL_ITV)), TOTAL_ITV)+1;
        a_votes(i, itv) = a_votes(i, itv) + mag;
    end
    % normalize according to total time
    a_votes(i, :) = a_votes(i, :) / (n/60);
end
fprintf('\n');

%figure;
%hist(a_var, 50);
%a_var = log(a_var);
%figure;
%hist(a_var);
[min_var, idx_min] = sort(a_var);
idx_max = wrev(idx_min);
max_var = wrev(min_var);
% scale to 0~180
a_var = 40+100*(1 - (a_var-min_var(3))/(max_var(3)-min_var(3)));
a_var(idx_min(1:2)) = [180 160];
a_var(idx_max(1:2)) = [0 20];
%figure;
%hist(a_var);

min_ave = min(v_ave);
max_ave = max(v_ave);
v_ave = 0.75*(v_ave-min_ave)/(max_ave-min_ave)+0.25;   % scale to 0.25~1

s = s/1000;                     % m -> km
l = l/1000;                     % m -> km

a_var = num2cell(a_var, N);
v_ave = num2cell(v_ave, N);
s = num2cell(s, N);
l = num2cell(l, N);
a_votes = mat2cell(a_votes, ones(1, N), TOTAL_ITV);

features = struct('l', l, 's', s, 'v_ave', v_ave, ...
        'a_var', a_var, 'a_hist', a_votes);
savejson('', ...
    features, ...
    strcat('json_', prefix, '.json'));

end