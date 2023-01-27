function [val] = obj_f(x,angles, pnts)
% objective function computation
gamma = 100;

[m,~] = size(pnts);
nr_changes = length(x)/6;
nr_pnts = 100; pos = zeros(nr_changes*nr_pnts,3); line = 1; 
for j=1:nr_changes
    destination = x((j-1)*6+1:j*6);
    for i=1:6
        a(:,i) = linspace(angles(i), (destination(i)), nr_pnts );
    end
    for i=1:nr_pnts 
        A = FK(a(i,:));
        pos(line,:) = A(1:3,end);
        line = line + 1;
    end   
    angles = destination;
end
len = sum(vecnorm(  diff( pos )  ,2,2));
for j=1:m
    dist(j) = min(vecnorm(  pos-pnts(j,:)  ,2,2));
end
val = gamma*max(dist)+len;
end
