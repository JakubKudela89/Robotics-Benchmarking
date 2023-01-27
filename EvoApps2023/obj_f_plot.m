function [val] = obj_f_plot(x,angles, pnts)
% objective function computation + visualization
clf;
[m,~] = size(pnts);
hold on;
nr_changes = length(x)/6;
nr_pnts = 100; pos = zeros(nr_changes*nr_pnts,3); line = 1;
for j=1:nr_changes
    %[destination,xyz, rpy, err] = IK2_mex(angles, pnts(j,:), x(3*(j-1)+1:3*j));
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
    plot3(pos(line-1,1),pos(line-1,2),pos(line-1,3),'gx','markersize',30,'linewidth',4);
end
plot3(pos(:,1),pos(:,2),pos(:,3),'b','LineWidth',2);
plot3(pos(1,1), pos(1,2), pos(1,3),'kx','markersize',30,'linewidth',4);
len = sum(vecnorm(  diff( pos )  ,2,2));
for j=1:m
    plot3(pnts(j,1),pnts(j,2),pnts(j,3),'rx','markersize',30,'linewidth',4);
    dist(j) = min(vecnorm(  pos-pnts(j,:)  ,2,2));
end
grid on;
val = 100*max(dist) + len;
ax = gca;
ax.View = [59.781022248308403,41.603587270998226];
ax.CameraPosition = [-4.5803   -2.7121    7.1002];

end

