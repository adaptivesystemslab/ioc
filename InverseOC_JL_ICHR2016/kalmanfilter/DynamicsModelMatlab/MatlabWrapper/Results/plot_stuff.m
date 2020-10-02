orthophoto = imread('E5Snapshot.png');
% unregistered = imread('figure8.png');
% figure, imshow(unregistered)
% cpselect(unregistered, orthophoto)

xlen_px = [574-52];
xlen_m = 50.6;
px_per_m = xlen_px/xlen_m;
figure, imshow(orthophoto)
axis on

start_px = [52,109];
traj_px = [ekf_state(:,1) -ekf_state(:,2)]*px_per_m;
traj_px = traj_px + repmat(start_px,size(traj_px,1),1);

hold on;
plot(traj_px(:,1),traj_px(:,2),'r','LineWidth',5)