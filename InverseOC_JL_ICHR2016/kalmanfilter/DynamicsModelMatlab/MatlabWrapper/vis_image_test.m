%Test Visualizer Image capabilities
clear vis
vis = rlVisualizer('vis',640,480);
vis.add3DModel('C:\aslab\projects\vjoukov\SLAM\3dScene\MRPTScene\MRPTScene\tree1_3ds\Tree1.3ds',[0 0 0 0 0 pi/2]);


tic;
for i=1:100
    vis.update();
    img = vis.takeImage();
    imgrbg = img;imgrbg(:,:,1) = img(:,:,3);imgrbg(:,:,3) = img(:,:,1);
    imshow(imgrbg);
end
toc