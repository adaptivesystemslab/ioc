%%
addpath('..\')
vis = rlVisualizer('vis',640,480);
vis.update();

for i=1:20000
   vis.update();
   [img,imgbool] = vis.getScreenshot();
   a = memory;
   if imgbool
      imshow(img); 
      title(num2str(a.MemUsedMATLAB));
      drawnow;
   end
   
end


%%
Z = peaks;
surf(Z); 
axis tight manual 
set(gca,'nextplot','replacechildren'); 
v = VideoWriter('peaks.avi');
open(v);
for k = 1:10000 
   surf(sin(2*pi*k/20)*Z,Z)
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v);