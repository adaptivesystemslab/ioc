function pennymelt(delta)
% PENNYMELT  Heat a penny.
% Initial value of the height is obtained from measurements
% made at the National Institute of Science and Technology
% of the depth of a mold for a U. S. one cent coin.
% What is the limiting value of the height as t -> inf ?
% pennymelt(delta) takes time steps of size delta.
% For what values of delta is the computation stable?

if nargin < 1, delta = .1; end

% Initial conditions

load penny
U = flipud(P);

% Lighted surface plot

shg
clf
surfu = surf(U);
daspect([1,1,128])
colormap(copper)
shading interp
material metal
lighting gouraud
view(2)
axis tight
axis off
set(gca,'zlimmode','auto','climmode','manual');
light('pos',[1,2,2000],'style','inf');
toggle = uicontrol('units','normal','pos',[.02 .02 .10 .05], ...
   'string','start','style','toggle');
while get(toggle,'value') == 0
   drawnow
end

% Finite difference indices

h = 1;
sigma = delta/h^2;
[p,q] = size(U);
n = [2:p p];
s = [1 1:p-1];
e = [2:q q];
w = [1 1:q-1];

% Diffusion

set(toggle,'value',0,'string','stop')
titl = title('0');
t = 0;
while get(toggle,'value') == 0
   U = U + sigma*(U(n,:)+U(s,:)+U(:,e)+U(:,w)-4*U);
   t = t + delta;
   set(surfu,'zdata',U,'cdata',U)
   set(titl,'string',sprintf('%10.3f',t))
   drawnow
end
set(toggle,'string','close','value',0,'callback','close(gcf)')
