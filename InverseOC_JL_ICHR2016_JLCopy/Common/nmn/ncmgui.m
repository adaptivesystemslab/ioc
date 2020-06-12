function ncmgui
%NCMGUI  Master GUI for Numerical Computing with MATLAB.
%  NCMGUI provides an interface to nineteen of the graphical demostrations
%  from "Numerical Computing with MATLAB".  Click on any graphic to access
%  the underlying function.
%     fern       floatgui   pivotgolf  eigshow    interpgui
%     fzerogui   censusgui  quadgui    lorenzgui  swinger
%     touchtone  fftgui     randgui    blackjack  walker
%     eigsvdgui  waves      pdegui     pennymelt  (help)
%
% See also BLACKJACK, CENSUSGUI, EIGSVDGUI, FERN, FFTGUI, FLOATGUI, FZEROGUI,
%  INTERPGUI, LORENZGUI, PIVOTGOLF, PDEGUI, PENNYMELT, QUADGUI, RANDGUI,
%  SWINGER, TOUCHTONE, WALKER, WAVES

% Initial screen

shg
clf reset
axes('pos',[0 0 1 1])
ncmlogo
pause(pi)

% Primary screen

clf
set(gcf,'color','white')
colormap([white(1); jet(13); gray(8); hot(16); copper(16)])
for q = 1:5
   for p = 1:5
      k = 5*(q-1)+p;
      if k <= 19 
         axk = [(p-1)/5 (4-q)/4 1/5 1/4];
      elseif k == 20
         axk = [0 .992 .008 .008];
      elseif k <= 23
         axk = [4/5 (23-k)/12 1/5 1/12];
      else
         break
      end
      ax = axes('pos',axk);
      switch k
         case  1, xfern; f = @fern; s = 'fern';
         case  2, xfloat; f = @floatgui; s = 'floatgui';
         case  3, xlu; f = @pivotgolf; s = 'pivotgolf';
         case  4, xeigshow; f = @eigshow; s = 'eigshow';
         case  5, xinterp; f = @interpgui; s = 'interpgui';
         case  6, xfzero; f = @fzerogui; s = 'fzerogui';
         case  7, xcensus; f = @censusgui; s = 'censusgui';
         case  8, xquad; f = @quadgui; s = 'quadgui';
         case  9, xlorenz; f = @lorenzgui; s = 'lorenzgui';
         case 10, xswinger; f = @swinger; s = 'swinger'; 
         case 11, xtouchtone; f = @touchtone; s = 'touchtone';
         case 12, xfft; f = @fftgui; s = 'fftgui';
         case 13, xrand; f = @randgui; s = 'randgui';
         case 14, xbj; f = @blackjack; s = 'blackjack';
         case 15, xwalker; f = @walker; s = 'walker';
         case 16, xeigsvdgui; f = @eigsvdgui; s = 'eigsvdgui';
         case 17, xwaves; f = @waves; s = 'waves'; 
         case 18, xpde; f = @pdegui; s = 'pdegui';
         case 19, xpenny; f = @pennymelt; s = 'pennymelt'; 
         case 20, f = @xfernf; box on
         case 21, f = 'helpwin ncmgui'; text(.1,.4,f), box on
         case 22, f = 'helpwin ncm'; text(.1,.4,f), box on
         case 23, f = 'close'; text(.1,.4,f), box on
      end
      set(ax,'userdata',f,'xtick',[],'ytick',[], ...
         'buttondownfcn','set(gcf,''userdata'',get(gca,''userdata''))');
      for h = get(ax,'children')', set(h,'buttondownfcn', ...
         'set(gcf,''userdata'',get(gca,''userdata''))');
      end
      if k <= 19
         uicontrol('style','text','back','w','units','normal', ...
            'pos',[(p-1)/5+1/64 (4-q)/4 1/5-1/32 1/32],'string',s);
      end
      drawnow
   end
end

set(gcf,'windowbuttonupfcn',['f = get(gca,''userdata'');' ...
   'if isa(f,''function_handle''), figure, feval(f),' ...
   'else eval(f), end'])


% --------------------------------

function ncmlogo
set(gcf,'color','black','numbertitle','off','menubar','none', ...
   'name','NCM GUI')
L = rot90(membranetx(1,50,20,20),2);
b = (1/16:1/8:15/16)';
[c,h] = contourf(L,b);
set(h,'linewidth',2,'edgecolor','w');
axis equal
axis off
set(gcf,'color',[0 0 1/4])
bluegreen = [0*b b 1/4+3/4*b];
colormap(bluegreen)
ax = axis;
axis(ax-[20 0 0 0]);
text(ax(1)-10,ax(3)+30,{'Numerical Computing','with MATLAB'}, ...
   'color','white','fontsize',2*get(0,'defaulttextfontsize'))

% --------------------------------

function xfern
spy(finitefern(15000))
set(get(gca,'child'),'color',[0 2/3 0])
axis off
set(gcf,'name','NCM GUI')

% --------------------------------

function xfernf
spy
i = get(get(gca,'child'),'ydata');
j = get(get(gca,'child'),'xdata');
S = sparse(500+i,290+j,1,768,1024);
clf
shg
set(gcf,'doublebuff','on')
axes('pos',[0 0 1 1])
F = imread('fern.png');
F = double((F(:,:,3)==255));
gci = image(2*F+1);
colormap([0 0 2/3; 1/2 1/2 2/3; 1 1 1])
axis off
F = xor(F,S);
set(gcf,'userd',F);
set(gca,'userd',S);
set(gci,'userd',1);
uicontrol('pos',[20 20 15 15],'style','push','callback', ...
   ['F=get(gcf,''userd''); S=get(gca,''userd''); gci=get(gca,''child'');' ...
    'b=1-get(gci,''userd''); set(gci,''userd'',b,''cdata'',2*xor(b*S,F)+1)']);
uicontrol('pos',[40 20 15 15],'style','push','callback', ...
   ['F=get(gcf,''userd''); S=get(gca,''userd''); gci=get(gca,''child'');' ...
    'b=1-get(gci,''userd''); set(gci,''userd'',b,''cdata'',max(2*b*S,F)+1)']);


% --------------------------------

function xfloat
f = (0:7)/8;
F = [];
for e = -1:2
   F = [F (1+f)*2^e];
end
for x = F
   text(x,0,'|')
end
axis([-1 9 -1 1])
text(1,-.2,'1','fontsize',8)
text(2,-.2,'2','fontsize',8)
text(4,-.2,'4','fontsize',8)

% --------------------------------

function xlu
n = ceil(3+3*rand);
A = floor(10*rand(n));
for j = 1:n
   for i = 1:n
      t(i,j) = text('units','norm','string',sprintf('%7.0f',A(i,j)), ...
         'fontname','courier','fontweight','bold','fontsize',2, ...
         'horiz','right','color','black','pos',[(j+1/2)/(n+1) (n+2-i)/(n+2)]);
   end
end
k = min(find(abs(A) == max(max(abs(A)))));
set(t(k),'color','magenta');

% --------------------------------

function xeigshow
A = [1 3; 4 2]/4;
p = 48;
x = exp(2*pi*i*(0:p)/p);
Ax = (A(1,1)+A(2,1)*i)*real(x)+(A(1,2)+i*A(2,2))*imag(x);
k = ceil(p*rand);
h = plot(real(x),imag(x),'.',real(Ax),imag(Ax),'.', ...
    [0 real(x(k))],[0 imag(x(k))], [0 real(Ax(k))],[0 imag(Ax(k))]);
set(h([1 3]),'color',[0 .6 0])
set(h([2 4]),'color',[0 0 .8])
set(h(3:4),'linewidth',2)
s = 1.3*max(1,norm(A));
set(gca,'xlim',[-s s],'ylim',[-s s])

% --------------------------------

function xinterp
x = 1:6;
y = 1/4 + rand(1,6)/2;
u = .75:.05:6.25;
plot(x,y,'o',u,[spline(x,y,u); pchip(x,y,u)])
set(gca,'xlim',[-1 8],'ylim',[0 1],'xtick',[],'ytick',[])

% --------------------------------

function xfzero
F = inline('x.^3-2*x-5');
a = 0;
b = 3;
fa = F(a);
fb = F(b);
c = a;
fc = fa;
xm = (a + b)/2.0;
s = fb/fc;
p = s*(a-b);
q = 1.0 - s;
if p > 0, q = -q; else p = -p; end;
xi = b + p/q;
t = (-11/10:1/100:11/10);
x = a + (b - a)*(1+t)/2;
y = F(x);
xl = min(x);
xr = max(x);
ym = 1.1*max(abs(y));
z = [xl xr];
w = fb+[xl-b xr-b]*(fc-fb)/(c-b);
hp = plot([xl xr xr xl xl],[-ym -ym ym ym -ym],'k-', ...
     x,y,'k:', [xl xr],[0 0],'k-', [a b c],[0 0 0],'kx', ...
     [a b c],[fa fb fc],'ko', [(a+b)/2 (a+b)/2],[-ym/8 ym/8],'r-', ...
     xm,0,'rx', z,w,'-', xi,0,'gx');
axis([a-1 b+1 -ym-10 ym+10])
set(hp([7 9]),'markersize',6,'linewidth',2);
set(hp(8),'color',[0 .75 0])

% --------------------------------

function xquad
x = [0:2:8 9:24 26:2:64]/64;
y = humps(x);
plot(x,y,'.')
hold on
for p = 22:23
   fill(x([p-1 p-1 p p]),[0 y([p-1 p]) 0],[0 .75 0])
end
hold off
axis([-.2 1.2 -20 120])

% --------------------------------

function xlorenz
sigma = 10;
rho = 28;
beta = 8/3;
eta = sqrt(72);
A = [ -beta    0     eta
         0  -sigma   sigma 
      -eta   rho    -1  ];
v0 = [rho-1; eta; eta];
y0 = v0 + [3; 2; -4] + randn(3,1);
tspan = [0 15];
rtol = 1.e-3;
[t,y] = ode23tx(@lorenzeqn, tspan, y0, rtol, A);
plot(y(:,2),y(:,3),'-','color',[0 .75 0])
axis([-30 30 -30 30]);

% ------------------------------

function ydot = lorenzeqn(t,y,A)
A(1,3) = y(2);
A(3,1) = -y(2);
ydot = A*y;

% --------------------------------

function xcensus
p = [ 75.995  91.972 105.711 123.203 131.669 150.697 ...
     179.323 203.212 226.505 249.633 281.422]';
t = (1900:10:2000)';
x = (1890:1:2019)';
w = 2010;
d = 3;
c = polyfit((t-1950)/50,p,d);
y = polyval(c,(x-1950)/50);
z = polyval(c,(w-1950)/50);
h = plot(t,p,'.',x,y,'-',w,z,'.');
axis([1870 2039 -80 480])
set(h(1),'markersize',12)
set(h(2),'linewidth',1)
set(h(3),'markersize',18,'color',[0 .75 0])
text(w-15,z+10,'?','fontweight','bold','color',[0 .75 0])

% --------------------------------

function xfft
n = 16;
x = 1:n;
y = rand(1,n)/2-1/4;
z = fft(y);
u = real(z);
plot([0 n+1],[0 0],'k-', [x;x],[0*u;u],'c-', x,u,'b.','markersize',16)
axis([-2 n+3 -2.5 2])

% --------------------------------

function xrand
n = 200;
X = 2*rand(3,n)-1;
r = sum(X.^2);
k = (r <= 1);
plot3(X(1,k),X(2,k),X(3,k),'r.',X(1,~k),X(2,~k),X(3,~k),'b.','markersize',8)
axis(1.5*[-1 1 -1 1 -1 1])
set(gca,'ztick',[])
axis off

% --------------------------------

function xpde
xv = [-1 -3 -3 1 1 3 1 -1 -1];
yv = [-3 -1 1 1 3 1 -1 -1 -3];
h = 1/8;
[x,y] = meshgrid(-3:h:3);
[in,on] = inregion(x,y,xv,yv);
p = find(in-on==1);
q = find(in-on==0);
G = zeros(size(x));
G(p) = 1:length(p);
A = delsq(G);
opts.disp = 0;
[v,d] = eigs(-A,9,0,opts);
U = zeros(size(x));
k = 4;
i = find(abs(v(:,k)) == max(abs(v(:,k))));
U(p) = v(:,k)/v(i,k);
U(q) = NaN;
s = 1/4;
c = [-1:s:-s -1.e-5 1.e-5 s:s:1];
contourf(U,c);
axis([-6 56 -6 52])
set(gca,'clim',[-1.25 6.7])

% --------------------------------

function xpenny
load penny
contour(flipud(P),0:32:255)
set(gca,'clim',[-511 255])
axis([-32 127+32 -48 127+16])

% --------------------------------

function xwaves
% Circular sector
m = 11;
mu = [3.37561065, 4.27534072, 5.13562230, 6.53025594]';
[r,theta] = meshgrid((0:m)/m,(3/4)*(0:2*m)/m*pi);
V{1} = besselj(2/3,mu(1)*r).*sin(2/3*theta);
V{2} = besselj(4/3,mu(2)*r).*sin(4/3*theta);
V{3} = besselj(2,mu(3)*r).*sin(2*theta);
V{4} = besselj(2/3,mu(4)*r).*sin(2/3*theta);
x = r.*cos(theta+pi);
y = r.*sin(theta+pi);
t = -1.234;
U = zeros(size(V{1}));
for k = 1:4
   U = U + 1/sqrt(k)*sin(mu(k)*t)*V{k};
end
surf(x,y,U);
axis([-1 1 -1 1 -1.5 1.5]);
set(gca,'clim',[-3 3])
view(225,30);
axis off

% --------------------------------

function xtouchtone
load touchtone
image(double(D)/32+15)
set(gca,'pos',get(gca,'pos')+[1  1 -2 -2]/20)

% --------------------------------

function xswinger
s = .9/sqrt(2);
plot([0 s 0],[0 -s -2*s],'o-')
axis([-2 2 -2.75 1.75])

% --------------------------------

function xeigsvdgui
A = diag(13:-1:3)+diag(11:-1:2,1)+diag(11:-1:2,-1)+1;
image(A);
axis([-2 14 -1 16])

% --------------------------------

function xbj
axis([-2.5 2.5 -2.5 2.5])
box on
x = [1/2 -1/2 -1/2 1/2 1/2];
y = [3/4 3/4 -3/4 -3/4 3/4];
offwhite = [.9 .9 .7];
patch(x-.7,y,offwhite)
text(-1.1,0,'A','color','black')
text(-0.65,0,char(170),'fontname','symbol','color','black')
patch(x+.7,y,offwhite)
text(0.5,0,'J','color','red')
text(0.75,0,char(169),'fontname','symbol','color','red')

% --------------------------------

function xwalker
load walkers
L = {[1 5 12],2:8,9:15};
X = reshape(F*[1;0;1;0;1],15,3);
ms = ceil(get(0,'defaultlinemarkersize')/2);
for k = 1:3
   line(X(L{k},1),X(L{k},2),X(L{k},3),'linestyle','-', ...
      'marker','o','markersize',ms);
end
axis([-600 600 -600 600 -200 1600])
set(gca,'xtick',[],'ytick',[],'ztick',[])
view(160,10)
axis off
