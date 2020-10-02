function scoreout = pivotgolf(course,pivotstrat)
% PIVOTGOLF  Pivot Pickin' Golf.
%    Your goal is to use LUGUI to compute the LU decompositions of nine
%    matrices with as little roundoff error as possible.  Your score for
%    each hole is norm(R,inf)+norm(Leps,inf)+norm(Ueps,inf) where
%    R = L*U - A(p,q) is the residual and Leps and Ueps are the nonzeros
%    that should be zero in L and U.
%
%    The six golf courses include:
%       magic -- magic squares and Golub matrices, some are rank deficient.
%       testmats -- various test matrices, Pascal, Hilbert, gallery, etc.
%       rand #s -- random integer matrices initialized by rand('state',s).
%
%    The four pivot strategies are:
%       pick -- use the mouse to pick the pivots.
%       diagonal -- pivot on the diagonal, it is possible to divide by zero.
%       partial -- pivot is the largest element in the current column.
%       complete -- pivot is the largest element in the unreduced matrix.
%
%    PIVOTGOLF(course,pivotstrat) bypasses the initial screen.
%       course = 'magic', 'testmats', or a numeric s to set rand('state',s).
%       pivotstrat = 'pick','diagonal','partial', or 'complete'.
%
%    Try to pick pivot elements that divide into the other elements exactly.
%    If you can choose a pivot element that is a power of two, there will
%    be no roundoff error at that step.  But if the pivot is small, the
%    elements in the unreduced matrix might grow larger and subsequent steps
%    might have larger roundoff errors.
%
%    You should know that "mulligan" is golf terminology for "OOPS, I want
%    to take that shot over".  Good luck!
%
%    See also LUGUI.

if nargin == 0

   % Initial screen

   shg
   clf
   axes('pos',[0 0 1 1]);
   axis off
   set(gcf,'double','on','name','Pivot Pickin'' Golf', ...
      'menu','none','numbertitle','off','color','white')
   text('units','norm','pos',[.25,.75],'color',[0 0.65 0], ...
      'fontweight','bold','fontsize',16, ...
      'string',sprintf('Welcome to Pivot Pickin'' Golf'))
   text('units','norm','pos',[.03 .60],'fontweight','bold','string','Course:')
   for k = 1:6
      switch k
         case 1, 
            str = 'magic';
         case 2, 
            str = 'testmats';
         otherwise,
            str = ['rand #' int2str(k-2)];
      end
      b(k) = uicontrol('units','norm','pos',[.14*k .57 .12 .06], ...
         'style','toggle','string',str,'back','w');
   end
   text('units','norm','pos',[.03 .45],'fontweight','bold','string','Strategy:')
   str = {'pick','diagonal','partial','complete'};
   for k = 1:4
      b(k+6) = uicontrol('units','norm','pos',[.14*k .42 .12 .06], ...
         'style','toggle','string',str{k},'back','w');
   end
   stop = uicontrol('style','toggle','string','X','fontweight','bold', ...
      'back','white','units','norm','pos',[.96 .96 .04 .04]);
   uicontrol('units','norm','pos',[.02 .02 .10 .05], ...
      'string','help','back','white','callback','helpwin pivotgolf')
   bvals = [];
   while sum(bvals) < 2
      if get(stop,'val'), break, end
      pause(.05)
      bvals = cell2mat(get(b,'val'));
   end
   if get(stop,'val')
      close(gcf)
      return
   end
   for k = find(bvals)'
      switch k
         case 1, course = 'magic';
         case 2, course = 'testmats';
         case {3,4,5,6}, course = k-2;
         otherwise, pivotstrat = get(b(k),'string');
      end
   end
elseif nargin == 1
   pivotstrat = 'pick';
end

if isnumeric(course)
   rand('state',course)
else
   course = lower(course);
end

% Loop over nine holes (matrices).

score = 0;
h = 1;
A = [];
while h <= 9
   if isempty(A)
      
      % Generate the matrix.

      switch course
         case 'magic'
            n = 2+ceil(2/3*h);
            switch h
               case {3,6,9}
                  A = golub(n);
               otherwise
                  A = magic(n);
            end
         case 'testmats'
            n = h;
            switch h
               case 1
                  n = 7;
                  e = ones(n,1);
                  A = full(spdiags([e (-3:3)' e],[-1 0 1],n,n));
               case 2
                  n = 5;
                  A = vander((-2:2)');
               case 3
                  A = gallery(n);
               case 4
                  A = hadamard(n);
               case 5
                  A = gallery(n);
               case 6
                  A = pascal(n);
               case 7
                  n = 6;
                  A = pascal(n);
                  A(n,n) = A(n,n)-1;
               case 8
                  n = 5;
                  A = 27720*hilb(n);
               case 9
                  n = 6;
                  U = eye(n,n) - triu(ones(n,n),1);
                  A = U'*U;
            end
         otherwise
            n = ceil(2+3*rand);
            A = round(10*(2*rand(n,n)-1));
      end
   end

   % Use LUGUI to compute the LU decomposition

   [L,U,p,q] = lugui(A,pivotstrat);
   pause(2);

   % Score

   R = (L*U - A(p,q));
   Leps = L.*(abs(L)<1000*norm(L,1)*eps);
   Ueps = U.*(abs(U)<1000*norm(U,1)*eps);
   show(abs(R)+abs(Leps)+abs(Ueps));
   if all(isfinite(R(:)))
      s = ceil(4*(norm(R(:),1)+norm(Leps(:),1)+norm(Ueps(:),1))/eps)/4;
   else
      s = Inf;
   end

   % Report the score and decide what to do next.

   set(gcf,'name','Pivot Pickin'' Golf')
   text('units','pixels','pos',[20+50*n,20+(n+2)*30], ...
      'fontweight','bold','fontsize',12, ...
      'color',[0 0.65 0],'string',sprintf('hole #%d',h))
   text('units','pixels','pos',[50*(n-1)+25,55], ...
      'fontweight','bold','fontsize',12,'color',[0 0 0.90], ...
      'string',sprintf('score = %s, total = %s',sph(s),sph(score+s)))
   stop = uicontrol('style','toggle','string','X','fontweight','bold', ...
      'back','w','pos',[100*n+75 30*n+65 25 25]);
   if isequal(pivotstrat,'pick')
      next = uicontrol('units','pixels','pos',[50*(n-1) 10 90 20], ...
         'style','toggle','fontweight','bold', ...
         'background','white','string','next');
         if h==9, set(next,'string','finish'), end
      mulligan = uicontrol('units','pixels','pos',[50*(n+1) 10 90 20], ...
         'style','toggle','fontweight','bold', ...
         'background','white','string','mulligan');
      uics = [stop next mulligan];
      while all(cell2mat(get(uics,'val'))==0)
         drawnow
      end
      if get(mulligan,'val'), continue, end
   else
      pause(3)
   end
   score = score + s;
   if get(stop,'val'), break, end
   h = h + 1;
   A = [];
end

% Final screen

clf
set(gcf,'double','on','name','Pivot Pickin'' Golf', ...
   'pos','default','menu','none','numbertitle','off','color','white')
axes('pos',[0 0 1 1]);
axis off
text('units','norm','pos',[.15,.60], ...
   'fontweight','bold','fontsize',16,'color',[0 0.65 0], ...
   'string',sprintf('Thanks for playing Pivot Pickin'' Golf.'))
if score == 0
   text('units','norm','pos',[.15,.50],...
      'fontweight','bold','fontsize',16,'color',[0 0.65 0], ...
      'string','Perfect score.  Congratulations!')
elseif score == Inf
   text('units','norm','pos',[.15,.50], ...
      'fontweight','bold','fontsize',16,'color',[0 0.65 0], ...
      'color',[0 0.65 0],'string','Your score was infinite.')
   text('units','norm','pos',[.15,.40], ...
      'fontweight','bold','fontsize',16,'color',[0 0.65 0], ...
      'color',[0 0.65 0],'string','You have to avoid dividing by zero.')
else
   text('units','norm','pos',[.15,.50], ...
      'fontweight','bold','fontsize',16,'color',[0 0.65 0], ...
      'string',sprintf('Your score was %s eps.',sph(score)))
   if score > 100
      text('units','norm','pos',[.15,.40], ...
         'fontweight','bold','fontsize',16,'color',[0 0.65 0], ...
         'string','Better luck next time.')
   end
end
uicontrol('style','toggle','units','norm','pos',[.96 .96 .04 .04], ...
   'string','X','fontweight','bold','callback','close(gcf)')
if nargout > 0
   scoreout = score;
   pause(3)
end


%------------------------------------------------------------

function show(A)
% Same code as LUGUI.
clf
axes('pos',[0 0 1 1]);
axis off
[m,n] = size(A);
dx = 100;
dy = 30;
Acolor = [0 0 0];
for j = 1:n
   for i = 1:m
      t(i,j) = text('units','pixels','string',spf(A(i,j)), ...
         'fontname','courier','fontweight','bold','fontsize',14, ...
         'horiz','right','color',Acolor, ...
         'pos',[20+j*dx 20+(m+2-i)*dy]);
   end
end

%------------------------------------------------------------

function s = spf(aij)
% Subfunction to format text strings
if aij == 0
   f = '%10.0f';
elseif (abs(aij) < 1.e-4) | (abs(aij) >= 1.e4) 
   f = '%10.1e';
else
   f = '%10.4f';
end
s = sprintf(f,aij);

%------------------------------------------------------------

function s = sph(x)
% Format text strings that are integer multiples of 1/4.
if x == 0
   f = '%d';
elseif x == round(x);
   f = '%1.0f';
elseif x == round(2*x)/2;
   f = '%2.1f';
else
   f = '%3.2f';
end
s = sprintf(f,x);
