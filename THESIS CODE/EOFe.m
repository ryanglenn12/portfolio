function [V,W,DV,EG,PV,E1,E2,E3,E4,E5,E6]=EOFe(Y,ioutpt,itest,t)
% [V,W,DV,EG,PV,E1,E2,E3,E4,E5,E6]=EOFe(Y,ioutpt,itest,t);
% V.951006
% Customized version of EOF for e
% Basic EOF program.  Assumes Y contains N observations of
% p variables (Y is Nxp, and p<N)
% Program standardizes the observations before conducting the analysis.
% Outputs are:
%   V = eigenvectors (columns) defining the principle components (pxp)
%   W = loading series which define the eofs (Nxp)
%   DV = Decomposition of variance (pXp)
%   EG = eigenvalues of correlation matrix (px1)
%   PV = Percentage of variance for each EOF mode (px1)
%   E1,E2, ... , E6 = m (m <= p) EOF matrices whose j th column contain the
%   1, 2, ..., mth EOF component of the jth species (jth column) of Y
%  NOTE: All but E1 have mean = 0.
%   User can ask program to output as many components (<= min(6,p) as req'd.
%  If ioutpt>0, eigenvectors and values are displayed, if ioutput=2, a
%  graphical summary is given.
%  If ioutpt is zero or missing, no output is displayed.
%  If itest=1, some approximate
% statistical tests of eigenvalues are conducted.
% Time values, t, are only required if ioutpt==2.
% Default values:output=0, itest=0
%  Check: If there is no mistake
%   Y=ones(N,1)*mnY+W*V'*diag(sdY)
%  (If you have problems with dimensions, this equation assumes that the
%  that the time series are in the columns of Y, otherwise use the
%   transpose, Y', of Y.)

if nargout>2,
  disp(' '),disp('Good Job! It Worked!!')
  %disp('IF YOU ARE GETTING FUCKED UP RESULTS, IT IS PROBABLY YOUR FAULT.')
  %disp('GET YOUR SHIT STRAIGHT, MY MAN!!!!')

end

if nargin==1,
  ioutpt=0;itest=0;
elseif nargin==2,
  itest=0;
end % if-block


[N,m]=size(Y);
if m>N,
  disp('More variables than samples.  Y replaced by Y''.'),
  Y=Y';
  [N,m]=size(Y);
end  % Assumes data series are in columns of Y

V=zeros(m,m);
mnY=mean(Y);sdY=std(Y);
eN=ones(N,1);
Z=(Y-eN*mnY)*diag(sdY.^(-1));

%if itest==0,clear Y;end

X=Z'*Z/(N-1);  % Covariance matrix of Y
dtx=det(X);
[v,d]=eig(X);
d=diag(d);
[ds,isrt]=sort(d);

% Order the eigenvectors according to pct. of variance.
for j=1:m
	eigval(j)=d(isrt(m-j+1));
	V(:,j)=v(:,isrt(m-j+1));
end %end for-loop

% Normalize signs by making the entry with largest abs. value positive
[jnk,id]=max(abs(V));
for j=1:m
  V(:,j)=V(:,j)*sign(V(id(j),j));
end

W=Z*V;

clear Z;

v=V; % Orthonormal eigenvectors

% Scale wts and evctrs
nrm=sqrt(diag(W'*W));
W=W*diag(nrm.^(-1));
V=V*diag(nrm);  % Scaled eigenvectors
Vn=diag(sdY)*V; % Unscaled eigenvectors

% Calculate desired EOF matrices
if nargout>5,
  nc=nargout-5;  % Number of desired components
  E1= eN*mnY +W(:,1)*V(:,1)'*diag(sdY);
  for j=2:nc,
     eval(['E' num2str(j) '= W(:,j)*V(:,j)''*diag(sdY);']);
  end % for-loop
end % if-block

pcv=100*eigval(:)/sum(d);

%if ioutpt>0,
%pcv=100*eigval(:)/sum(d);
%np=9;
%mp=min(m,np);
%rmp=rem(m,np);
%kp=max(1,fix(m/np)*np-1);
% Set up print formats
%frmt1=['%6.1f'];
%frmt2=frmt1;
%for j=2:mp,
%  frmt1=[frmt1,' %6.1f'];
%end
%for j=2:rmp;
%  frmt2=[frmt2,' %6.1f'];
%end
%frmt1=[frmt1,'\n'];
%frmt2=[frmt2,'\n'];
%disp('Eigenvalues')
%disp(eigval)
%disp(' ')
%disp('Percent of total variance explained by each mode =  ')
%for p1=1:np:kp
%  disp('EOF Mode Number')
%  fprintf(frmt1,(p1:p1+mp-1)');,disp(' ')
%  fprintf(frmt1,pcv(p1:p1+mp-1))
%end
%if p1+mp-1<m,
%  disp(' '),disp('EOF Mode Number')
%  fprintf(frmt2,(np*fix(m/np)+1:m)');,disp(' ')
%  fprintf(frmt2,pcv(np*fix(m/np)+1:m))
%end

%disp(' ')
%disp(' Hit Return for normalized eigenmodes')
%pause
%disp('Corresponding normalized eigenmodes (the columns) = ')
%frmt1=['%7.2f'];
%frmt2=frmt1;
%for j=2:mp,
%  frmt1=[frmt1,' %7.2f'];
%end
%for j=2:rmp;
%  frmt2=[frmt2,' %7.2f'];
%end
%frmt1=[frmt1,'\n'];
%frmt2=[frmt2,'\n'];
%for p1=1:mp:kp
%  disp('EOF Mode Number')
%  fprintf(frmt1,(p1:p1+mp-1)');,disp(' ')
%  fprintf(frmt1,V(:,p1:p1+mp-1)')
%end
%if p1+mp-1<m,
%  disp(' '),disp('EOF Mode Number')
%  fprintf(frmt2,(np*fix(m/np)+1:m)');,disp(' ')
% fprintf(frmt2,V(:,np*fix(m/np)+1:m)')
%end

%disp(' ')
%disp(' Hit Return for variance decomposition')
%pause
%disp('Variance Decomposition')
%disp('Row entries represent percent of row variance assoc. with column mode')
%Var=100*(v.^2)*diag(eigval);
%for p1=1:mp:kp
%  disp('EOF Mode Number')
%  fprintf(frmt1,(p1:p1+mp-1)');,disp(' ')
%  fprintf(frmt1,Var(:,p1:p1+mp-1)')
%end
%if p1+np-1<m,
%  disp(' '),disp('EOF Mode Number')
%  fprintf(frmt2,(np*fix(m/np)+1:m)');,disp(' ')
%  fprintf(frmt2,Var(:,np*fix(m/np)+1:m)')
%end

%if ioutpt==2,  % Graphical summary
%  figure(1), clg
%  for j=1:min(m,4),
%   subplot(2,2,j)
%   bar(V(:,j));
%    ax=axis;
%    axis([.5 m+.5,min(0,1.2*min(V(:,j))) 1.2*max(V(:,j))])
%   title(['EOF ',int2str(j),' :  ',num2str(pcv(j)),'% Total Var.'])
%    if j==3 | j==4,
%      xlabel('Variate Number');
%    end
%  end
%  figure(2),clg
%  if nargin<4
%    t=(1:N)';
%    xlbl='Indices';
%  else
%    xlbl='Time';
%  end
%  for j=1:min(m,4),
%   subplot(4,1,j),
%    plot(t,W(:,j))
%   ax=axis;
%    axis([min(t) max(t) ax(3:4)]);
%   title(['EOF ',int2str(j),' Time Series'])
%  end
% xlabel(xlbl)
%end

%end %output-if-block

if itest==1,
% Test for equality of eigenvalues (See, Kendall, "Course in MV Anal."pg 97)
disp(['Press "return" for Bartlett''s test for equality of eigenvalues']),pause
%clc
df=.5*m*(m-1);
Nf=N-(2*m+11)/6;
chisq=-Nf*log(dtx);
pdlmb=1;smlmb=0;
for k=1:m-2,
   df=[df;.5*(m-k)*(m-k-1)];
   mlt=Nf-.667*k;
   pdlmb=pdlmb*eigval(k);smlmb=smlmb+eigval(k);
   den=pdlmb*((m-smlmb)/(m-k))^(m-k);
   chisq=[chisq;-mlt*log(dtx/den)];
end  %for-loop
disp([' Test of equality of: lmb(1)=...lmb(m);lmb(2)=...lmb(m),...'])
disp('   deg.f.   Chisq. '),
disp([df,chisq]),
disp('Press "return" for confidence intervals for the covariance matrix eigenvalues.'),pause
%clc
[V1,D]=eig(Y'*Y/(N-1));
D=-sort(diag(-D));
lwr=D/(1+2*sqrt(2/(N-1)));upp=D/(1-2*sqrt(2/(N-1)));
disp('Approximate 95% individual confidence intervals.'),
format long
disp('Lwr Bnd.      eig.val       Uppr. Bnd')
disp([lwr,D,upp]),
zval=input('Input z-value for joint confidence intervals. ("0" exits) ');
if zval>0,
    lwr=D/(1+zval*sqrt(2/(N-1)));upp=D/(1-zval*sqrt(2/(N-1)));
    disp('Lwr Bnd.      eig.val       Uppr. Bnd')
    disp([lwr,D,upp]),
end %end-if
format short
end %end-if
DV=100*(v.^2)*diag(eigval);
EG=eigval;
PV=pcv;

PV=PV';
% save eoftable.txt PV V DV /ascii /tabs
save('eoftable.txt','PV','V', 'DV', '-ascii')

