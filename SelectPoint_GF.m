function [x3, Pq, GFD] = SelectPoint_GF( s, n, s1, s2)
%% Codes for point selection based on GF-discrepancy
% Select representative points for Normal and Uniform random variables

% By G.H. Chen 
% Date: 2022-04-24
% References:
% [1] G.H. Chen, D.X. Yang, Direct probability integral method for 
%     stochastic response analysis of static and dynamic structural systems,
%     Computer Methods in Applied Mechanics and Engineering, 357 (2019) 112612.
% [2] G.H. Chen, D.X. Yang, A unified analysis framework of static and 
%     dynamic structural reliabilities based on direct probability integral method,
%     Mechanical Systems and Signal Processing, 158 (2021) 107783.
% [3] G.H. Chen, D.X. Yang, Y.H. Liu, H.C. Guo. System reliability analysis 
%     of static and dynamic structures via direct probability integral method.
%     Computer Methods in Applied Mechanics and Engineering, 388 (2022) 114262
% =========================================================================

% Inputs£º
% s  ¡ª¡ª Number of random variables£»
% n  ¡ª¡ª Number of representative points£»
% s1 ¡ª¡ª Number of Normally distributed variables£»
% s2 ¡ª¡ª Number of Uniformly distributed variables within [0,1], s1 + s2 = s£»
% s=9;  n=2^9;

% Outputs£º
% x3  ¡ª¡ª Representative points£»
% Pq  ¡ª¡ª Assigned probability£»
% GFD ¡ª¡ª GF-discrepancy.

%% Step 1£ºgenerate basic point set 
p=sobolset(s,'Skip',2^3);
x0=net(p,n);

for kk=1:1
if s1~=0
    for i=1:s1
        x1(:,i)=norminv(x0(:,i));    
    end
end
for i=s1+1:s
    x1(:,i)=unifinv(x0(:,i));          
end

%% Step 2: Point set rearrangement
for i=1:n
    for j=1:s
        y1(i,j)=(length(find(x1(:,j)<x1(i,j)))+1/2)/n;
    end
end

if s1~=0
    for j=1:s1
        x2(:,j)=norminv(y1(:,j));
    end
end

for j=s1+1:s
    x2(:,j)=unifinv(y1(:,j));    
end


% Define pbound of Voroni cell
for i=1:s
    band_max(i)=ceil(max(x2(:,i)));
    band_min(i)=floor(min(x2(:,i)));
    band(i,:)=[band_min(i),band_max(i)];
end
x0=[]; x1=[]; y1=[];

%% Step 3:  Calculation of point set using QMCS
Np=1e5;
p2=sobolset(s,'Skip',2^4);
P=net(p2,Np);
inteval=(band(:,2)-band(:,1));
V_Tol=prod(inteval);

xx=bsxfun(@plus,bsxfun(@times,P,inteval'),band(:,1)');
sX=sum(xx.^2,2);
sX2=sum(x2.^2,2);

d=(bsxfun(@plus,bsxfun(@plus,-2*xx*x2',sX),sX2')).^(1/2);
[~,I]=min(d,[],2);

parfor i=1:n
    p_s=[];
    II=find(I==i);
    M(i)=length(find(I==i));
    
    if M(i)==0
        p_s=zeros(1,s);
    else
        for j=1:s1
            p_s(:,j)=normpdf(xx(II,j));
        end
        
        for j=s1+1:s
            p_s(:,j)=unifpdf(xx(II,j));
        end
    end
    Pq1(i)=sum(V_Tol/Np*(prod(p_s,2))); 
end
xx=[];

id1=find(Pq1==0);
x2(id1,:)=[];
Pq1(id1)=[];

Pq=Pq1/sum(Pq1);
sum(Pq)

%% Step 4: Point set rearrangement according to assigned probability
for i=1:size(x2,1)
    for j=1:s
        y2(i,j)=sum(Pq(x2(:,j)<x2(i,j)))+1/2*Pq(i);
    end
end

if s1~=0
    for i=1:s1
        x3(:,i)=norminv(y2(:,i));
    end
end
for i=s1+1:s
    x3(:,i)=unifinv(y2(:,i));
end
x2=x3;

for i=1:size(x2,1)
    for j=1:s
        y2(i,j)=sum(Pq(x2(:,j)<x2(i,j)))+1/2*Pq(i);
    end
end

if s1~=0
    for i=1:s1
        x3(:,i)=norminv(y2(:,i));
    end
end
for i=s1+1:s
    x3(:,i)=unifinv(y2(:,i));
end

end
%% Step 5: Calculation of GF discrepancy
if s1~=0
    for i=1:s1
        F(:,i)=normcdf(x3(:,i));
    end
end
for i=s1+1:s
    F(:,i)=unifcdf(x3(:,i));
end

for j=1:s
    for i=1:size(x2,1)
        Fn(i,j)=sum(Pq(x3(:,j)<x3(i,j)));
    end
    FD(:,j)=abs(Fn(:,j)-F(:,j));
end

GFD=max(max(FD));
disp(['GF-discrepancy=',num2str(GFD)])
x3=x3';

end