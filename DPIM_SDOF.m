%% This code is for stochastic response analysis of SDOF system using DPIM
% The details of this example can be refered in Ref. [1] 
% This code shows the results by direct probability integral method (DPIM)and analytical method  
%
% By G.H. Chen
% Email:chengh_414@dlut.edu.cn
% Date: 2022-07-20
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
% 

clear;clc
close all

%% Initial Parameters
w1=5*pi/4; w2=7*pi/4;

% Number of presentative points
Nsel=500;  

% Number of random parameters
Num_var=1;
Norm_var=0;   % Number of uniformly distributed random parameters
Unif_var=1;   % Number of normly distributed random parameters

% Number of representative points
Num_point=Nsel+1;

% Generate representative points and assigned probability
[ z, Pq, GFD] = SelectPoint_GF( Num_var, Num_point, Norm_var, Unif_var);

if Norm_var==0
    w=w1+pi/2*z(1,:);         % Uniform distribution
elseif Norm_var==1
    cov = 0.2;
    w=3*pi*(1+cov*z(1,:));   %Norm distribution 
end

% Generate computational mesh within [x, t]
Nx=501; xend=0.15; dx=2*xend/(Nx-1); x=-xend:dx:xend;
Nt=2001; tend=2; dt=tend/(Nt-1); t=[dt:dt:tend]';

% Initial condition of SDOF free vibration
x0=0.1;

%% Analytical solution
ZA=zeros(length(t),length(x));
tic
for ix=1:length(x)
    xi=x(ix);
    X1=2*pi*(0:100)+2*pi-acos(xi/x0);
    X2=2*pi*(0:100)+acos(xi/x0);
    for it=1:length(t)
        ti=t(it);
        if Norm_var==0
            p1=1/ti*unifpdf(X1/ti,w1,w2);
            p2=1/ti*unifpdf(X2/ti,w1,w2);  
        elseif Norm_var==1
            p1=1/ti*normpdf(X1/ti,3*pi,cov*3*pi);
            p2=1/ti*normpdf(X2/ti,3*pi,cov*3*pi); 
        end
        p(it,:) = p1+p2;
    end
    pp = sum(p,2);
    p3 = (1/(sqrt(x0^2-xi^2))).*(abs(xi)<=abs(x0))+0.*(abs(xi)>abs(x0));
    ZA(:,ix) = p3*pp;
end
T1=toc;
disp(['CPU Time taken by Analytical method',num2str(T1)])


%% Direct probability integral method (DPIM)
tic
for i=1:Nsel+1
    g(:,i) = x0*cos(w(i)*t);
end

% Smoothing parameter
sigma=0.02*min(std(g,[],2),iqr(g,2)/1.34)*(Nsel)^(-1/5);

ZDPIM = zeros(length(t),length(x));
for i=1:Nsel+1
    pp1 = Pq(i)*normpdf(x,g(:,i),sigma);
    ZDPIM = ZDPIM+pp1;
end
T2=toc;
disp(['CPU time taken by DPIM',num2str(T2)])

t1=1.1;
It11=abs(t-t1);
t11=find(It11==min(min(It11)));

figure(11)
hold on
plot(x,ZA(t11,:),'r-','LineWidth',1.5)
plot(x,ZDPIM(t11,:),'g--','LineWidth',1.5)
legend('Analytical method', 'DPIM')
xlabel('x')
ylabel('PDF')

