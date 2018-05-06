%% System Identificaiton Midterm P3
% By: Sam Nazari
% Last Modified: Today
clear all
clc

%% init
% initialization
M  = 350
N  = 300;
Ts = 1/100;
A = tf([1 0 0],[1 -0.5 0.5],Ts);
BA = tf([1 0 0],[1 -0.5 0.5],Ts);
AC = tf([1 0 0 0],[1 0.35 0.075 0.425],Ts);
TStop = N*Ts;
sim('P3SysID','StopTime','TStop')
y = P3Data.OutputData;
u = P3Data.InputData;
e = 0;
mse = zeros(M,1);

%%
yn=awgn(y(1:N),0); % add white gaussian noise
Y1=[0;yn(1:N-1)];
Y2=[0;0;yn(1:N-2)];
X1=[0;u(1:N-1)];
% X2=[0;0;u(1:N-2)];
% U=[Y1 Y2 X1 X2];
U = [Y1 Y2 X1];
theta=(U'*U)\(U'*yn);
a1=-theta(1);
a2=-theta(2);
b1=theta(3);
%b2=theta(4);
n=0;
%% eHat
for n=1:M
    ehat=zeros(N,1);
    for k=1:1
        ehat(k)=yn(k);
    end
    for k=2:2
        ehat(k)=yn(k)+a1*yn(k-1)-b1*u(k-1);
    end
    for k=3:N
        ehat(k)=yn(k)+a1*yn(k-1)+a2*yn(k-2)-b1*u(k-1);
    end
    Ue1=[0;-ehat(1:N-1)];
    Ue2=[0;0;-ehat(1:N-2)];
    Ue3=[0;0;0;-ehat(1:N-3)];
    Ue=[Ue1 Ue2 Ue3];
    d=(Ue'*Ue)\(Ue'*ehat);
    
%% x_t
    
    unew=zeros(N,1);
    for k=1:1
        unew(k)=u(k);
    end
    for k=2:2
        unew(k)=u(k)+d(1)*u(k-1);
    end
    for k=3:N
        unew(k)=u(k)+d(1)*u(k-1)+d(2)*u(k-2);
    end
    
    %% y_t
    ynew=zeros(N,1);
    for k=1:1
        ynew(k)=yn(k);
    end
    for k=2:2
        ynew(k)=yn(k)+d(1)*yn(k-1);
    end
    for k=3:N
        ynew(k)=yn(k)+d(1)*yn(k-1)+d(2)*yn(k-2);
    end

    Ynew1=[0;-ynew(1:N-1)];
    Ynew2=[0;0;-ynew(1:N-2)];
    Xnew1=[0;unew(1:N-1)];
    Unew=[Ynew1 Ynew2 Xnew1];
    theta=(Unew'*Unew)\(Unew'*ynew);
    a1=theta(1);
    a2=theta(2);
    b1=theta(3);
    sys = idpoly([1 theta(1:2)'],[theta(3)],[1],[1 d'],[1],0.64,Ts);
    yhat = sim(sys,P3Data.InputData);
    e = P3Data.OutputData-yhat;
    mse(n) = (e'*e)/N;
end % while

theta
d

%% Using SystemID Toolbox
% ARARX Model Estimation

na = 2;
nb = 1;
nc = 0;
nd = 3;
nf = 0;
nk = 0;

opt = polyestOptions('SearchMethod','auto','Display','on')
opt.SearchOption.MaxIter = 100
opt.SearchOption.Tolerance= 1e-6
sys = polyest(P3Data,[na nb nc nd nf nk],opt)

%% Set up for P3 b

init_sys = init(sys);
sysPEM = pem(P3Data,init_sys,opt)

%% Plot GLS MSE
plot(1:M,mse),xlim([1 50]),title('MSE vs Iteration'),xlabel('Nth Step'),ylabel('MSE')
