%clc;
%clear all;

phi_1=-5;phi_2=-6;phi_3=-7;phi_4=-8;

A1=[-1.46 0 2.428;
     0.1643+0.5*phi_1 -0.4+phi_1 -0.3788;
     0.3107 0 -2.23];
B1=[0.4;
    0.2;
    -0.1];
Bw1=[-0.1;-0.1;0.1];
C11=[0.2 0.1 0];
C21=[0.3 0.3 0];
T1=B1;

A2=[-1.46 0 2.428;
     0.1643+0.5*phi_2 -0.4+phi_2 -0.3788;
     0.3107 0 -2.23];
B2=[0.2;
    0.1;
    -0.1];
Bw2=[0.1;0;0];
C12=[0.2 0.0 0.2];
C22=[0.1 0.1 0];
T2=B2;

A3=[-1.46 0 2.428;
     0.1643+0.5*phi_3 -0.4+phi_3 -0.3788;
     0.3107 0 -2.23];
B3=[0.2;
    0.1;
    -0.2];
Bw3=[0.1;0.1;0];
C13=[0.1 0.1 0];
C23=[0.2 0.1 0];
T3=B3;

A4=[-1.46 0 2.428;
     0.1643+0.5*phi_4 -0.4+phi_4 -0.3788;
     0.3107 0 -2.23];
B4=[0.2;
    0.1;
    -0.2];
Bw4=[0.1;0.1;0.1];
C14=[0.2 0.2 0];
C24=[0.1 0.1 0];
T4=B4;

S=0.1;

F1=[B1*S Bw1];
F2=[B2*S Bw2];
F3=[B3*S Bw3];
F4=[B4*S Bw4];
N=[1;2;3];

lo1=(1-1.5*B1(1,1)-1.4*B1(2,1))/(B1(3,1));
L_o1=[1.5 1.4 lo1;
      3.0 2.8 2*lo1;
      4.5 4.2 3*lo1];
lo2=(1-1.4*B2(1,1)-1.3*B2(2,1))/(B2(3,1));
L_o2=[1.4 1.3 lo2;
      2.8 2.6 2*lo2;
      4.2 3.9 3*lo2];  
lo3=(1-1.6*B3(1,1)-1.5*B3(2,1))/(B3(3,1));
L_o3=[1.6 1.5 lo3;
      3.2 3.0 2*lo3;
      4.8 4.5 3*lo3];  
lo4=(1-1.3*B4(1,1)-1.2*B4(2,1))/(B4(3,1));
L_o4=[1.3 1.2 lo4;
      2.6 2.4 2*lo4;
      3.9 3.6 3*lo4];  
  
  
M=[-10 -5 0;
-10 -15 -15;
-5 0 -15];


% A=[-0.1 0.0;
%    0 -0.1];
% eig(A)
% B=[1;1];
sigma=1.1;
gamma=9;
delta=0.1;
setlmis([]);

P1=lmivar(1,[3,1]);
P2=lmivar(1,[3,1]);
P3=lmivar(1,[3,1]);
P4=lmivar(1,[3,1]);

L1=lmivar(2,[1,3]);
L2=lmivar(2,[1,3]);
L3=lmivar(2,[1,3]);
L4=lmivar(2,[1,3]);

U1=lmivar(1,[1,0]);
U2=lmivar(1,[1,0]);
U3=lmivar(1,[1,0]);
U4=lmivar(1,[1,0]);

J1=lmivar(1,[3,1]);
J2=lmivar(1,[3,1]);
J3=lmivar(1,[3,1]);
J4=lmivar(1,[3,1]);

T112=lmivar(1,[3,1]);
T111=lmivar(1,[3,1]);
T131=lmivar(1,[3,1]);
T133=lmivar(1,[3,1]);

T212=lmivar(1,[3,1]);
T211=lmivar(1,[3,1]);
T231=lmivar(1,[3,1]);
T233=lmivar(1,[3,1]);

Pq1=lmivar(1,[3,1]);
Pq2=lmivar(1,[3,1]);
Pq3=lmivar(1,[3,1]);
Pq4=lmivar(1,[3,1]);
% W2=lmivar(2,[2,2]);
%% MODE1
lmiterm([1 1 1 P1],1,A1,'s');
lmiterm([1 1 1 L1],T1,1,'s');
lmiterm([1 1 1 J1],1,1);
lmiterm([1 1 1 P1],-2.4,1);
lmiterm([1 1 1 P3],2.4,1);
lmiterm([1 1 1 P4],2.4,1);
lmiterm([1 1 1 P2],2.2,1);
lmiterm([1 1 1 P3],-2.2,1);
lmiterm([1 1 1 P4],-2.2,1);
lmiterm([1 1 1 T112],0.005,1);

lmiterm([1 1 3 P1],1,F1);
lmiterm([1 1 4 P1],1,B1);
lmiterm([1 1 5 0],C11');
lmiterm([1 1 6 P2],1,1);
lmiterm([1 1 6 P1],-1,1);
lmiterm([1 1 6 P3],-1,1);
lmiterm([1 1 6 P4],-1,1);
lmiterm([1 1 7 P3],1,1);
lmiterm([1 1 7 P4],1,1);
lmiterm([1 1 10 -L1],-1,1);

lmiterm([1 2 2 Pq1],1,M,'s');
lmiterm([1 2 2 0],delta*sigma*eye(3));
lmiterm([1 2 2 Pq1],-2.4,1);
lmiterm([1 2 2 Pq3],2.4,1);
lmiterm([1 2 2 Pq4],2.4,1);
lmiterm([1 2 2 Pq2],2.2,1);
lmiterm([1 2 2 Pq3],-2.2,1);
lmiterm([1 2 2 Pq4],-2.2,1);
lmiterm([1 2 2 T212],0.005,1);


lmiterm([1 2 3 Pq1],1,L_o1*F1);
lmiterm([1 2 5 0],C21');
lmiterm([1 2 8 Pq2],1,1);
lmiterm([1 2 8 Pq1],-1,1);
lmiterm([1 2 8 Pq3],-1,1);
lmiterm([1 2 8 Pq4],-1,1);
lmiterm([1 2 9 Pq3],1,1);
lmiterm([1 2 9 Pq4],1,1);

lmiterm([1 3 3 0],-gamma*eye(2));
lmiterm([1 4 4 0],-sigma*eye(1));
lmiterm([1 5 5 0],-eye(1));
lmiterm([1 6 6 T112],-1,1);
lmiterm([1 7 7 T111],-1,1);
lmiterm([1 8 8 T212],-1,1);
lmiterm([1 9 9 T211],-1,1);


lmiterm([1 10 10 U1],-1,1,'s');

lmiterm([1 10 11 -P1],B1',1);
lmiterm([1 10 11 -U1],1,-T1');
lmiterm([1 11 11 J1],-1,1);

%% MODE2
lmiterm([2 1 1 P2],1,A2,'s');
lmiterm([2 1 1 L2],T2,1,'s');
lmiterm([2 1 1 J2],1,1);
lmiterm([2 1 1 P2],-1.9,1);
lmiterm([2 1 1 P1],1.9,1);
lmiterm([2 1 1 P4],1.9,1);
lmiterm([2 1 1 P3],1.7,1);
lmiterm([2 1 1 P1],-1.7,1);
lmiterm([2 1 1 P4],-1.7,1);

lmiterm([2 1 3 P2],1,F2);
lmiterm([2 1 4 P2],1,B2);
lmiterm([2 1 5 0],C12');

lmiterm([2 1 6 -L2],-1,1);

lmiterm([2 2 2 Pq2],1,M,'s');
lmiterm([2 2 2 0],delta*sigma*eye(3));
lmiterm([2 2 2 Pq2],-1.9,1);
lmiterm([2 2 2 Pq1],1.9,1);
lmiterm([2 2 2 Pq4],1.9,1);
lmiterm([2 2 2 Pq3],1.7,1);
lmiterm([2 2 2 Pq1],-1.7,1);
lmiterm([2 2 2 Pq4],-1.7,1);


lmiterm([2 2 3 Pq2],1,L_o2*F2);
lmiterm([2 2 5 0],C22');


lmiterm([2 3 3 0],-gamma*eye(2));
lmiterm([2 4 4 0],-sigma*eye(1));
lmiterm([2 5 5 0],-eye(1));


lmiterm([2 6 6 U2],-1,1,'s');

lmiterm([2 6 7 -P2],B2',1);
lmiterm([2 6 7 -U2],1,-T2');
lmiterm([2 7 7 J2],-1,1);

%% MODE3
lmiterm([3 1 1 P3],1,A3,'s');
lmiterm([3 1 1 L3],T3,1,'s');
lmiterm([3 1 1 J3],1,1);
lmiterm([3 1 1 P1],1.2,1);
lmiterm([3 1 1 P2],-1.2,1);
lmiterm([3 1 1 P4],-1.2,1);
lmiterm([3 1 1 P3],-1.8,1);
lmiterm([3 1 1 P2],1.8,1);
lmiterm([3 1 1 P4],1.8,1);
lmiterm([3 1 1 T131],0.005,1);

lmiterm([3 1 3 P3],1,F3);
lmiterm([3 1 4 P3],1,B3);
lmiterm([3 1 5 0],C13');
lmiterm([3 1 6 P1],1,1);
lmiterm([3 1 6 P2],-1,1);
lmiterm([3 1 6 P3],-1,1);
lmiterm([3 1 6 P4],-1,1);
lmiterm([3 1 7 P2],1,1);
lmiterm([3 1 7 P4],1,1);
lmiterm([3 1 10 -L3],-1,1);

lmiterm([3 2 2 Pq3],1,M,'s');
lmiterm([3 2 2 0],delta*sigma*eye(3));
lmiterm([3 2 2 Pq1],1.2,1);
lmiterm([3 2 2 Pq2],-1.2,1);
lmiterm([3 2 2 Pq4],-1.2,1);
lmiterm([3 2 2 Pq3],-1.8,1);
lmiterm([3 2 2 Pq2],1.8,1);
lmiterm([3 2 2 Pq4],1.8,1);
lmiterm([3 2 2 T231],0.005,1);


lmiterm([3 2 3 Pq3],1,L_o3*F3);
lmiterm([3 2 5 0],C23');
lmiterm([3 2 8 Pq1],1,1);
lmiterm([3 2 8 Pq2],-1,1);
lmiterm([3 2 8 Pq3],-1,1);
lmiterm([3 2 8 Pq4],-1,1);
lmiterm([3 2 9 Pq2],1,1);
lmiterm([3 2 9 Pq4],1,1);

lmiterm([3 3 3 0],-gamma*eye(2));
lmiterm([3 4 4 0],-sigma*eye(1));
lmiterm([3 5 5 0],-eye(1));
lmiterm([3 6 6 T131],-1,1);
lmiterm([3 7 7 T133],-1,1);
lmiterm([3 8 8 T231],-1,1);
lmiterm([3 9 9 T233],-1,1);


lmiterm([3 10 10 U3],-1,1,'s');

lmiterm([3 10 11 -P3],B3',1);
lmiterm([3 10 11 -U3],1,-T3');
lmiterm([3 11 11 J3],-1,1);
%% MODE4
lmiterm([4 1 1 P4],1,A4,'s');
lmiterm([4 1 1 L4],T4,1,'s');
lmiterm([4 1 1 J4],1,1);
lmiterm([4 1 1 P1],2.6,1);
lmiterm([4 1 1 P2],-2.6,1);
lmiterm([4 1 1 P3],-2.6,1);
lmiterm([4 1 1 P4],-2.8,1);
lmiterm([4 1 1 P2],2.8,1);
lmiterm([4 1 1 P3],2.8,1);

lmiterm([4 1 3 P4],1,F4);
lmiterm([4 1 4 P4],1,B4);
lmiterm([4 1 5 0],C14');

lmiterm([4 1 6 -L4],-1,1);

lmiterm([4 2 2 Pq4],1,M,'s');
lmiterm([4 2 2 0],delta*sigma*eye(3));
lmiterm([4 2 2 Pq1],2.6,1);
lmiterm([4 2 2 Pq1],-2.6,1);
lmiterm([4 2 2 Pq4],-2.6,1);
lmiterm([4 2 2 Pq4],-2.8,1);
lmiterm([4 2 2 Pq1],2.8,1);
lmiterm([4 2 2 Pq4],2.8,1);

lmiterm([4 2 3 Pq4],1,L_o4*F4);
lmiterm([4 2 5 0],C24');

lmiterm([4 3 3 0],-gamma*eye(2));
lmiterm([4 4 4 0],-sigma*eye(1));
lmiterm([4 5 5 0],-eye(1));

lmiterm([4 6 6 U4],-1,1,'s');

lmiterm([4 6 7 -P4],B4',1);
lmiterm([4 6 7 -U4],1,-T4');
lmiterm([4 7 7 J4],-1,1);
%% LIMITED
lmiterm([-5 1 1 P1],1,1);
lmiterm([-6 1 1 P2],1,1);
lmiterm([-7 1 1 P3],1,1);
lmiterm([-8 1 1 P4],1,1);

lmiterm([-9 1 1 T112],1,1);
lmiterm([-10 1 1 T111],1,1);

lmiterm([-11 1 1 T131],1,1);
lmiterm([-12 1 1 T133],1,1);

lmiterm([-14 1 1 T212],1,1);
lmiterm([-15 1 1 T211],1,1);

lmiterm([-16 1 1 T231],1,1);
lmiterm([-17 1 1 T233],1,1);

lmiterm([-18 1 1 J1],1,1);
lmiterm([-19 1 1 J2],1,1);
lmiterm([-20 1 1 J3],1,1);
lmiterm([-21 1 1 J4],1,1);

lmiterm([-22 1 1 U1],1,1);
lmiterm([-23 1 1 U2],1,1);
lmiterm([-24 1 1 U3],1,1);
lmiterm([-25 1 1 U4],1,1);

lmiterm([-26 1 1 Pq1],1,1);
lmiterm([-27 1 1 Pq2],1,1);
lmiterm([-28 1 1 Pq3],1,1);
lmiterm([-29 1 1 Pq4],1,1);
% lmiterm([12 1 1 P13],1,1);
% lmiterm([-12 1 1 P11],1,1);

% lmiterm([2 1 1 0],1*eye(2));
% lmiterm([4 1 1 P2],1,1);
% lmiterm([-4 1 1 0],2*eye(2));

% lmiterm([31 1 1 U1],1,1);
% lmiterm([31 1 1 0],100);
% lmiterm([32 1 1 U2],1,1);
% lmiterm([32 1 1 0],100);
% lmiterm([33 1 1 U3],1,1);
% lmiterm([33 1 1 0],100);
% lmiterm([34 1 1 U4],1,1);
% lmiterm([34 1 1 0],100);

lmisys=getlmis;
options=[0,0,0,0,0];
[tmin,xfeas]=feasp(lmisys,options);

P1=dec2mat(lmisys,xfeas,P1);
P2=dec2mat(lmisys,xfeas,P2);
P3=dec2mat(lmisys,xfeas,P3);
P4=dec2mat(lmisys,xfeas,P4);
U1=dec2mat(lmisys,xfeas,U1);
U2=dec2mat(lmisys,xfeas,U2);
U3=dec2mat(lmisys,xfeas,U3);
U4=dec2mat(lmisys,xfeas,U4);
L1=dec2mat(lmisys,xfeas,L1);
L2=dec2mat(lmisys,xfeas,L2);
L3=dec2mat(lmisys,xfeas,L3);
L4=dec2mat(lmisys,xfeas,L4);
% Pq=dec2mat(lmisys,xfeas,Pq);

K1=inv(U1)*L1;
K2=inv(U2)*L2;
K3=inv(U3)*L3;
K4=inv(U4)*L4;
if tmin<0
eig(A1+B1*K1) 
eig(A2+B2*K2)
eig(A3+B3*K3)
eig(A4+B4*K4)
end
