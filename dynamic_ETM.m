tm=40;h=0.001;t=0:h:tm;
n=1;  
G=[0 0 1.0;
   0 -1.0 0;
   -1.0 0 0];
C=[-4 4 4];
OBSER=obsv(A1,C);
rank(OBSER);
% mode list
List=[2,2,2,1,4,3,3,1,3,4,4,1,3,1,1,2,2,3,3,1,1,4,4,2,2,3,1,4,4,2,2,2,2,1,1,2,3,3,1,1,1,2,2,2,1,4,3,3,1,3,4,4,1,3,1,1,2,2,3,3,1,1,4,4,2,2,3,1,4,4,2,2,2,2,1,1,2,3,3,1,1,1,2,2,2,1,4,3,3,1,3,4,4,1,3,1,1,2,2,3,3,1,1,4,4,2,2,3,1,4,4,2,2,2,2,1,1,2,3,3,1,1,1];
k=2;
% threshold parameter
max_phi1=0.225;
min_phi1=0.075;
topm1=0.271;%topm1>=max_phi
%other parameter
vare=1/5;
Count_1=0;%the flag of triggering times
r1=1;
Gamma=1;
a=10;
b=1;
while n<=(tm/h)+1;
    if n==1
       %% Initialization
        x1(1,:)=[-4,3.5,3.0];
        D_es1(1)=0.0;D_es2(1)=0.0;D_es3(1)=0.0;%The estimate of D
        d1_es(1)=0;%The estimate of disturbance
        rho_es(1,:)=[0.0 0.0 0.0];
        vd(1,:)=[0.0 0.0 0.0];
        u1(1)=0;
        d1(1)=0;
        intervals1(1)=0;
        z(1,:)=[0.35 0.35 0.35];
        w1(1)=0;
        phi1(1)=1.0;
        A=A2;B=B2;C1=C12;C2=C22;P=P2;K=K2;L_o=L_o2;Bw=Bw2;F=F2;       
    else
        t1=h*(n-1);
        t2=h*n;
        sigma(n)=0.01*exp(-5*t1);
     %% Mode
        if mod(n,1000)==0
            flag=List(k);
        switch flag
            case 1
                A=A1;B=B1;C1=C11;C2=C21;P=P1;K=K1;L_o=L_o1;Bw=Bw1;F=F1;
            case 2
                A=A2;B=B2;C1=C12;C2=C22;P=P2;K=K2;L_o=L_o2;Bw=Bw2;F=F2;
            case 3
                A=A3;B=B3;C1=C13;C2=C23;P=P3;K=K3;L_o=L_o3;Bw=Bw3;F=F3;
            case 4
                A=A4;B=B4;C1=C14;C2=C24;P=P4;K=K4;L_o=L_o4;Bw=Bw4;F=F4;
        end
        k=k+1;
        end
     %%  ADOB 
        vd(n,:)=h*(M*rho_es(n-1,:)'+L_o*(A1*x1(n-1,:)'+B*w1(n-1)))'+vd(n-1,:);
        rho_es(n,:)=vd(n-1,:)-(L_o*x1(n-1,:)')';
        D_es1(n)=h*(-Gamma*sigma(n)*D_es1(n-1)+Gamma*B'*P*x1(n-1,:)'*rho_es(n-1,1))+D_es1(n-1);
        D_es2(n)=h*(-Gamma*sigma(n)*D_es2(n-1)+Gamma*B'*P*x1(n-1,:)'*rho_es(n-1,2))+D_es2(n-1);
        D_es3(n)=h*(-Gamma*sigma(n)*D_es3(n-1)+Gamma*B'*P*x1(n-1,:)'*rho_es(n-1,3))+D_es3(n-1);
        D_es(n-1,:)=[D_es1(n-1)  D_es2(n-1) D_es3(n-1)];
        d1_es(n)=D_es(n-1,:)*rho_es(n-1,:)';
     %% Disturance
        z(n,:)=h*(G*z(n-1,:)')'+z(n-1,:);
        d1(n)=C*z(n-1,:)';
        d2(n,:)=[1/(8+t1),1/(5+2*t1)];
     %% Controller
        tanhf=[tanh(topm1*x1(n-1,:)*P*B(:,1)/vare)];
        w1(n)=K*x1(n-1,:)'-d1_es(n);
     %% ET   
        e1=w1(n)-u1(n-1);
        phi1(n)=max_phi1-(max_phi1-min_phi1)*tanh((a*e1*e1)/b*(u1(n-1)*u1(n-1)));
        if phi1(n)>=max_phi1
            phi1(n)=max_phi1;
        end
        if phi1(n)<=min_phi1
            phi1(n)=min_phi1;
        end
       if norm(e1,inf)>=phi1(n)
        u1(n)=w1(n);
        Count_1=Count_1+1;
        Lount_1(n)=1;
        intervals1(n)=(n-r1)*h;
        r1=n;
       else
        u1(n)=u1(n-1);
        intervals1(n)=-1;
        Lount_1(n)=-1;
       end 
       %% Plant
        s=quadv(@(v) expm(A*(t2-v))*(B*u1(n)+B*d1(n-1)+F*d2(n,:)'),t1,t2);
        x1(n,:)=expm(A*h)*x1(n-1,:)'+s;
    end
    n=n+1;
end
