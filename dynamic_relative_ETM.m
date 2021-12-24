tm=40;h=0.001;t=0:h:tm;
n=1;  
% Design_K_P;
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
max_phi2=0.125;
min_phi2=0.05;
max_delta_r=0.12;
min_delta_r=0.08;
topm2=0.145;
%other parameter
vare=1/5;
Count_2=0;%the flag of triggering times
r2=1;
Gamma=1;
a=10;
b=1;
while n<=(tm/h)+1;
    if n==1
     %% Initialization
        x2(1,:)=[-4,3.5,3.0];
        D_es12(1)=0.0;D_es22(1)=0.0;D_es32(1)=0.0;%The estimate of D
        d1_es2(1)=0;%The estimate of disturbance
        rho_es2(1,:)=[0.0 0.0 0.0];
        vd2(1,:)=[0.0 0.0 0.0];
        u2(1)=0;
        d12(1)=0;intervals2(1)=0;
        z2(1,:)=[0.35 0.35 0.35];
        w2(1)=0;phi2(1)=1.0;
        A=A2;B=B2;C1=C12;C2=C22;P=P2;K=K2;L_o=L_o2;Bw=Bw2;F=F2;
        delta_r(1)=0.08;
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
     %%  ADOB Gamma=1;
        vd2(n,:)=h*(M*rho_es2(n-1,:)'+L_o*(A1*x2(n-1,:)'+B*w2(n-1)))'+vd2(n-1,:);
        rho_es2(n,:)=vd2(n-1,:)-(L_o*x2(n-1,:)')';
        D_es12(n)=h*(-Gamma*sigma(n)*D_es12(n-1)+Gamma*B'*P*x2(n-1,:)'*rho_es2(n-1,1))+D_es12(n-1);
        D_es22(n)=h*(-Gamma*sigma(n)*D_es12(n-1)+Gamma*B'*P*x2(n-1,:)'*rho_es2(n-1,2))+D_es22(n-1);
        D_es32(n)=h*(-Gamma*sigma(n)*D_es12(n-1)+Gamma*B'*P*x2(n-1,:)'*rho_es2(n-1,3))+D_es32(n-1);
        D_es(n-1,:)=[D_es12(n-1)  D_es22(n-1) D_es32(n-1)];
        d1_es2(n)=D_es(n-1,:)*rho_es2(n-1,:)';
     %% Disturance
        z2(n,:)=h*(G*z2(n-1,:)')'+z2(n-1,:);
        d12(n)=C*z2(n-1,:)';
        d2(n,:)=[1/(8+t1),1/(5+2*t1)];

     %% Controller
        uc=K*x2(n-1,:)'-d1_es2(n-1);
        tanhf=[tanh(topm2*x2(n-1,:)*P*B(:,1)/vare)];
        K_f2=x2(n-1,:)*P*B*uc;
        w2(n)=-(1+delta_r(n-1))*(uc*tanh(K_f2/vare)+topm2*tanhf);
        
     %% ET     
        e=w2(n)-u2(n-1);
        phi2(n)=max_phi2-(max_phi2-min_phi2)*tanh((a*e*e)/b*(u2(n-1)*u2(n-1)));
        if phi2(n)>=max_phi2
            phi2(n)=max_phi2;
        end
        if phi2(n)<=min_phi2
            phi2(n)=min_phi2;
        end
        delta_r(n)=max_delta_r-(max_delta_r-min_delta_r)*tanh((1*e*e)/(u2(n-1)*u2(n-1)));
        if delta_r(n)>=max_delta_r
            delta_r(n)=max_delta_r;
        end
        if delta_r(n)<=min_delta_r
            delta_r(n)=min_delta_r;
        end
        if norm(e,inf)>(delta_r(n)*norm(u2(n-1),inf)+phi2(n))
        u2(n)=w2(n);
        Count_2=Count_2+1;
        Lount_1(n)=1;
        intervals2(n)=(n-r2)*h;
        r2=n;
       else
        u2(n)=u2(n-1);
        intervals2(n)=-1;
        Lount_1(n)=-1;
       end 
       %% Plant
        A21=A;
        s=quadv(@(v) expm(A21*(t2-v))*(B*u2(n-1)+B*d12(n-1)+F*d2(n,:)'),t1,t2);
        x2(n,:)=expm(A21*h)*x2(n-1,:)'+s;
    end
    n=n+1;
end
