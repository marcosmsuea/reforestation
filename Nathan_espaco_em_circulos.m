% EDP do crescimento da Mata
%
clc
clear all
%
%% entrada de parametros do problema
% Esp¨¦cie a: GOIABA DE ANTA, vermelha
% Esp¨¦cie b: LACRE, preta
% Esp¨¦cie c: UCUUBA, Azul
% Esp¨¦cie d: ANDIROBA, verde
alfa_A = 0.01; alfa_B = 0.015; alfa_C = 0.02; alfa_D = 0.025; % difusão
u_A = 0.005; u_B = 0.005; u_C = 0.01; u_D = 0.01; % advecção horizontal
v_A = 0.001; v_B = 0.002; v_C = 0.003; v_D = 0.004; % advecção vertical
mu_A = 0.0; mu_B = 0.0; mu_C = 0.0; mu_D = 0.0; % decaimentos
lambda_A = 0.8; lambda_B = 0.6; lambda_C = 0.8; lambda_D = 0.7; % taxas de crescimento intrínseco
gama_A = 0.01; gama_B = 0.02; gama_C = 0.03; gama_D = 0.04; % competição interespecíficas
sigma_A = 0.02; sigma_B = 0.03; sigma_C = 0.04; sigma_D = 0.05; % competição interespecíficas
delta_A = 0.1; delta_B = 0.2; delta_C = 0.15; delta_D = 0.1; % competição interespecíficas
% capacidade de total de suporte
K=20; 
% parametros do dominio
l=6.0; h=3; tf=0;
%
% parametros da discretizacao
nx=32;  ny=25; nnx=nx+1; nny=ny+1; nn=nnx*nny;
dx=l/nx; dy=h/ny; 
npt=80; dt=tf/npt; dt2=dt/2;

[dx dy dt];

% Número de Péclet
u_Ax=u_A*dx/alfa_A; v_Ay=v_A*dy/alfa_A;
u_Bx=u_B*dx/alfa_B; v_By=v_B*dy/alfa_B;
u_Cx=u_C*dx/alfa_C; v_Cy=v_C*dy/alfa_C;
u_Dx=u_D*dx/alfa_D; v_Dy=v_D*dy/alfa_D;
%
[u_Ax v_Ay u_Bx v_By u_Cx v_Cy u_Dx v_Dy];

%
% cálculos auxiliares
ddx=dx*dx; ddy=dy*dy;
dd2x=2*ddx; dd2y= 2*ddy;

atx_A=alfa_A*dt/dd2x; aty_A=alfa_A*dt/dd2y;
atx_B=alfa_B*dt/dd2x; aty_B=alfa_B*dt/dd2y;
atx_C=alfa_C*dt/dd2x; aty_C=alfa_C*dt/dd2y;
atx_D=alfa_D*dt/dd2x; aty_D=alfa_D*dt/dd2y;
d4x=4*dx; d4y= 4*dy;

u_Adx_A= u_A*dt/d4x; v_Ady_A=v_A*dt/d4y;
u_Bdx_B= u_B*dt/d4x; v_Bdy_B=v_B*dt/d4y;
u_Cdx_C= u_C*dt/d4x; v_Cdy_C=v_C*dt/d4y;
u_Ddx_D= u_D*dt/d4x; v_Ddy_D=v_D*dt/d4y;
dAle=(mu_A-lambda_A)*dt2;
dBle=(mu_B-lambda_B)*dt2;
dCle=(mu_C-lambda_C)*dt2;
dDle=(mu_D-lambda_D)*dt2;
dAld=(-mu_A+lambda_A)*dt2;
dBld=(-mu_B+lambda_B)*dt2;
dCld=(-mu_C+lambda_C)*dt2;
dDld=(-mu_D+lambda_D)*dt2;
Ae=zeros(nn); Ad=zeros(nn);
Be=zeros(nn); Bd=zeros(nn);
Ce=zeros(nn); Cd=zeros(nn);
De=zeros(nn); Dd=zeros(nn);
%
% montagem das matrizes - parte linear apenas
for i=1:nn
Ae(i,i)=1+2*(atx_A+aty_A)+dAle;
Ad(i,i)=1-2*(atx_A+aty_A)+dAld;
Be(i,i)=1+2*(atx_B+aty_B)+dBle;
Bd(i,i)=1-2*(atx_B+aty_B)+dBld;
Ce(i,i)=1+2*(atx_C+aty_C)+dCle;
Cd(i,i)=1-2*(atx_C+aty_C)+dCld;
De(i,i)=1+2*(atx_D+aty_D)+dDle;
Dd(i,i)=1-2*(atx_D+aty_D)+dDld;
end 
for i= 1:nn-1;
Ae(i+1,i)=-aty_A-v_Ady_A;
Ad(i+1,i)=aty_A+v_Ady_A;
Ae(i,i+1)=-aty_A+v_Ady_A;
Ad(i,i+1)=aty_A-v_Ady_A;
Be(i+1,i)=-aty_B-v_Bdy_B;
Bd(i+1,i)=aty_B+v_Bdy_B;
Be(i,i+1)=-aty_B+v_Bdy_B;
Bd(i,i+1)=aty_B-v_Bdy_B;
Ce(i+1,i)=-aty_C-v_Cdy_C;
Cd(i+1,i)=aty_C+v_Cdy_C;
Ce(i,i+1)=-aty_C+v_Cdy_C;
Cd(i,i+1)=aty_C-v_Cdy_C;
De(i+1,i)=-aty_D-v_Ddy_D;
Dd(i+1,i)=aty_D+v_Ddy_D;
De(i,i+1)=-aty_D+v_Ddy_D;
Dd(i,i+1)=aty_D-v_Ddy_D;
end 
for i= 1:nn-nny;
Ae(i+nny,i)=-atx_A-u_Adx_A;
Ad(i+nny,i)=atx_A+u_Adx_A;
Ae(i,i+nny)=-atx_A+u_Adx_A;
Ad(i,i+nny)=atx_A-u_Adx_A;
Be(i+nny,i)=-atx_B-u_Bdx_B;
Bd(i+nny,i)=atx_B+u_Bdx_B;
Be(i,i+nny)=-atx_B+u_Bdx_B;
Bd(i,i+nny)=atx_B-u_Bdx_B;
Ce(i+nny,i)=-atx_C-u_Cdx_C;
Cd(i+nny,i)=atx_C+u_Cdx_C;
Ce(i,i+nny)=-atx_C+u_Cdx_C;
Cd(i,i+nny)=atx_C-u_Cdx_C;
De(i+nny,i)=-atx_D-u_Ddx_D;
Dd(i+nny,i)=atx_D+u_Ddx_D;
De(i,i+nny)=-atx_D+u_Ddx_D;
Dd(i,i+nny)=atx_D-u_Ddx_D;
end 
%adequar a borda horizontal inferior 
for i=1:nnx-1;
ind=i*nny+ny;
Ae(ind,ind-1)=0;
Ad(ind,ind-1)=0;
Ae(ind,ind+1)=-2*aty_A;
Ad(ind,ind+1)=2*aty_A;
Be(ind,ind-1)=0;
Bd(ind,ind-1)=0;
Be(ind,ind+1)=-2*aty_B;
Bd(ind,ind+1)=2*aty_B;
Ce(ind,ind-1)=0;
Cd(ind,ind-1)=0;
Ce(ind,ind+1)=-2*aty_C;
Cd(ind,ind+1)=2*aty_C;
De(ind,ind-1)=0;
Dd(ind,ind-1)=0;
De(ind,ind+1)=-2*aty_D;
Dd(ind,ind+1)=2*aty_D;
end
%
%falta ind=1
Ae(1,2)=-2*aty_A;
Ad(1,2)=2*aty_A;
Be(1,2)=-2*aty_B;
Bd(1,2)=2*aty_B;
Ce(1,2)=-2*aty_C;
Cd(1,2)=2*aty_C;
De(1,2)=-2*aty_D;
Dd(1,2)=2*aty_D;

%
%adequar a borda horizontal superior 
for i=1:nnx-1;
ind=i*nny;
Ae(ind,ind+1)=0;
Ad(ind,ind+1)=0;
Ae(ind,ind-1)=-2*aty_A;
Ad(ind,ind-1)=2*aty_A;
Be(ind,ind+1)=0;
Bd(ind,ind+1)=0;
Be(ind,ind-1)=-2*aty_B;
Bd(ind,ind-1)=2*aty_B;
Ce(ind,ind+1)=0;
Cd(ind,ind+1)=0;
Ce(ind,ind-1)=-2*aty_C;
Cd(ind,ind-1)=2*aty_C;
De(ind,ind+1)=0;
Dd(ind,ind+1)=0;
De(ind,ind-1)=-2*aty_D;
Dd(ind,ind-1)=2*aty_D;
end
%falta ind=nn
Ae(nn,nn-1)=-2*aty_A;
Ad(nn,nn-1)=2*aty_A;
Be(nn,nn-1)=-2*aty_B;
Bd(nn,nn-1)=2*aty_B;
Ce(nn,nn-1)=-2*aty_C;
Cd(nn,nn-1)=2*aty_C;
De(nn,nn-1)=-2*aty_D;
Dd(nn,nn-1)=2*aty_D;
%
%adequar a borda vertical da esquerda 
for i=1:nny
Ae(i,i+nny)=-2*atx_A;
Ad(i,i+nny)=2*atx_A;
Be(i,i+nny)=-2*atx_B;
Bd(i,i+nny)=2*atx_B;
Ce(i,i+nny)=-2*atx_C;
Cd(i,i+nny)=2*atx_C;
De(i,i+nny)=-2*atx_D;
Dd(i,i+nny)=2*atx_D;
end
%
%adequar a borda vertical da direita
for i=1:nny
in=nx*nny+i;
Ae(in,in-nny)=-2*atx_A;
Ad(in,in-nny)=2*atx_A;
Be(in,in-nny)=-2*atx_B;
Bd(in,in-nny)=2*atx_B;
Ce(in,in-nny)=-2*atx_C;
Cd(in,in-nny)=2*atx_C;
De(in,in-nny)=-2*atx_D;
Dd(in,in-nny)=2*atx_D;
end
%
% condições iniciais por espécie
Az=zeros(nn,1); Bz=zeros(nn,1); Cz=zeros(nn,1); Dz=zeros(nn,1);
% espécie A, um círculo no centro do domínio com raio 2
for i=1:nnx
  for j=1:nny
    ind=(i-1)*nny+j;
    if ((i-nnx/2)^2+(j-nny/2)^2 < 2^2)
        Az(ind)=K/7;
    end
  end
end

% espécie B, um círculo com raio 4, centrado em (nnx/2, nny/2)
for i=1:nnx
  for j=1:nny
    ind=(i-1)*nny+j;
    if ((i-nnx/2)^2+(j-nny/2)^2 >= 2^2 && (i-nnx/2)^2+(j-nny/2)^2 < 4^2)
        Bz(ind)=K/7;
    end
  end
end

% espécie C, um círculo com raio 6, centrado em (nnx/2, nny/2)
for i=1:nnx
  for j=1:nny
    ind=(i-1)*nny+j;
    if ((i-nnx/2)^2+(j-nny/2)^2 >= 4^2 && (i-nnx/2)^2+(j-nny/2)^2 < 6^2)
        Cz(ind)=K/7;
    end
  end
end

% espécie D, um círculo com raio 8, centrado em (nnx/2, nny/2)
for i=1:nnx
  for j=1:nny
    ind=(i-1)*nny+j;
    if ((i-nnx/2)^2+(j-nny/2)^2 >= 6^2 && (i-nnx/2)^2+(j-nny/2)^2 < 8^2)
        Dz(ind)=K/7;
    end
  end
end
% Inicialização
A=Az; B=Bz; C=Cz; D=Dz;
%
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resolução repetida do sistema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nAe=Ae; nAd=Ad;   
nBe=Be; nBd=Bd;
nCe=Ce; nCd=Cd;
nDe=De; nDd=Dd;   
%para visualizar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verA=zeros(nny,nnx);
verB=zeros(nny,nnx);
verC=zeros(nny,nnx);
verD=zeros(nny,nnx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nAe=Ae; nAd=Ad;
nBe=Be; nBd=Bd;
nCe=Ce; nCd=Cd;
nDe=De; nDd=Dd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it=1:npt;
%%%%%%%%%%%%%%%%%%%%%%%%%
contador(it)=it; %Incremento para colcoar o crecimento da população D - andiroba
%%%%%%%%%%%%%%%%%%%%%%%%%
Aint=(A+Az)/2;
Bint=(B+Bz)/2;
Cint=(C+Cz)/2;
Dint=(D+Dz)/2;
for k=1:nn
nAe(k,k)=Ae(k,k)+(lambda_A/(K))*(Aint(k)+gama_A*Bint(k)+sigma_A*Cint(k)+delta_A*Dint(k))*dt/2;
nAd(k,k)=Ad(k,k)-(lambda_A/(K))*(Aint(k)+gama_A*Bint(k)+sigma_A*Cint(k)+delta_A*Dint(k))*dt/2;
end
A=nAe\(nAd*Az);
Aint=(A+Az)/2;
%%%%%%%%%%%%
for k=1:nn
nBe(k,k)=Be(k,k)+(lambda_B/(K))*((Bint(k)+gama_B*Aint(k)+sigma_B*Cint(k)+delta_B*Dint(k)))*dt/2;
nBd(k,k)=Bd(k,k)-(lambda_B/(K))*((Bint(k)+gama_B*Aint(k)+sigma_B*Cint(k)+delta_B*Dint(k)))*dt/2;
end
B=nBe\(nBd*Bz);
Bint=(B+Bz)/2;
%%%%%%%%%%%%
for k=1:nn
nCe(k,k)=Ce(k,k)+(lambda_C/(K))*((Cint(k)+gama_C*Bint(k)+sigma_C*Aint(k)+delta_C*Dint(k)))*dt/2;
nCd(k,k)=Cd(k,k)-(lambda_C/(K))*((Cint(k)+gama_C*Bint(k)+sigma_C*Aint(k)+delta_C*Dint(k)))*dt/2;
end
C=nCe\(nCd*Cz);
Cint=(C+Cz)/2;
%%%%%%%%%%%%%%%%%%%
for k=1:nn
nDe(k,k)=De(k,k)+(lambda_D)*Dint(k)/((K + gama_D*Aint(k)+sigma_D*Bint(k)+delta_D*Cint(k)))*dt/2;
nDd(k,k)=Dd(k,k)-(lambda_D)*Dint(k)/((K + gama_D*Aint(k)+sigma_D*Bint(k)+delta_D*Cint(k)))*dt/2;
end
D=nDe\(nDd*Dz);
Dint=(D+Dz)/2;
%%%%%%%%%%%%
%%% Final (real)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Az=A; Bz=B; Cz=C; Dz=D; 
Alt(it)=A(1);
Blt(it)=B(100);
Clt(it)=C(480);
Dlt(it)=D(300);
end
for i=1:nny
 for j=1:nnx
 ind=(j-1)*nny+i;
 verA(nny+1-i,j)=A(ind);
 end
end
 for i=1:nny
 for j=1:nnx
 ind=(j-1)*nny+i;
 verB(nny+1-i,j)=B(ind);
 end
 end
 for i=1:nny
 for j=1:nnx
 ind=(j-1)*nny+i;
 verC(nny+1-i,j)=C(ind);
 end
 end
 for i=1:nny
 for j=1:nnx
 ind=(j-1)*nny+i;
 verD(nny+1-i,j)=D(ind);
 end
end
 
%%%%%%%%%%%%%%%%%%%%%%
%figure(1)
 %hold on
  %surf(verA ), grid on
  %surf(verB), grid on
  %surf(verC), grid on
  %surf(verD), grid on
 %hold off
 %pause(0.7)
 %figure(2)
%hold on
% plot(Alt,'dr')
% plot(Blt,'pk')
% plot(Clt,'sb')
% plot(Dlt,'vg')
 %plot(contador(5:300),Dlt(5:300),'vg')
%hold off
 %%%%%%%%%%%%%%%%%%%%%%
figure(6)
 hold on
subplot(3,2,1)
surf(verA)
title('SPECIES A - GOIABA DE ANTA')
xlabel("LENGTH")
ylabel("WIDTH")
axis([0 40 0 30 0 20])
zlim([0 20])
subplot(3,2,2)
surf(verB)
title('SPECIES B - LACRE')
xlabel("LENGTH")
ylabel("WIDTH")
axis([0 40 0 30 0 20])
zlim([0 20])

subplot(3,2,3)
surf(verC)
title('SPECIES C - UCUUBA')
xlabel("LENGTH")
ylabel("WIDTH")
axis([0 40 0 30 0 20])
zlim([0 20])

subplot(3,2,4)
surf(verD)
title('SPECIES D - ANDIROBA')
xlabel("LENGTH")
ylabel("WIDTH")
axis([0 40 0 30 0 20])
zlim([0 20])
 hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %figure(2)
%hold on
  %caxis([0 K]) % Atribui a cor azul ao valor 0 e a cor vermelha ao valor K
  %axis([1 21 1 17 0 33.43]) % Fixa os eixos
  %contour(verA,6,'dr') %,grid
  %contour(verB,6,'pk') %,grid
  %contour(verC,6,'sb') %,grid
  %contour(verD,6,'vg') %,grid
 %hold off
%figure (3)
 %subplot(2,2,1),plot(Alt,'linewidth',2)
 %subplot(2,2,2),plot(Blt,'linewidth',2)
 %subplot(2,2,3),plot(Clt,'linewidth',2)
 %subplot(2,2,4),plot(Dlt,'linewidth',2)
 figure(1);
verA_norm = verA / max(abs(verA(:)));
%subplot(3,2,1);
surf(verA_norm*20)
%zlim([0 20])
%axis([0 40 0 30 0 8]) % defina os limites do eixo Z entre 0 e 1
figure(2);
verB_norm = verB / max(abs(verB(:)));
%subplot(3,2,2); 
surf(verB_norm*20)
%axis([0 40 0 30 0 8]) % defina os limites do eixo Z entre 0 e 1
figure(3);
verC_norm = verC / max(abs(verC(:)));
%subplot(3,2,3);
surf(verC_norm*20)
%axis([0 40 0 30 0 8]) % defina os limites do eixo Z entre 0 e 1
figure(4);
verD_norm = verD / max(abs(verD(:)));
%subplot(3,2,4);
surf(verD_norm*20)
%axis([0 40 0 30 0 8]) % defina os limites do eixo Z entre 0 e 1
 