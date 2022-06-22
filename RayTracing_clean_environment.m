close all;
clc;
clear all;
tic
%% Dados de Entrada
c0 = 3*10^8;
f = 1*10^9;
W= 2*pi*f;
lambda = c0/f;
B = 2*pi/lambda;
ti = 0;
t = ti;
Ns = 4; % numero de Superficies
Nr = 2; %número de Reflexões
Nmax = (Ns-1)^(Nr);
Elementos = 4;
aux = 1;

%Dados da condição de contorno
cond = 4.73*10^-2;
e0 = 8.8541*10^-12;
er = 7;
im = sqrt(-1);
ec = e0*er - im*((cond)/W);
%ec = er*(e0 - im*((cond)/(W*e0)));
ur = 1;
u0 = 4*pi*10^-7;
uc = ur*u0;
%Ninc = sqrt(uc/ec);
Ninc = sqrt((im*W*ur*u0)/(cond +(im*W*e0*er)));
%N0 = sqrt(u0/e0);
N0 = 120*pi;
alfamat = W*sqrt((uc*ec/2)*[sqrt(1+((cond/(W*ec))^2)-1)]);
betamat = W*sqrt((uc*ec/2)*[sqrt(1+((cond/(W*ec))^2)+1)]);
Yc = alfamat + j*betamat;
d = 0.2;

% matrizes de inicialização
Qx = zeros(Ns,Elementos);
Qy = zeros(Ns,Elementos);
Ix = zeros(Ns,Elementos);
Iy = zeros(Ns,Elementos);
alfa = zeros(Ns,Elementos);
w = zeros(Ns,Elementos);

%% Pontos da área de cálculo
% PontoS1
xs1 = 1;
ys1 = 1;
S1 = [xs1,ys1]; %Ponto correspondente
%Ponto S2
xs2 = 1;
ys2 = 4;
S2 = [xs2,ys2]; %Ponto correspondente
%Ponto S3
xs3 = 4;
ys3 = 4;
S3 = [xs3,ys3]; %Ponto correspondente
%Ponto S4
xs4 = 4;
ys4 = 1;
S4 = [xs4,ys4]; %Ponto correspondente

Ptx = [xs1; xs2; xs3; xs4];
Pty = [ys1; ys2; ys3; ys4];

%%Criação de segmentos da área de cálculo
Segm1 = line([xs1 xs2], [ys1 ys2],'Color','black','LineStyle','-'); %Segmento 1
Segm2 = line([xs2 xs3], [ys2 ys3],'Color','black','LineStyle','-'); %Segmento 2
Segm3 = line([xs3 xs4], [ys3 ys4],'Color','black','LineStyle','-'); %Segmento 3
Segm4 = line([xs4 xs1], [ys4 ys1],'Color','black','LineStyle','-'); %Segmento 4
 xlim([0 5]);
 ylim([0 5]);

%% Localização da antena receptora e transmissora

%Coordenada ponto de transmissão
xTi = 2.5;
yTi = 2.5;
xtf = 4;
hold on
plot(xTi,yTi,'k*');
txt1 = '\ Ti';
plot(xtf,yTi,'k*');
txt2 = '\ Ty';
text(xTi,yTi,txt1,'HorizontalAlignment','left');
text(xtf,yTi,txt2,'HorizontalAlignment','left');
hold off

%Coordenada ponto de recepção
s0 = distance(xTi,yTi,xs4,yTi);
points = s0/100;
xTx = [xTi+points:points:xs4];
N = length(xTx);
yTx = 2.5;
hold on
plot(xTx,yTx,'--');
hold off

%% Cálculos das imagens

 for z = 1:N
    for i = 1:Ns
       j = (i+1);
       if i == Ns
        j = 1;
       end   
            ind = 1;
%           if  (Ptx(j,aux) >= Ptx(i,aux)) || (Pty(j,aux) >= Pty(i,aux))   
             alfa(i,ind) = [(xTi-Ptx(i,aux))*(Ptx(j,aux)-Ptx(i,aux))+(yTi-Pty(i,aux))*(Pty(j,aux)-Pty(i,aux))]/[(Ptx(j,aux)-Ptx(i,aux))^2+(Pty(j,aux)-Pty(i,aux))^2];
             Qx(i,ind) = Ptx(i,aux) + alfa(i,aux)*(Ptx(j,aux)-Ptx(i,aux));
             Qy(i,ind) = Pty(i,aux) + alfa(i,aux)*(Pty(j,aux)-Pty(i,aux)); 
% 
%          else
%              alfa(i,ind) = [(xTi-Ptx(j,aux))*(Ptx(i,aux)-Ptx(j,aux))+(yTi-Pty(j,aux))*(Pty(i,aux)-Pty(j,aux))]/[(Ptx(i,aux)-Ptx(j,aux))^2+(Pty(i,aux)-Pty(j,aux))^2];
%              Qx(i,ind) = Ptx(j,aux) + alfa(i,aux)*(Ptx(i,aux)-Ptx(j,aux));
%              Qy(i,ind) = Pty(j,aux) + alfa(i,aux)*(Pty(i,aux)-Pty(j,aux));
%          end
         
         if (abs(alfa(i,ind)) <= 1 && abs(alfa(i,ind)) >= 0)
             Ix(i,ind) = 2*Qx(i,aux) - xTi;
             Iy(i,ind) = 2*Qy(i,aux) - yTi;
         else
              Ix(i,ind) = NaN;
              Iy(i,ind) = NaN;
         end
%%%%%%%%%%%%%%%%%%%%%%%Coodernada da imagem do ponto Ti em relação ao segmento %i
             
             
% hold on
% plot(Ix(4,ind),Iy(4,ind),'.');
% hold off            
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cálculo dos pontos de reflexões   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            k(i,aux) = (Pty(j,aux)-Pty(i,aux))/(Ptx(j,aux)-Ptx(i,aux));
            k(~isfinite(k))=0
            w(i,z) = ((yTx-Iy(i,aux))/(xTx(z)-Ix(i,aux)));
            beta(i,z) = (-k(i,aux)*(xTx(z)+Ix(i,aux))+(yTx+Iy(i,aux)) -(2*(Pty(i,aux)-(k(i,aux)*Ptx(i,aux)))))/[(k(i,aux)*((xTx(z)-Ix(i,aux))))-(yTx-Iy(i,aux))];
            beta(isnan(beta))= 0
            beta(~isfinite(beta))= 0;
            ALFA(i,z) = (-w(i,z)*(Ptx(j,aux)+Ptx(i,aux))+(Pty(j,aux)+Pty(i,aux)) - (2*(Iy(i,aux)-(w(i,z)*Ix(i,aux)))))/(w(i,z)*(Ptx(j,aux)-Ptx(i,aux))-(Pty(j,aux)-Pty(i,aux)));
            Px(i,z) = 1/2*[(xTx(z) +Ix(i,aux))+ beta(i,z)*(xTx(z)-Ix(i,aux))];
            Px(~isfinite(Px))= 0;
            Py(i,z) = 1/2*[(yTx+Iy(i,aux))+ beta(i,z)*(yTx-Iy(i,aux))];
            Py(~isfinite(Py))= 0;
            Py(isnan(Py))= 0;
            Sx(i,z) = 1/2*((Ptx(j,aux)+Ptx(i,aux))+ ALFA(i,z)*(Ptx(j,aux)-Ptx(i,aux)));
            Sy(i,z) = 1/2*((Pty(j,aux)+Pty(i,aux))+ ALFA(i,z)*(Pty(j,aux)-Pty(i,aux)));
            
             if (abs(beta(i,z)) >= 1 )
                Px(i,z) = NaN;
                Py(i,z) = NaN;
             end
                if (abs(ALFA(i,z)) >= 1  )  
                Sx(i,z) = NaN;
                Sy(i,z) = NaN;
                end
            if (Px(i,z) == Sx(i,z) & Py(i,z)== Sy(i,z) & abs(beta(i,z)) <= 1 &  abs(ALFA(i,z)) <= 1)
                Qx(i,z) = Px(i,z);
                Qy(i,z) = Py(i,z);
            else
                Qx(i,z) = NaN;
                Qy(i,z) = NaN;
            end
         
            
hold on
% plot(Px(i,z), Py(i,z),'.'); 
%plot(Sx(i,z), Sy(i,z),'.'); 
% plot(Qx(:,:), Qy(:,:),'.');
% plot(Ix(4,aux),Iy(4,aux),'.');
 hold off

    end
 end   



 
%% Calculo da imagem para 2reflexão
if Nr > 1 
for z = 1:N
    aux = 1;
 for i = 1:Ns
    for m = 2:Elementos
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1º SEGMENTO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       if i == 1
       p = m;
       n = m + 1;
       if p == 4
        n = 1;
       end
         if  (Ptx(n,aux) >= Ptx(p,aux)) && (Pty(n,aux) >= Pty(p,aux))
             alfa(i,m) = ((Ix(i,1)-Ptx(p,aux))*(Ptx(n,aux)-Ptx(p,aux))+(Iy(i,1)-Pty(p,aux))*(Pty(n,aux)-Pty(p,aux)))/((Ptx(n,aux)-Ptx(p,aux))^2+(Pty(n,aux)-Pty(p,aux))^2);  
             Qx(i,m) = Ptx(p,aux) + alfa(i,m)*(Ptx(n,aux)-Ptx(p,aux));
             Qy(i,m) = Pty(p,aux) + alfa(i,m)*(Pty(n,aux)-Pty(p,aux));
       else 
             alfa(i,m) = ((Ix(i,1)-Ptx(n,aux))*(Ptx(p,aux)-Ptx(n,aux))+(Iy(i,1)-Pty(n,aux))*(Pty(p,aux)-Pty(n,aux)))/((Ptx(p,aux)-Ptx(n,aux))^2+(Pty(p,aux)-Pty(n,aux))^2);
             Qx(i,m) = Ptx(n,aux) + alfa(i,m)*(Ptx(p,aux)-Ptx(n,aux));
             Qy(i,m) = Pty(n,aux) + alfa(i,m)*(Pty(p,aux)-Pty(n,aux));
       end
       
       if   (abs(alfa(i,m)) <= 1 && abs(alfa(i,m)) >= 0)
           
             Ix(i,m) = 2*Qx(i,m) - Ix(i,1);
             Iy(i,m) = 2*Qy(i,m) - Iy(i,1);
       else
              Ix(i,m) = NaN;
              Iy(i,m) = NaN;
         end
             
       
     
    
% hold on
% plot(Ix(i,m),Iy(i,m),'.');
% hold off      

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%% calculos dos pontos de reflexões DA 1 IMAGEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da imagem de segunda ordem para ponto de recepção %%%%%%%%%%%%%%%%%%%
            k(i,m) = (Pty(n,aux)-Pty(p,aux))/(Ptx(n,aux)-Ptx(p,aux));
            k(~isfinite(k))= 0;
            w(i,4*z-4+m) = ((yTx-Iy(i,m))/(xTx(z)-Ix(i,m)));
            Beta(i,4*z-4+m) = ((-k(i,m)*(xTx(z)+Ix(i,m)))+(yTx+Iy(i,m)) -(2*(Pty(p,aux)-(k(i,m)*Ptx(p,aux)))))/((k(i,m)*((xTx(z)-Ix(i,m))))-(yTx-Iy(i,m)));
            Beta(isnan(Beta))= 0;
            Beta(~isfinite(Beta))= 0;
            Alfa(i,4*z-4+m) = ((-w(i,4*z-4+m)*(Ptx(n,aux)+Ptx(p,aux)))+(Pty(n,aux)+Pty(p,aux)) - (2*(Iy(i,m)-(w(i,4*z-4+m)*Ix(i,m)))))/((w(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)))-(Pty(n,aux)-Pty(p,aux)));    
            Pontox(i,4*z-4+m) = 1/2*((xTx(z) + Ix(i,m))+ Beta(i,4*z-4+m)*(xTx(z)-Ix(i,m)));
            Pontox(~isfinite(Pontox))= 0;
            Pontoy(i,4*z-4+m) = 1/2*((yTx + Iy(i,m))+ Beta(i,4*z-4+m)*(yTx-Iy(i,m)));
            Rx(i,4*z-4+m) = 1/2*((Ptx(n,aux)+Ptx(p,aux))+ Alfa(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)));
            Ry(i,4*z-4+m) = 1/2*((Pty(n,aux)+Pty(p,aux))+ Alfa(i,4*z-4+m)*(Pty(n,aux)-Pty(p,aux)));

               
             if (abs(Alfa(i,4*z-4+m)) >= 1 ) 
                Rx(i,4*z-4+m) = NaN;
                Ry(i,4*z-4+m) = NaN;
             end
             if (abs(Beta(i,4*z-4+m)) >= 1 )
                Pontoy(i,4*z-4+m) = NaN;
                Pontox(i,4*z-4+m) = NaN;
             end
            
             if (Pontox(i,4*z-4+m) == Rx(i,4*z-4+m) & Pontoy(i,4*z-4+m)== Ry(i,4*z-4+m))
                 Qx1(i,4*z-4+m) = Pontox(i,4*z-4+m);
                 Qy1(i,4*z-4+m) = Pontoy(i,4*z-4+m);
             else
                 Qx1(i,4*z-4+m) = NaN;
                 Qy1(i,4*z-4+m) = NaN;
             end

   hold on
   %plot(Qx1(i,4*z-4+m),Qy1(i,4*z-4+m),'.');
  % plot(Rx(1,4*z-4+m),Ry(1,4*z-4+m),'.');
   %plot(Pontox(1,4*z-4+m), Pontoy(1,4*z-4+m),'.');% Ponto de reflexão em cada segmento
   hold off 
 
 
 px = 1;
 py = 2;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da imagem de Primeira ordem para ponto de transmissao %%%%%%%%%%%%%%%%%%%
            k1(i,m) = (Pty(py,aux)-Pty(px,aux))/(Ptx(py,aux)-Ptx(px,aux));
            k1(~isfinite(k1))= 0;
            w1(i,4*z-4+m) = ((Iy(i,1)- Ry(i,4*z-4+m))/(Ix(i,1)-Rx(i,4*z-4+m)));
            Beta1(i,4*z-4+m) = ((-k1(i,m)*(Ix(i,1)+Rx(i,4*z-4+m)))+(Iy(i,1)+Ry(i,4*z-4+m)) -(2*(Pty(px,aux)-(k1(i,m)*Ptx(px,aux)))))/((k1(i,m)*(Ix(i,1)-Rx(i,4*z-4+m)))-(Iy(i,1)-Ry(i,4*z-4+m)));
            Beta1(isnan(Beta1))= 0;
            Beta1(~isfinite(Beta1))= 0;
            Alfa1(i,4*z-4+m) = ((-w1(i,4*z-4+m)*(Ptx(py,aux)+Ptx(px,aux)))+(Pty(py,aux)+Pty(px,aux))-(2*(Ry(i,4*z-4+m)-(w1(i,4*z-4+m)*Rx(i,4*z-4+m)))))/((w1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)))-(Pty(py,aux)-Pty(px,aux)));
            Pontox1(i,4*z-4+m) = 1/2*((Ix(i,1)+Rx(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Ix(i,1)-Rx(i,4*z-4+m)));
            Pontox1(~isfinite(Pontox1))= 0;
            Pontoy1(i,4*z-4+m) = 1/2*((Iy(i,1)+Ry(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Iy(i,1)-Ry(i,4*z-4+m)));
            Rx1(i,4*z-4+m) = 1/2*((Ptx(py,aux)+Ptx(px,aux))+ Alfa1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)));
            Ry1(i,4*z-4+m) = 1/2*((Pty(py,aux)+Pty(px,aux))+ Alfa1(i,4*z-4+m)*(Pty(py,aux)-Pty(px,aux)));
            
            if (abs(Alfa1(i,4*z-4+m)) >= 1 )
                Rx1(i,4*z-4+m) = NaN;
                Ry1(i,4*z-4+m) = NaN;
            end
            if (abs(Beta1(i,4*z-4+m)) >= 1)
                Pontoy1(i,4*z-4+m) = NaN;
                Pontox1(i,4*z-4+m) = NaN;
             end
            
             if (Pontox1(i,4*z-4+m) == Rx1(i,4*z-4+m) & Pontoy1(i,4*z-4+m) == Ry1(i,4*z-4+m))
                 Qx11(i,4*z-4+m) = Pontox1(i,4*z-4+m);
                 Qy11(i,4*z-4+m) = Pontoy1(i,4*z-4+m);
                 
             else
                 Qx11(i,4*z-4+m) = NaN;
                 Qy11(i,4*z-4+m) = NaN;
             end
             
   hold on
   %plot(Qx11(1,4*z-4+m),Qy11(1,4*z-4+m),'.');
   %plot(Pontox1(1,4*z-4+m),Pontoy1(1,4*z-4+m),'.'); % Ponto de reflexão em cada segmento
  % plot(Rx1(1,4*z-4+m),Ry1(1,4*z-4+m),'.');
   hold off
 
 %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%%
          
           
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2º SEGMENTO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    elseif i == 2
            p = m + 1;
            n = p + 1;
       if p == 4
            n = 1;
       elseif p == 5
           p = 1;
           n = 2;
       end

             alfa(i,m) = ((Ix(i,1)-Ptx(p,aux))*(Ptx(n,aux)-Ptx(p,aux))+(Iy(i,1)-Pty(p,aux))*(Pty(n,aux)-Pty(p,aux)))/((Ptx(n,aux)-Ptx(p,aux))^2+(Pty(n,aux)-Pty(p,aux))^2);  
             Qx(i,m) = Ptx(p,aux) + alfa(i,m)*(Ptx(n,aux)-Ptx(p,aux));
             Qy(i,m) = Pty(p,aux) + alfa(i,m)*(Pty(n,aux)-Pty(p,aux));

           
       if (abs(alfa(i,m)) <= 1 && abs(alfa(i,m)) >= 0)
           
             Ix(i,m) = 2*Qx(i,m) - Ix(i,1);
             Iy(i,m) = 2*Qy(i,m) - Iy(i,1);
       else
              Ix(i,m) = NaN;
              Iy(i,m) = NaN;
        end
             
          

       
% if abs(alfa(i,m)) <= 1
% hold on
% plot(Ix(2,m),Iy(2,m),'.');
% hold off      
% end  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculos dos pontos de reflexões DA 2 IMAGEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            k(i,m) = (Pty(n,aux)-Pty(p,aux))/(Ptx(n,aux)-Ptx(p,aux));
            k(~isfinite(k))= 0;
            w(i,4*z-4+m) = ((yTx-Iy(i,m))/(xTx(z)-Ix(i,m)));
            Beta(i,4*z-4+m) = ((-k(i,m)*(xTx(z)+Ix(i,m)))+(yTx+Iy(i,m)) -(2*(Pty(p,aux)-(k(i,m)*Ptx(p,aux)))))/((k(i,m)*((xTx(z)-Ix(i,m))))-(yTx-Iy(i,m)));
            Beta(isnan(Beta))= 0;
            Beta(~isfinite(Beta))= 0;
            Alfa(i,4*z-4+m) = (-w(i,4*z-4+m)*(Ptx(n,aux)+Ptx(p,aux))+(Pty(n,aux)+Pty(p,aux)) - (2*(Iy(i,m)-(w(i,4*z-4+m)*Ix(i,m)))))/(w(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux))-(Pty(n,aux)-Pty(p,aux)));
              
            Pontox(i,4*z-4+m) = 1/2*[(xTx(z) + Ix(i,m))+ Beta(i,4*z-4+m)*(xTx(z)-Ix(i,m))];
            Pontox(~isfinite(Pontox))= 0;
            Pontoy(i,4*z-4+m) = 1/2*[(yTx+Iy(i,m))+ Beta(i,4*z-4+m)*(yTx-Iy(i,m))];
            Rx(i,4*z-4+m) = 1/2*((Ptx(n,aux)+Ptx(p,aux))+ Alfa(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)));
            Ry(i,4*z-4+m) = 1/2*((Pty(n,aux)+Pty(p,aux))+ Alfa(i,4*z-4+m)*(Pty(n,aux)-Pty(p,aux)));
            
            if (abs(Alfa(i,4*z-4+m)) >= 1 ) 
                Rx(i,4*z-4+m) = NaN;
                Ry(i,4*z-4+m) = NaN;
             end
             if (abs(Alfa(i,4*z-4+m)) >= 1 && abs(Beta(i,4*z-4+m)) >= 1)
                Pontoy(i,4*z-4+m) = NaN;
                Pontox(i,4*z-4+m) = NaN;
            end
            
            if (Pontox(i,4*z-4+m) == Rx(i,4*z-4+m) & Pontoy(i,4*z-4+m)== Ry(i,4*z-4+m))
                 Qx1(i,4*z-4+m) = Pontox(i,4*z-4+m);
                 Qy1(i,4*z-4+m) = Pontoy(i,4*z-4+m);
             else
                 Qx1(i,4*z-4+m) = NaN;
                 Qy1(i,4*z-4+m) = NaN;
             end

   hold on
   %plot(Qx1(2,4*z-4+m),Qy1(2,4*z-4+m),'.');
  % plot(Rx(2,4*z-4+m),Ry(2,4*z-4+m),'.');
   %plot(Pontox(1,4*z-4+m), Pontoy(1,4*z-4+m),'.');% Ponto de reflexão em cada segmento
   hold off   
            
 px = 2;
 py = 3;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da imagem de Primeira ordem para ponto de recepção %%%%%%%%%%%%%%%%%%%
 
            k1(i,m) = (Pty(py,aux)-Pty(px,aux))/(Ptx(py,aux)-Ptx(px,aux));
            k1(~isfinite(k1))= 0;
            w1(i,4*z-4+m) = ((Iy(i,1)-Ry(i,4*z-4+m))/(Ix(i,1)-Rx(i,4*z-4+m)));
            Beta1(i,4*z-4+m) = ((-k1(i,m)*(Ix(i,1)+Rx(i,4*z-4+m)))+(Iy(i,1)+Ry(i,4*z-4+m)) -(2*(Pty(px,aux)-(k1(i,m)*Ptx(px,aux)))))/((k1(i,m)*(Ix(i,1)-Rx(i,4*z-4+m)))-(Iy(i,1)-Ry(i,4*z-4+m)));
            Beta1(isnan(Beta1))= 0;
            Beta1(~isfinite(Beta1))= 0;
            Alfa1(i,4*z-4+m) = (-w1(i,4*z-4+m)*(Ptx(py,aux)+Ptx(px,aux))+(Pty(py,aux)+Pty(px,aux)) - (2*(Ry(i,4*z-4+m)-(w1(i,4*z-4+m)*Rx(i,4*z-4+m)))))/((w1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)))-(Pty(py,aux)-Pty(px,aux)));
            Pontox1(i,4*z-4+m) = 1/2*[(Ix(i,1)+Rx(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Ix(i,1)-Rx(i,4*z-4+m))];
            Pontox1(~isfinite(Pontox))= 0;
            Pontoy1(i,4*z-4+m) = 1/2*[(Iy(i,1)+Ry(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Iy(i,1)-Ry(i,4*z-4+m))];
            Rx1(i,4*z-4+m) = 1/2*((Ptx(py,aux)+Ptx(px,aux))+ Alfa1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)));
            Ry1(i,4*z-4+m) = 1/2*((Pty(py,aux)+Pty(px,aux))+ Alfa1(i,4*z-4+m)*(Pty(py,aux)-Pty(px,aux)));
                        
            if (abs(Alfa1(i,4*z-4+m)) >= 1 )
                Rx1(i,4*z-4+m) = NaN;
                Ry1(i,4*z-4+m) = NaN;
            end
            if (abs(Beta1(i,4*z-4+m)) >= 1)
                Pontoy1(i,4*z-4+m) = NaN;
                Pontox1(i,4*z-4+m) = NaN;
             end
            
             if (Pontox1(i,4*z-4+m) == Rx1(i,4*z-4+m) & Pontoy1(i,4*z-4+m) == Ry1(i,4*z-4+m))
                 Qx11(i,4*z-4+m) = Pontox1(i,4*z-4+m);
                 Qy11(i,4*z-4+m) = Pontoy1(i,4*z-4+m);
                 
             else
                 Qx11(i,4*z-4+m) = NaN;
                 Qy11(i,4*z-4+m) = NaN;
             end
             
    hold on
%    %plot(Qx11(1,4*z-4+m),Qy11(1,4*z-4+m),'.');
%    %plot(Pontox1(1,4*z-4+m),Pontoy1(1,4*z-4+m),'.'); % Ponto de reflexão em cada segmento
   % plot(Rx1(i,4*z-4+m),Ry1(i,4*z-4+m),'.');
   hold off
 

 %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3º SEGMENTO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   elseif i == 3
        p = 1;
        n = 4;
      if m == 3
           p = 1;
           n = 2;
      elseif m == 4
           p = 2;
           n = 3;
      end
             
%         if  (Ptx(n,aux) >= Ptx(p,aux)) && (Pty(n,aux) >= Pty(p,aux))
             alfa(i,m) = ((Ix(i,1)-Ptx(p,aux))*(Ptx(n,aux)-Ptx(p,aux))+(Iy(i,1)-Pty(p,aux))*(Pty(n,aux)-Pty(p,aux)))/((Ptx(n,aux)-Ptx(p,aux))^2+(Pty(n,aux)-Pty(p,aux))^2);  
             Qx(i,m) = Ptx(p,aux) + alfa(i,m)*(Ptx(n,aux)-Ptx(p,aux));
             Qy(i,m) = Pty(p,aux) + alfa(i,m)*(Pty(n,aux)-Pty(p,aux));
%        else 
%              alfa(i,m) = ((Ix(i,1)-Ptx(n,aux))*(Ptx(p,aux)-Ptx(n,aux))+(Iy(i,1)-Pty(n,aux))*(Pty(p,aux)-Pty(n,aux)))/((Ptx(p,aux)-Ptx(n,aux))^2+(Pty(p,aux)-Pty(n,aux))^2);
%              Qx(i,m) = Ptx(n,aux) + alfa(i,m)*(Ptx(p,aux)-Ptx(n,aux));
%              Qy(i,m) = Pty(n,aux) + alfa(i,m)*(Pty(p,aux)-Pty(n,aux));
%        end

                       
        if (abs(alfa(i,m)) <= 1  && abs(alfa(i,m)) >= 0)
           
             Ix(i,m) = 2*Qx(i,m) - Ix(i,1);
             Iy(i,m) = 2*Qy(i,m) - Iy(i,1);
       else
              Ix(i,m) = NaN;
              Iy(i,m) = NaN;
        end
             

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculos dos pontos de reflexões DA 3 IMAGEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            k(i,m) = (Pty(n,aux)-Pty(p,aux))/(Ptx(n,aux)-Ptx(p,aux));
            k(~isfinite(k))= 0;
            w(i,4*z-4+m) = ((yTx-Iy(i,m))/(xTx(z)-Ix(i,m)));
            Beta(i,4*z-4+m) = ((-k(i,m)*(xTx(z)+Ix(i,m)))+(yTx+Iy(i,m)) -(2*(Pty(p,aux)-(k(i,m)*Ptx(p,aux)))))/((k(i,m)*((xTx(z)-Ix(i,m))))-(yTx-Iy(i,m)));
            Beta(isnan(Beta))= 0;
            Beta(~isfinite(Beta))= 0;
            Alfa(i,4*z-4+m) = (-w(i,4*z-4+m)*(Ptx(n,aux)+Ptx(p,aux))+(Pty(n,aux)+Pty(p,aux)) - (2*(Iy(i,m)-(w(i,4*z-4+m)*Ix(i,m)))))/((w(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)))-(Pty(n,aux)-Pty(p,aux)));
            Pontox(i,4*z-4+m) = 1/2*[(xTx(z) + Ix(i,m))+ Beta(i,4*z-4+m)*(xTx(z)-Ix(i,m))];
            Pontox(~isfinite(Pontox))= 0;
            Pontoy(i,4*z-4+m) = 1/2*[(yTx+Iy(i,m))+ Beta(i,4*z-4+m)*(yTx-Iy(i,m))];
            Rx(i,4*z-4+m) = 1/2*((Ptx(n,aux)+Ptx(p,aux))+ Alfa(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)));
            Ry(i,4*z-4+m) = 1/2*((Pty(n,aux)+Pty(p,aux))+ Alfa(i,4*z-4+m)*(Pty(n,aux)-Pty(p,aux)));

            
            if (abs(Alfa(i,4*z-4+m)) >= 1) 
                Rx(i,4*z-4+m) = NaN;
                Ry(i,4*z-4+m) = NaN;
             end
             if (abs(Beta(i,4*z-4+m)) >= 1)
                Pontoy(i,4*z-4+m) = NaN;
                Pontox(i,4*z-4+m) = NaN;
            end
            
            if (Pontox(i,4*z-4+m) == Rx(i,4*z-4+m) & Pontoy(i,4*z-4+m)== Ry(i,4*z-4+m))
                 Qx1(i,4*z-4+m) = Pontox(i,4*z-4+m);
                 Qy1(i,4*z-4+m) = Pontoy(i,4*z-4+m);
             else
                 Qx1(i,4*z-4+m) = NaN;
                 Qy1(i,4*z-4+m) = NaN;
             end

                hold on
%                plot(Qx1(3,4*z-4+m),Qy1(3,4*z-4+m),'.');
               % plot(Rx(3,4*z-4+m),Ry(3,4*z-4+m),'.');
%                plot(Pontox(3,4*z-4+m), Pontoy(3,4*z-4+m),'.');% Ponto de reflexão em cada segmento
                hold off      
%               
 
 px = 4;
 py = 3;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da imagem de Primeira ordem para ponto de recepção %%%%%%%%%%%%%%%%%%%
 
            k1(i,m) = (Pty(py,aux)-Pty(px,aux))/(Ptx(py,aux)-Ptx(px,aux));
            k1(~isfinite(k1))= 0;
            w1(i,4*z-4+m) = ((Iy(i,1)-Ry(i,4*z-4+m))/(Ix(i,1)-Rx(i,4*z-4+m)));
            Beta1(i,4*z-4+m) = ((-k1(i,m)*(Ix(i,1)+Rx(i,4*z-4+m)))+(Iy(i,1)+Ry(i,4*z-4+m)) -(2*(Pty(px,aux)-(k1(i,m)*Ptx(px,aux)))))/((k1(i,m)*(Ix(i,1)-Rx(i,4*z-4+m)))-(Iy(i,1)-Ry(i,4*z-4+m)));
            Beta1(isnan(Beta1))= 0;
            Beta1(~isfinite(Beta1))= 0;
            Alfa1(i,4*z-4+m) = (-w1(i,4*z-4+m)*(Ptx(py,aux)+Ptx(px,aux))+(Pty(py,aux)+Pty(px,aux)) - (2*(Ry(i,4*z-4+m)-(w1(i,4*z-4+m)*Rx(i,4*z-4+m)))))/((w1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)))-(Pty(py,aux)-Pty(px,aux)));
            Pontox1(i,4*z-4+m) = 1/2*[(Ix(i,1)+Rx(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Ix(i,1)-Rx(i,4*z-4+m))];
            Pontox1(~isfinite(Pontox1))= 0;
            Pontoy1(i,4*z-4+m) = 1/2*[(Iy(i,1)+Ry(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Iy(i,1)-Ry(i,4*z-4+m))];
            Rx1(i,4*z-4+m) = 1/2*((Ptx(py,aux)+Ptx(px,aux))+ Alfa1(i,4*z-4+m)*(Pty(py,aux)-Ptx(px,aux)));
            Ry1(i,4*z-4+m) = 1/2*((Pty(py,aux)+Pty(px,aux))+ Alfa1(i,4*z-4+m)*(Pty(py,aux)-Ptx(px,aux)));
  
 
                        if (abs(Alfa1(i,4*z-4+m)) >= 1)
                            Rx1(i,4*z-4+m) = NaN;
                            Ry1(i,4*z-4+m) = NaN;
                        end
                        if (abs(Beta1(i,4*z-4+m)) >= 1)
                            Pontoy1(i,4*z-4+m) = NaN;
                            Pontox1(i,4*z-4+m) = NaN;
                         end

                         if (Pontox1(i,4*z-4+m) == Rx1(i,4*z-4+m) & Pontoy1(i,4*z-4+m) == Ry1(i,4*z-4+m))
                             Qx11(i,4*z-4+m) = Pontox1(i,4*z-4+m);
                             Qy11(i,4*z-4+m) = Pontoy1(i,4*z-4+m);

                         else
                             Qx11(i,4*z-4+m) = NaN;
                             Qy11(i,4*z-4+m) = NaN;
                         end

                hold on
%                plot(Qx11(3,4*z-4+m),Qy11(3,4*z-4+m),'.');
%                plot(Pontox1(3,4*z-4+m),Pontoy1(3,4*z-4+m),'.'); % Ponto de reflexão em cada segmento
               %plot(Rx1(3,4*z-4+m),Ry1(3,4*z-4+m),'.');
                hold off

 %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%%        

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4º SEGMENTO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%
   elseif i == 4
        p = 1;
        n = 2;
       if m == 3
        p = 2;
        n = 3;
       elseif m == 4
        p = 3;
        n = 4;
       end
           
       if  (Ptx(n,aux) >= Ptx(p,aux)) && (Pty(n,aux) >= Pty(p,aux))
             alfa(i,m) = ((Ix(i,1)-Ptx(p,aux))*(Ptx(n,aux)-Ptx(p,aux))+(Iy(i,1)-Pty(p,aux))*(Pty(n,aux)-Pty(p,aux)))/((Ptx(n,aux)-Ptx(p,aux))^2+(Pty(n,aux)-Pty(p,aux))^2);  
             Qx(i,m) = Ptx(p,aux) + alfa(i,m)*(Ptx(n,aux)-Ptx(p,aux));
             Qy(i,m) = Pty(p,aux) + alfa(i,m)*(Pty(n,aux)-Pty(p,aux));
       else 
             alfa(i,m) = ((Ix(i,1)-Ptx(n,aux))*(Ptx(p,aux)-Ptx(n,aux))+(Iy(i,1)-Pty(n,aux))*(Pty(p,aux)-Pty(n,aux)))/((Ptx(p,aux)-Ptx(n,aux))^2+(Pty(p,aux)-Pty(n,aux))^2);
             Qx(i,m) = Ptx(n,aux) + alfa(i,m)*(Ptx(p,aux)-Ptx(n,aux));
             Qy(i,m) = Pty(n,aux) + alfa(i,m)*(Pty(p,aux)-Pty(n,aux));
       end

        
             Ix(i,m) = 2*Qx(i,m) - Ix(i,1);
             Iy(i,m) = 2*Qy(i,m) - Iy(i,1);
             
         if (abs(alfa(i,m)) <= 1  && abs(alfa(i,m)) >= 0)
           
             Ix(i,m) = 2*Qx(i,m) - Ix(i,1);
             Iy(i,m) = 2*Qy(i,m) - Iy(i,1);
         else
              Ix(i,m) = NaN;
              Iy(i,m) = NaN;
        end
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da imagem de segunda ordem para ponto de recepção %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            k(i,m) = (Pty(n,aux)-Pty(p,aux))/(Ptx(n,aux)-Ptx(p,aux));
            k(~isfinite(k))= 0;
            w(i,4*z-4+m) = ((yTx-Iy(i,m))/(xTx(z)-Ix(i,m)));
            Beta(i,4*z-4+m) = ((-k(i,m)*(xTx(z)+Ix(i,m)))+(yTx+Iy(i,m)) -(2*(Pty(p,aux)-(k(i,m)*Ptx(p,aux)))))/((k(i,m)*((xTx(z)-Ix(i,m))))-(yTx-Iy(i,m)));
            Beta(isnan(Beta))= 0;
            Beta(~isfinite(Beta))= 0;
            Alfa(i,4*z-4+m) = (-w(i,4*z-4+m)*(Ptx(n,aux)+Ptx(p,aux))+(Pty(n,aux)+Pty(p,aux)) - (2*(Iy(i,m)-(w(i,4*z-4+m)*Ix(i,m)))))/(w(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux))-(Pty(n,aux)-Pty(p,aux)));
              
            Pontox(i,4*z-4+m) = 1/2*[(xTx(z) + Ix(i,m))+ Beta(i,4*z-4+m)*(xTx(z)-Ix(i,m))];
            Pontox(~isfinite(Pontox))= 0;
            Pontoy(i,4*z-4+m) = 1/2*[(yTx+Iy(i,m))+ Beta(i,4*z-4+m)*(yTx-Iy(i,m))];
            Rx(i,4*z-4+m) = 1/2*((Ptx(n,aux)+Ptx(p,aux))+ Alfa(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)));
            Ry(i,4*z-4+m) = 1/2*((Pty(n,aux)+Pty(p,aux))+ Alfa(i,4*z-4+m)*(Pty(n,aux)-Pty(p,aux)));
               
            
            if (abs(Alfa(i,4*z-4+m)) >= 1 ) 
                Rx(i,4*z-4+m) = NaN;
                Ry(i,4*z-4+m) = NaN;
             end
             if (abs(Beta(i,4*z-4+m)) >= 1)
                Pontoy(i,4*z-4+m) = NaN;
                Pontox(i,4*z-4+m) = NaN;
            end
            
            if (Pontox(i,4*z-4+m) == Rx(i,4*z-4+m) & Pontoy(i,4*z-4+m)== Ry(i,4*z-4+m))
                 Qx1(i,4*z-4+m) = Pontox(i,4*z-4+m);
                 Qy1(i,4*z-4+m) = Pontoy(i,4*z-4+m);
             else
                 Qx1(i,4*z-4+m) = NaN;
                 Qy1(i,4*z-4+m) = NaN;
             end

                hold on
%                plot(Qx1(4,4*z-4+m),Qy1(3,4*z-4+m),'.');
                %plot(Rx(4,4*z-4+m),Ry(4,4*z-4+m),'.');
%                plot(Pontox(4,4*z-4+m), Pontoy(3,4*z-4+m),'.');% Ponto de reflexão em cada segmento
                hold off     
 
 px = 1;
 py = 4;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da imagem de Primeira ordem para ponto de recepção %%%%%%%%%%%%%%%%%%%
 
            k1(i,m) = (Pty(py,aux)-Pty(px,aux))/(Ptx(py,aux)-Ptx(px,aux));
            k1(~isfinite(k1))= 0;
            w1(i,4*z-4+m) = ((Iy(i,1)-Ry(i,4*z-4+m))/(Ix(i,1)-Rx(i,4*z-4+m)));
            Beta1(i,4*z-4+m) = ((-k1(i,m)*(Ix(i,1)+Rx(i,4*z-4+m)))+(Iy(i,1)+Ry(i,4*z-4+m)) -(2*(Pty(px,aux)-(k1(i,m)*Ptx(px,aux)))))/((k1(i,m)*(Ix(i,1)-Rx(i,4*z-4+m)))-(Iy(i,1)-Ry(i,4*z-4+m)));
            Beta1(isnan(Beta1))= 0;
            Beta1(~isfinite(Beta1))= 0;
            Alfa1(i,4*z-4+m) = (-w1(i,4*z-4+m)*(Ptx(py,aux)+Ptx(px,aux))+(Pty(py,aux)+Pty(px,aux)) - (2*(Ry(i,4*z-4+m)-(w1(i,4*z-4+m)*Rx(i,4*z-4+m)))))/((w1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)))-(Pty(py,aux)-Pty(px,aux)));
            Pontox1(i,4*z-4+m) = 1/2*[(Ix(i,1)+Rx(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Ix(i,1)-Rx(i,4*z-4+m))];
            Pontox1(~isfinite(Pontox))= 0;
            Pontoy1(i,4*z-4+m) = 1/2*[(Iy(i,1)+Ry(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Iy(i,1)-Ry(i,4*z-4+m))];
            Rx1(i,4*z-4+m) = 1/2*((Ptx(py,aux)+Ptx(px,aux))+ Alfa1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)));
            Ry1(i,4*z-4+m) = 1/2*((Pty(py,aux)+Pty(px,aux))+ Alfa1(i,4*z-4+m)*(Pty(py,aux)-Pty(px,aux)));
  
 

                        if (abs(Alfa1(i,4*z-4+m)) >= 1 )
                            Rx1(i,4*z-4+m) = NaN;
                            Ry1(i,4*z-4+m) = NaN;
                        end
                        if (abs(Beta1(i,4*z-4+m)) >= 1)
                            Pontoy1(i,4*z-4+m) = NaN;
                            Pontox1(i,4*z-4+m) = NaN;
                         end

                         if (Pontox1(i,4*z-4+m) == Rx1(i,4*z-4+m) & Pontoy1(i,4*z-4+m) == Ry1(i,4*z-4+m))
                             Qx11(i,4*z-4+m) = Pontox1(i,4*z-4+m);
                             Qy11(i,4*z-4+m) = Pontoy1(i,4*z-4+m);

                         else
                             Qx11(i,4*z-4+m) = NaN;
                             Qy11(i,4*z-4+m) = NaN;
                         end

               hold on
               %plot(Qx11(4,4*z-4+m),Qy11(3,4*z-4+m),'.');
               %plot(Pontox1(4,4*z-4+m),Pontoy1(4,4*z-4+m),'.'); % Ponto de reflexão em cada segmento
               %plot(Rx1(4,4*z-4+m),Ry1(4,4*z-4+m),'.');
               hold off 

       
       
            end  
        end
    end
  end 
end
%% Para 3 reflexões 
if Nr > 2 
  pt1 = [1 2 3 4];
  for z = 1:N
    for i = 1:Ns
        G(i,1) = i;
         G(i,2:4) = find(pt1 ~= i)     
            for m = 2:4
                N3 = find(pt1 ~= G(i,m))
                     for q = 2*(m-2)+(m-1) : 2*(m-2)+ m + 1
                         G1(i,q) = G(i,m);
                         G2(i,q) = G(i,1);
                         for j = 1:3
                            if (rem(q,j) == 0)
                              point(i,q) = N3(3);
                            elseif (rem((q+1),3) == 0)
                               point(i,q) = N3(2);
                            else
                              point(i,q) = N3(1);
                            end
                            
             p1(i,q) = point(i,q);
             p2(i,q) = point(i,q)+1;
             if p2(i,q) == 5
                 p2(i,q)  = 1;
             end
             
             alfa3(i,q) = ((Ix(i,m)-Ptx(p1(i,q),aux))*(Ptx(p2(i,q),aux)-Ptx(p1(i,q),aux))+(Iy(i,m)-Pty(p1(i,q),aux))*(Pty(p2(i,q),aux)-Pty(p1(i,q),aux)))/((Ptx(p2(i,q),aux)-Ptx(p1(i,q),aux))^2+(Pty(p2(i,q),aux)-Pty(p1(i,q),aux))^2);  
             Qx3(i,q) = Ptx(p1(i,q),aux) + alfa3(i,q)*(Ptx(p2(i,q),aux)-Ptx(p1(i,q),aux));
             Qy3(i,q) = Pty(p1(i,q),aux) + alfa3(i,q)*(Pty(p2(i,q),aux)-Pty(p1(i,q),aux));

        if (abs(alfa3(i,q)) <= 1  && abs(alfa3(i,q)) >= 0)
           
             Ix3(i,q) = 2*Qx3(i,q) - Ix(i,m);
             Iy3(i,q) = 2*Qy3(i,q) - Iy(i,m);
        else
              Ix3(i,q) = NaN;
              Iy3(i,q) = NaN;
             
        end
                         end
 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da imagem de terceira ordem para ponto de recepção %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aux = 1
            k3(i,q) = (Pty(p2(i,q),aux)-Pty(p1(i,q),aux))/(Ptx(p2(i,q),aux)-Ptx(p1(i,q),aux));
            k3(~isfinite(k3))= 0;
            w3(i,9*z-9+q) = ((Iy3(i,q)-yTx)/(Ix3(i,q)-xTx(z)));
            Beta3(i,9*z-9+q) = ((-k3(i,q)*(xTx(z)+Ix3(i,q)))+(yTx+Iy3(i,q)) -(2*(Pty(p1(i,q),aux)-(k3(i,q)*Ptx(p1(i,q),aux)))))/((k3(i,q)*((Ix3(i,q)-xTx(z))))-(Iy3(i,q)-yTx));
            Beta3(isnan(Beta3))= 0;
            Beta3(~isfinite(Beta3))= 0;
            Alfa3(i,9*z-9+q) = (-w3(i,9*z-9+q)*(Ptx(p2(i,q),aux)+Ptx(p1(i,q),aux))+(Pty(p2(i,q),aux)+Pty(p1(i,q),aux)) - (2*(yTx-(w3(i,9*z-9+q)*xTx(z)))))/((w3(i,9*z-9+q)*(Ptx(p2(i,q),aux)-Ptx(p1(i,q),aux)))-(Pty(p2(i,q),aux)-Pty(p1(i,q),aux)));
            Pontox3(i,9*z-9+q) = 1/2*[(xTx(z) + Ix3(i,q))+ Beta3(i,9*z-9+q)*(Ix3(i,q)-xTx(z))];
            Pontox3(~isfinite(Pontox3))= 0;
            Pontoy3(i,9*z-9+q) = 1/2*[(yTx+Iy3(i,q))+ Beta3(i,9*z-9+q)*(Iy3(i,q)-yTx)];
            Rx3(i,9*z-9+q) = 1/2*((Ptx(p2(i,q),aux)+Ptx(p1(i,q),aux))+ Alfa3(i,9*z-9+q)*(Ptx(p2(i,q),aux)-Ptx(p1(i,q),aux)));
            Ry3(i,9*z-9+q) = 1/2*((Pty(p2(i,q),aux)+Pty(p1(i,q),aux))+ Alfa3(i,9*z-9+q)*(Pty(p2(i,q),aux)-Pty(p1(i,q),aux))); 
            
            if (abs(Alfa3(i,9*z-9+q)) >= 1)
                Rx3(i,9*z-9+q) = NaN;
                Ry3(i,9*z-9+q) = NaN;
            end
            if (abs(Beta3(i,9*z-9+q)) >= 1 )
                Pontox3(i,9*z-9+q) = NaN;
                Pontoy3(i,9*z-9+q) = NaN;
            end
hold on
%plot(Rx3(i,9*z-9+q),Ry3(i,9*z-9+q),'.')
%plot(Pontox3(i,9*z-9+q),Pontoy3(i,9*z-9+q),'.');
% plot(Rx31(1,1),Ry31(1,1),'.');
% plot(Rx32(:,:),Ry32(:,:),'.');
 hold off            
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% do ponto de reflexão da terceira imagem de terceira ordem para ponto de segunda reflexão %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p01(i,q) = G1(i,q);
            p11(i,q) = p01(i,q) + 1;
            if p11(i,q) == 5
                 p11(i,q)  = 1;
            end
  
            k32(i,q) = (Pty(p11(i,q),aux)-Pty(p01(i,q),aux))/(Ptx(p11(i,q),aux)-Ptx(p01(i,q),aux));
            k32(~isfinite(k32))= 0;
            w32(i,9*z-9+q) = ((Ry3(i,9*z-9+q)-Iy(i,m))/(Rx3(i,9*z-9+q)-Ix(i,m)));
            Beta32(i,9*z-9+q) = ((-k32(i,q)*(Ix(i,m)+ Rx3(i,9*z-9+q)))+(Iy(i,m)+Ry3(i,9*z-9+q)) -(2*(Pty(p01(i,q),aux)-(k32(i,q)*Ptx(p01(i,q),aux)))))/((k32(i,q)*((Ix(i,m)-Rx3(i,9*z-9+q))))-((Iy(i,m)-Ry3(i,9*z-9+q))));
            Beta32(isnan(Beta32))= 0;
            Beta32(~isfinite(Beta32))= 0;
            Alfa32(i,9*z-9+q) = (-w32(i,9*z-9+q)*(Ptx(p11(i,q),aux)+Ptx(p01(i,q),aux))+(Pty(p11(i,q),aux)+Pty(p01(i,q),aux)) - (2*(Iy(i,m)-(w32(i,9*z-9+q)*Ix(i,m)))))/((w32(i,9*z-9+q)*(Ptx(p11(i,q),aux)-Ptx(p01(i,q),aux)))-(Pty(p11(i,q),aux)-Pty(p01(i,q),aux)));
            Pontox32(i,9*z-9+q) = 1/2*[(Ix(i,m)+ Rx3(i,9*z-9+q))+ Beta32(i,9*z-9+q)*(Rx3(i,9*z-9+q)-Ix(i,m))];
            Pontox32(~isfinite(Pontox32))= 0;
            Pontoy32(i,9*z-9+q) = 1/2*[(Iy(i,m)+Ry3(i,9*z-9+q))+ Beta32(i,9*z-9+q)*(Ry3(i,9*z-9+q)-Iy(i,m))];
            Rx32(i,9*z-9+q) = 1/2*((Ptx(p11(i,q),aux)+Ptx(p01(i,q),aux))+ Alfa32(i,9*z-9+q)*(Ptx(p11(i,q),aux)-Ptx(p01(i,q),aux)));
            Ry32(i,9*z-9+q) = 1/2*((Pty(p11(i,q),aux)+Pty(p01(i,q),aux))+ Alfa32(i,9*z-9+q)*(Pty(p11(i,q),aux)-Pty(p01(i,q),aux)));
            
            if (abs(Alfa32(i,9*z-9+q)) >= 1)
                Rx32(i,9*z-9+q) = NaN;
                Ry32(i,9*z-9+q) = NaN;
            end
            if (abs(Beta32(i,9*z-9+q)) >= 1)
                Pontox32 = NaN;
                Pontoy32 = NaN;
            end
            
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% do ponto de reflexão da segunda imagem de terceira ordem para ponto de primeira reflexão %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q0(i,q) = G2(i,q);
            p0(i,q) = q0(i,q) + 1;
             if p0(i,q) == 5
                 p0(i,q)  = 1;
            end
            
            k31(i,q) = (Pty(p0(i,q),aux)-Pty(q0(i,q),aux))/(Ptx(p0(i,q),aux)-Ptx(q0(i,q),aux));
            k31(~isfinite(k31))= 0;
            w31(i,9*z-9+q) = (Ry32(i,9*z-9+q)-Iy(i,aux))/(Rx32(i,9*z-9+q)-Ix(i,aux));
            Beta31(i,9*z-9+q) = ((-k31(i,q)*(Rx32(i,9*z-9+q)+Ix(i,aux)))+(Ry32(i,9*z-9+q)+Iy(i,aux)) -(2*(Pty(q0(i,q),aux)-(k31(i,q)*Ptx(q0(i,q),aux)))))/((k31(i,q)*((Rx32(i,9*z-9+q)-Ix(i,aux))))-((Ry32(i,9*z-9+q)-Iy(i,aux))));
            Beta31(isnan(Beta31))= 0;
            Beta31(~isfinite(Beta31))= 0;
            Alfa31(i,9*z-9+q) = (-w31(i,9*z-9+q)*(Ptx(p0(i,q),aux)+Ptx(q0(i,q),aux))+(Pty(p0(i,q),aux)+Pty(q0(i,q),aux)) - (2*(Iy(i,aux)-(w31(i,9*z-9+q)*Ix(i,aux)))))/((w31(i,9*z-9+q)*(Ptx(p0(i,q),aux)-Ptx(q0(i,q),aux)))-(Pty(p0(i,q),aux)-Pty(q0(i,q),aux)));
            Pontox31(i,9*z-9+q) = 1/2*[(Rx32(i,9*z-9+q)+Ix(i,aux))+ Beta31(i,9*z-9+q)*(Rx32(i,9*z-9+q)-Ix(i,aux))];
            Pontox31(~isfinite(Pontox31))= 0;
            Pontoy31(i,9*z-9+q) = 1/2*[(Ry32(i,9*z-9+q)+Iy(i,aux))+ Beta31(i,9*z-9+q)*(Ry32(i,9*z-9+q)-Iy(i,aux))];
            Rx31(i,9*z-9+q) = 1/2*((Ptx(p0(i,q),aux)+Ptx(q0(i,q),aux))+ Alfa31(i,9*z-9+q)*(Ptx(p0(i,q),aux)-Ptx(q0(i,q),aux)));
            Ry31(i,9*z-9+q) = 1/2*((Pty(p0(i,q),aux)+Pty(q0(i,q),aux))+ Alfa31(i,9*z-9+q)*(Pty(p0(i,q),aux)-Pty(q0(i,q),aux)));
            
            if (abs(Alfa31(i,9*z-9+q)) >= 1)
                Rx31(i,9*z-9+q) = NaN;
                Ry31(i,9*z-9+q) = NaN;
            end
            if (abs(Beta31(i,9*z-9+q)) >= 1)
                Pontox31(i,9*z-9+q) = NaN;
                Pontoy31(i,9*z-9+q) = NaN;
            end
            
                     end
            end                   
    end
 end
end
 hold on
% plot(Rx3(:,:),Ry3(:,:),'.')
%plot(Pontox3(:,:),Pontoy3(:,:),'.');
% plot(Rx31(1,1),Ry31(1,1),'.');
% plot(Rx32(:,:),Ry32(:,:),'.');
 hold off 

%% Para 4 reflexões 
if Nr > 3
  pt2 = [1 2 3 4];
  Refl = Ns*(Ns-1)*(Ns-1)*(Ns-1)
  for z = 1:N
    for i = 1:Ns 
            for m = 2:4
                     for q = 2*(m-2)+(m-1) : 2*(m-2)+ m + 1
                            N4 = find(pt2 ~= point(i,q));
                             for c = 2*(q-1) + q : 2*(q-1) + q + 2
                                 
                                 if (rem(c,3) == 0)
                              refx4(i,c) = N4(3);
                                    elseif (rem((c+1),3) == 0)
                               refx4(i,c) = N4(2);
                                    else
                              refx4(i,c) = N4(1);
                                 end
                                  
                              ptx0(i,c) = refx4(i,c);
                              ptx1(i,c) = (refx4(i,c)+1);
                             if ptx0(i,c) ==  4
                                 ptx1(i,c)  = 1;
                             end
             
        alfa4(i,c) = ((Ix3(i,q)-Ptx(ptx0(i,c),aux))*(Ptx(ptx1(i,c),aux)-Ptx(ptx0(i,c),aux))+(Iy3(i,q)-Pty(ptx0(i,c),aux))*(Pty(ptx1(i,c),aux)-Pty(ptx0(i,c),aux)))/((Ptx(ptx1(i,c),aux)-Ptx(ptx0(i,c),aux))^2+(Pty(ptx1(i,c),aux)-Pty(ptx0(i,c),aux))^2);  
        Qx4(i,c) = Ptx(ptx0(i,c),aux) + alfa4(i,c)*(Ptx(ptx1(i,c),aux)-Ptx(ptx0(i,c),aux));
        Qy4(i,c) = Pty(ptx0(i,c),aux) + alfa4(i,c)*(Pty(ptx1(i,c),aux)-Pty(ptx0(i,c),aux));

         if (abs(alfa4(i,c)) <= 1  && abs(alfa4(i,c)) >= 0)
           
             Ix4(i,c) = 2*Qx4(i,c) - Ix3(i,q);
             Iy4(i,c) = 2*Qy4(i,c) - Iy3(i,q);
        else
              Ix4(i,c) = NaN;
              Iy4(i,c) = NaN;
             
        end
        
        
            k4(i,c) = (Pty(ptx1(i,c),aux)-Pty(ptx0(i,c),aux))/(Ptx(ptx1(i,c),aux)-Ptx(ptx0(i,c),aux));
            k4(~isfinite(k4))= 0;
            w4(i,27*z-27+c) = ((Iy4(i,c)-yTx)/(Ix4(i,c)-xTx(z)));
            Beta4(i,27*z-27+c) = ((-k4(i,c)*(xTx(z)+Ix4(i,c)))+(yTx+Iy4(i,c)) -(2*(Pty(ptx0(i,c),aux)-(k4(i,c)*Ptx(ptx0(i,c),aux)))))/((k4(i,c)*((Ix4(i,c)-xTx(z))))-(Iy4(i,c)-yTx));
            Beta4(isnan(Beta4))= 0;
            Beta4(~isfinite(Beta4))= 0;
            Alfa4(i,27*z-27+c) = (-w4(i,27*z-27+c)*(Ptx(ptx1(i,c),aux)+Ptx(ptx0(i,c),aux))+(Pty(ptx1(i,c),aux)+Pty(ptx0(i,c),aux)) - (2*(yTx-(w4(i,27*z-27+c)*xTx(z)))))/((w4(i,27*z-27+c)*(Ptx(ptx1(i,c),aux)-Ptx(ptx0(i,c),aux)))-(Pty(ptx1(i,c),aux)-Pty(ptx0(i,c),aux)));
            Pontox4(i,27*z-27+c) = 1/2*[(xTx(z) + Ix4(i,c))+ Beta4(i,27*z-27+c)*(Ix4(i,c)-xTx(z))];
            Pontox4(~isfinite(Pontox4))= 0;
            Pontoy4(i,27*z-27+c) = 1/2*[(yTx+Iy4(i,c))+ Beta4(i,27*z-27+c)*(Iy4(i,c)-yTx)];
            Rx4(i,27*z-27+c) = 1/2*((Ptx(ptx1(i,c),aux)+Ptx(ptx0(i,c),aux))+ Alfa4(i,27*z-27+c)*(Ptx(ptx1(i,c),aux)-Ptx(ptx0(i,c),aux)));
            Ry4(i,27*z-27+c) = 1/2*((Pty(ptx1(i,c),aux)+Pty(ptx0(i,c),aux))+ Alfa4(i,27*z-27+c)*(Pty(ptx1(i,c),aux)-Pty(ptx0(i,c),aux))); 
            
            if (abs(Alfa4(i,27*z-27+c)) >= 1)
                Rx4(i,27*z-27+c) = NaN;
                Ry4(i,27*z-27+c) = NaN;
            end
            
            
% hold on
%         plot(Ix4(2,1),Iy4(2,1),'.');
%          plot(Rx4(i,27*z-27+c),Ry4(i,27*z-27+c),'.');
%        plot(Rx43(:,:),Ry43(:,:),'.');
%  hold off 
            
            limt01(i,c) = point(i,q);
            limt0(i,c) = (limt01(i,c)+1);
            if limt01(i,c) == 4
                 limt0(i,c)  = 1;
            end  
            
            aux = 1;
            k43(i,c) = (Pty(limt0(i,c),aux)-Pty(limt01(i,c),aux))/(Ptx(limt0(i,c),aux)-Ptx(limt01(i,c),aux));
            k43(~isfinite(k43))= 0;
            w43(i,27*z-27+c) = ((Ry4(i,27*z-27+c)-Iy3(i,q))/(Rx4(i,27*z-27+c)-Ix3(i,q)));
            Beta43(i,27*z-27+c) = ((-k43(i,c)*(Ix3(i,q)+ Rx4(i,27*z-27+c)))+(Ry4(i,27*z-27+c)+Iy3(i,q)) -(2*(Pty(limt01(i,c),aux)-(k43(i,c)*Ptx(limt01(i,c),aux)))))/((k43(i,c)*(Rx4(i,27*z-27+c)-Ix3(i,q)))-(Ry4(i,27*z-27+c)-Iy3(i,q)));
            Beta43(isnan(Beta43))= 0;
            Beta43(~isfinite(Beta43))= 0;
            Alfa43(i,27*z-27+c) = ((-w43(i,27*z-27+c)*(Ptx(limt0(i,c),aux)+Ptx(limt01(i,c),aux)))+(Pty(limt0(i,c),aux)+Pty(limt01(i,c),aux)) - (2*(Iy3(i,q)-(w43(i,27*z-27+c)*Ix3(i,q)))))/((w43(i,27*z-27+c)*(Ptx(limt0(i,c),aux)-Ptx(limt01(i,c),aux)))-(Pty(limt0(i,c),aux)-Pty(limt01(i,c),aux)));
            Pontox43(i,27*z-27+c) = 1/2*[(Ix3(i,q)+ Rx4(i,27*z-27+c))+ Beta43(i,27*z-27+c)*(Ix3(i,q)-Rx4(i,27*z-27+c))];
            Pontox43(~isfinite(Pontox43))= 0;
            Pontoy43(i,27*z-27+c) = 1/2*[(Ry4(i,27*z-27+c)+Iy3(i,q))+ Beta43(i,27*z-27+c)*(Iy3(i,q)-Ry4(i,27*z-27+c))];
            Rx43(i,27*z-27+c) = 1/2*((Ptx(limt0(i,c),aux)+Ptx(limt01(i,c),aux))+ Alfa43(i,27*z-27+c)*(Ptx(limt0(i,c),aux)-Ptx(limt01(i,c),aux)));
            Ry43(i,27*z-27+c) = 1/2*((Pty(limt0(i,c),aux)+Pty(limt01(i,c),aux))+ Alfa43(i,27*z-27+c)*(Pty(limt0(i,c),aux)-Pty(limt01(i,c),aux))); 
            
            if (abs(Alfa43(i,27*z-27+c)) >= 1)
                Rx43(i,27*z-27+c) = NaN;
                Ry43(i,27*z-27+c) = NaN;
            end
 
hold on
%         plot(Ix4(2,1),Iy4(2,1),'.');
%         plot(Rx43(i,27*z-27+c),Ry43(i,27*z-27+c),'.');
%        plot(Rx43(:,:),Ry43(:,:),'.');
 hold off
 
            ptx2(i,c) = G1(i,q);
            ptx3(i,c) = (ptx2(i,c)+1);
               if ptx2(i,c) == 4
                 ptx3(i,c)  = 1;
               end  

            k44(i,c) = (Pty(ptx3(i,c),aux)-Pty(ptx2(i,c),aux))/(Ptx(ptx3(i,c),aux)-Ptx(ptx2(i,c),aux));
            k44(~isfinite(k44))= 0;
            w44(i,27*z-27+c) = ((Ry43(i,27*z-27+c)-Iy(i,m))/(Rx43(i,27*z-27+c)-Ix(i,m)));
            Beta44(i,27*z-27+c) = ((-k44(i,c)*(Ix(i,m)+ Rx43(i,27*z-27+c)))+(Ry43(i,27*z-27+c)+Iy(i,m)) -(2*(Pty(ptx2(i,c),aux)-(k44(i,c)*Ptx(ptx2(i,c),aux)))))/((k44(i,c)*(Rx43(i,27*z-27+c)-Ix(i,m)))-(Ry43(i,27*z-27+c)-Iy(i,m)));
            Beta44(isnan(Beta44))= 0;
            Beta44(~isfinite(Beta44))= 0;
            Alfa44(i,27*z-27+c) = ((-w44(i,27*z-27+c)*(Ptx(ptx3(i,c),aux)+Ptx(ptx2(i,c),aux)))+(Pty(ptx3(i,c),aux)+Pty(ptx2(i,c),aux)) - (2*(Iy(i,m)-(w44(i,27*z-27+c)*Ix(i,m)))))/((w44(i,27*z-27+c)*(Ptx(ptx3(i,c),aux)-Ptx(ptx2(i,c),aux)))-(Pty(ptx3(i,c),aux)-Pty(ptx2(i,c),aux)));
            Pontox44(i,27*z-27+c) = 1/2*[(Ix(i,m)+ Rx43(i,27*z-27+c))+ Beta44(i,27*z-27+c)*(Ix(i,m)-Rx43(i,27*z-27+c))];
            Pontox44(~isfinite(Pontox44))= 0;
            Pontoy44(i,27*z-27+c) = 1/2*[(Ry43(i,27*z-27+c)+Iy(i,m))+ Beta44(i,27*z-27+c)*(Iy(i,m)-Ry43(i,27*z-27+c))];
            Rx44(i,27*z-27+c) = 1/2*((Ptx(ptx3(i,c),aux)+Ptx(ptx2(i,c),aux))+ Alfa44(i,27*z-27+c)*(Ptx(ptx3(i,c),aux)-Ptx(ptx2(i,c),aux)));
            Ry44(i,27*z-27+c) = 1/2*((Pty(ptx3(i,c),aux)+Pty(ptx2(i,c),aux))+ Alfa44(i,27*z-27+c)*(Pty(ptx3(i,c),aux)-Pty(ptx2(i,c),aux))); 
            
            if (abs(Alfa44(i,27*z-27+c)) >= 1)
                Rx44(i,27*z-27+c) = NaN;
                Ry44(i,27*z-27+c) = NaN;
            end
hold on
%         plot(Ix4(2,1),Iy4(2,1),'.');
%         plot(Rx44(i,27*z-27+c),Ry44(i,27*z-27+c),'.');
%        plot(Rx43(:,:),Ry43(:,:),'.');
 hold off
 
 
            ptx4(i,c) = G2(i,q);
            ptx5(i,c) = (ptx4(i,c)+1);
               if ptx4(i,c) == 4
                 ptx5(i,c)  = 1;
               end  

            k40(i,c) = (Pty(ptx5(i,c),aux)-Pty(ptx4(i,c),aux))/(Ptx(ptx5(i,c),aux)-Ptx(ptx4(i,c),aux));
            k40(~isfinite(k40))= 0;
            w40(i,27*z-27+c) = (( Ry44(i,27*z-27+c)-Iy(i,1))/( Rx44(i,27*z-27+c)-Ix(i,1)));
            Beta40(i,27*z-27+c) = ((-k40(i,c)*(Ix(i,1)+  Rx44(i,27*z-27+c)))+( Ry44(i,27*z-27+c)+Iy(i,1)) -(2*(Pty(ptx4(i,c),aux)-(k40(i,c)*Ptx(ptx4(i,c),aux)))))/((k40(i,c)*( Rx44(i,27*z-27+c)-Ix(i,1)))-( Ry44(i,27*z-27+c)-Iy(i,1)));
            Beta40(isnan(Beta40))= 0;
            Beta40(~isfinite(Beta40))= 0;
            Alfa40(i,27*z-27+c) = ((-w40(i,27*z-27+c)*(Ptx(ptx5(i,c),aux)+Ptx(ptx4(i,c),aux)))+(Pty(ptx5(i,c),aux)+Pty(ptx4(i,c),aux)) - (2*(Iy(i,1)-(w40(i,27*z-27+c)*Ix(i,1)))))/((w40(i,27*z-27+c)*(Ptx(ptx5(i,c),aux)-Ptx(ptx4(i,c),aux)))-(Pty(ptx5(i,c),aux)-Pty(ptx4(i,c),aux)));
            Pontox40(i,27*z-27+c) = 1/2*[(Ix(i,1)+  Rx44(i,27*z-27+c))+ Beta40(i,27*z-27+c)*(Ix(i,1)- Rx44(i,27*z-27+c))];
            Pontox40(~isfinite(Pontox40))= 0;
            Pontoy40(i,27*z-27+c) = 1/2*[( Ry44(i,27*z-27+c)+Iy(i,1))+ Beta40(i,27*z-27+c)*(Iy(i,1)- Ry44(i,27*z-27+c))];
            Rx40(i,27*z-27+c) = 1/2*((Ptx(ptx5(i,c),aux)+Ptx(ptx4(i,c),aux))+ Alfa40(i,27*z-27+c)*(Ptx(ptx5(i,c),aux)-Ptx(ptx4(i,c),aux)));
            Ry40(i,27*z-27+c) = 1/2*((Pty(ptx5(i,c),aux)+Pty(ptx4(i,c),aux))+ Alfa40(i,27*z-27+c)*(Pty(ptx5(i,c),aux)-Pty(ptx4(i,c),aux))); 
            
            if (abs(Alfa40(i,27*z-27+c)) >= 1)
                Rx40(i,27*z-27+c) = NaN;
                Ry40(i,27*z-27+c) = NaN;
            end
            
% hold on
% %         plot(Ix4(2,1),Iy4(2,1),'.');
%          plot(Rx40(i,27*z-27+c),Ry40(i,27*z-27+c),'.');
% %        plot(Rx43(:,:),Ry43(:,:),'.');
%  hold off
 
 
 
          end
        end
    end
  end
  end
end
 hold on
%         plot(Ix4(2,1),Iy4(2,1),'.');
%         plot(Rx4(:,:),Ry4(:,:),'.');
%        plot(Rx43(:,:),Ry43(:,:),'.');
 hold off 
                         

%% PLOTAGEM DOS PERCURSSOS ENVOLVIDOS
hold on
Tst = 100;
line ([xTi,xTx(Tst)],[yTi,yTx],'Color','m');
 
for i = 1:Ns
    line ([xTi,Sx(i,Tst),xTx(Tst)],[yTi,Sy(i,Tst),yTx]);
    
    if Nr > 1
        for m = 2:Elementos
            if (abs(Alfa(i,4*Tst-4+m))< 1   && abs(Alfa1(i,4*Tst-4+m))< 1  )
            line ([xTi,Rx1(i,4*Tst-4+m),Rx(i,4*Tst-4+m),xTx(Tst)],[yTi,Ry1(i,4*Tst-4+m),Ry(i,4*Tst-4+m),yTx]);
            end
                    if Nr > 2
                     for q = 2*(m-2)+(m-1) : 2*(m-2)+ m + 1 
                             if (abs(Alfa31(i,9*Tst-9+q))< 1  && abs(Alfa32(i,9*Tst-9+q)) < 1  && abs(Alfa3(i,9*Tst-9+q))< 1 )
                             line ([xTi,Rx31(i,9*Tst-9+q),Rx32(i,9*Tst-9+q),Rx3(i,9*Tst-9+q),xTx(Tst)],[yTi,Ry31(i,9*Tst-9+q),Ry32(i,9*Tst-9+q),Ry3(i,9*Tst-9+q),yTx]);
                             end
                             if Nr > 3
                                        for c = 2*(q-1) + q : 2*(q-1) + q + 2
                                            if (abs(Alfa4(i,27*Tst-27+c))< 1  && abs(Alfa43(i,27*Tst-27+c)) < 1  && abs(Alfa40(i,27*Tst-27+c))< 1 && abs(Alfa44(i,27*Tst-27+c))< 1 )
                                           line ([xTi,Rx40(i,27*Tst-27+c),Rx44(i,27*Tst-27+c),Rx43(i,27*Tst-27+c),Rx4(i,27*Tst-27+c),xTx(Tst)],[yTi,Ry40(i,27*Tst-27+c),Ry44(i,27*Tst-27+c),Ry43(i,27*Tst-27+c),Ry4(i,27*Tst-27+c),yTx])
                                            end
                                        end
                             end
                     end
               end
        end
    end
end

hold off

%% -------------------------------------------- CÁLCULO DE CAMPOS -----------------------------------------
%Dados de entrada
rho = 0.15;
Pt = 1600;
Gt = 1;
E0 = sqrt((N0*Pt*Gt)/(2*pi))*1; % Ponto de Excitação
%E0 = 180
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Raio direto
 for z = 1:N
    im = sqrt(-1);
    dist(z) = distance(xTi,yTi,xTx(z),yTi);
    A(z) = sqrt(rho/dist(z)); %Fator de espalhamento
    
    if dist(z)<= rho 
        Ein(z) = 0;
    else   
    Ein(z) = E0*A(z)*(exp(-im*B*dist(z))); %Campo elétrico referente na antena receptora
    Ein(~isfinite(Ein))= 0;

    end
 end
 
% Campo para primeira reflexão

 for z = 1:N
  for i = 1:Ns
      
     im = sqrt(-1);
     dist1(i,z) = distance(xTi,yTi,Sx(i,z),Sy(i,z)); %s
     %dist1(i,z) = distance(Ix(i,1),Iy(i,1),Sx(i,z),Sy(i,z)); %s
     %distI(i,z) = distance(Sx(i,z),Sy(i,z),xTx(z),yTx); %s'
     Esp(i,z) = sqrt(rho/dist1(i,z)) %Fator de espalhamento
    
     
     Ein0(i,z) = E0*Esp(i,z)*(exp(-im*B*dist1(i,z)));
     
%      %%%%%%%%%%%%%%%%%%%angulos
     
     m1(i,z) = (Sy(i,z) - yTi)/(Sx(i,z) - xTi);
     m2(i,z) = (Sy(i,z)-yTx)/(Sx(i,z) - xTx(z));
     m2(isnan(m2))= 0;
     m2(~isfinite(m2))= 0;
     an0(i,z)= abs((m2(i,z)-m1(i,z))/(1 + m1(i,z)*m2(i,z)));
     an0(isnan(an0))= 0;
     an1(i,z)= atand(an0(i,z));
     tetai(i,z) = an1(i,z)/2;
     an2(i,z) = cosd(tetai(i,z));
     
     cosTheta(i,z) = (sqrt(1-(((u0*e0)/(uc*ec))^2)*sind(tetai(i,z))*sind(tetai(i,z))));
     an2(isnan(an2))= 0;

% %      %%%%%%%%%%%%%%%%%%%%
     
%      Nincn(i,z) = (Ninc/cosTheta(i,z))*[[(N0/an2(i,z))+(Ninc/cosTheta(i,z))*tanh(Yc*d*cosTheta(i,z))]/((Ninc/cosTheta(i,z))+(N0/an2(i,z))*tanh(Yc*d*cosTheta(i,z)))]; %Componente normal da impedância de entrada
%      Nincp(i,z) = (Ninc*cosTheta(i,z))*[[(N0*an2(i,z))+(Ninc*cosTheta(i,z))*tanh(Yc*d*cosTheta(i,z))]/[(Ninc*cosTheta(i,z))+(N0*an2(i,z))*tanh(Yc*d*cosTheta(i,z))]]; %Componente paralela da impedância de entrada
%      coefref1(i,z) =  (Nincn(i,z) - (N0*an2(i,z)))/(Nincp(i,z) + (N0*an2(i,z)));
     
     coefref1(i,z) = (((Ninc*cosTheta(i,z)- N0*an2(i,z))/(Ninc*cosTheta(i,z) + N0*an2(i,z))));
%      coefrefl(i,z) = ((N0*cosTheta(i,z)- Ninc*an2(i,z))/(N0*cosTheta(i,z) + Ninc*an2(i,z)));
     coefref10(i,z)=  abs(coefref1(i,z));
     
     Eir(i,z) = coefref10(i,z)*(Ein0(i,z));
     Eir(isnan(Eir))= 0;
     dist2(i,z) = distance(Sx(i,z),Sy(i,z),xTx(z),yTx);
     pIm(i,z) = distance(Ix(i,1),Iy(i,1),Sx(i,z),Sy(i,z));
     Esp1(i,z) = sqrt(pIm(i,z)/(dist2(i,z)+ pIm(i,z))); %Fator de espalhamento
     Esp1(~isfinite(Esp1))= 0;
     distTx(z) = distance(xTi,yTi,xTx(z),yTi);
     
     Ein1(i,z) = (Eir(i,z))*Esp1(i,z)*exp(-im*B*dist2(i,z));
     Ein1(~isfinite(Ein1))= 0;
     
     
     
      if (distTx(z)<= rho || abs(ALFA(i,z)) > 1 )
          
        Ein1(i,z) = 0;
       % Ein0(i,z) = 0;
        
      end
  end
 end
 
%  CAMPO PARA 2 REFLEXÕES     
for z = 1:N
    for i = 1:Ns   
        for m = 2:Elementos                    
               
                     
                            im = sqrt(-1);
                            distR(i,4*z-4+m) = distance(xTi,yTi,Rx1(i,4*z-4+m),Ry1(i,4*z-4+m)); %s
                           % EspR(i,4*z-4+m) = sqrt(distR(i,4*z-4+m)/(distRI(i,4*z-4+m)+distR(i,4*z-4+m))); %Fator de espalhamento
                            EspR(i,4*z-4+m) = sqrt(rho/distR(i,4*z-4+m));
                            Einput(i,4*z-4+m) = E0*EspR(i,4*z-4+m)*exp(-im*B*distR(i,4*z-4+m));
                            
                     %%%%%%%%%%%%%%%%%%%%%%%%%angulos
                            m11(i,4*z-4+m) = (Ry1(i,4*z-4+m) - yTi)/(Rx1(i,4*z-4+m) - xTi);
                            m3(i,4*z-4+m) = (Ry(i,4*z-4+m) - Ry1(i,4*z-4+m))/(Rx(i,4*z-4+m) - Rx1(i,4*z-4+m));
                            ANG0(i,4*z-4+m)= abs((m3(i,4*z-4+m)-  m11(i,4*z-4+m))/(1 + m11(i,4*z-4+m)*m3(i,4*z-4+m)));
                            ANG0(isnan(ANG0))= 0;
                            TETA0(i,4*z-4+m)= atand(ANG0(i,4*z-4+m));
                            aux2(i,4*z-4+m) = TETA0(i,4*z-4+m)/2;
                            aux2(isnan(aux2))= 0;
                            cosTheta1(i,4*z-4+m) = (sqrt(1-((((u0*e0)/(uc*ec))^2)*sind(aux2(i,4*z-4+m))*sind(aux2(i,4*z-4+m)))));
                            angulos1(i,4*z-4+m) = cosd(aux2(i,4*z-4+m));
                            angulos1(isnan(angulos1))= 0;
                             
                     %%%%%%%%%%%%%%%%%%%%%%%%%%
                     
%                           Nincn1(i,4*z-4+m) = (Ninc/cosTheta1(i,4*z-4+m))*[(N0/angulos1(i,4*z-4+m))+(Ninc/cosTheta1(i,4*z-4+m))*tanh(Yc*d*cosTheta(i,z))]/((Ninc/cosTheta1(i,4*z-4+m))+(N0/angulos1(i,4*z-4+m))*tanh(Yc*d*cosTheta1(i,4*z-4+m))); %Componente normal da impedância de entrada
%                           Nincp1(i,4*z-4+m) = (Ninc*cosTheta1(i,4*z-4+m))*[(N0*angulos1(i,4*z-4+m))+(Ninc*cosTheta1(i,4*z-4+m))*tanh(Yc*d*cosTheta(i,z))]/[(Ninc*cosTheta1(i,4*z-4+m))+(N0*angulos1(i,4*z-4+m))*tanh(Yc*d*cosTheta1(i,4*z-4+m))]; %Componente paralela da impedância de entrada
%                           coefref20(i,4*z-4+m) =  (Nincn1(i,4*z-4+m) - (N0*angulos1(i,4*z-4+m)))/(Nincp1(i,4*z-4+m) + (N0*angulos1(i,4*z-4+m)));

                     coefref20(i,4*z-4+m) = ([Ninc*cosTheta1(i,4*z-4+m)- N0*angulos1(i,4*z-4+m)]/[Ninc*cosTheta1(i,4*z-4+m) + N0*angulos1(i,4*z-4+m)]);
                     coefref2(i,4*z-4+m)=  abs(coefref20(i,4*z-4+m));
                     
                     Eir1(i,4*z-4+m) = coefref2(i,4*z-4+m)*(Einput(i,4*z-4+m));
                     pIm1(i,4*z-4+m) = distance(Ix(i,1),Iy(i,1),Rx1(i,4*z-4+m),Ry1(i,4*z-4+m));
                     dist3(i,4*z-4+m) = distance(Rx1(i,4*z-4+m),Ry1(i,4*z-4+m),Rx(i,4*z-4+m),Ry(i,4*z-4+m));
                     Esp2(i,4*z-4+m) = sqrt(pIm1(i,4*z-4+m)/(dist3(i,4*z-4+m)+ pIm1(i,4*z-4+m))); %Fator de espalhamento
                     Eir2(i,4*z-4+m) = Eir1(i,4*z-4+m)*Esp2(i,4*z-4+m)*exp(-im*B*dist3(i,4*z-4+m));
                      
                          
                    

                          %%%%%%%%%%%%%%%%%%%angulos
                            m4(i,4*z-4+m) = (yTx-Ry(i,4*z-4+m))/(xTx(z)-Rx(i,4*z-4+m));
                            m4(~isfinite(m4))= 0;
                            m4(isnan(m4))= 0;
                            ANG1(i,4*z-4+m)= abs((m3(i,4*z-4+m)-m4(i,4*z-4+m))/(1 + m4(i,4*z-4+m)*m3(i,4*z-4+m)));
                            ANG1(isnan(ANG1))= 0;
                            TETA1(i,4*z-4+m)= atand(ANG1(i,4*z-4+m));
                            TETA1(isnan(TETA1))= 0;
                            aux3(i,4*z-4+m) = TETA1(i,4*z-4+m)/2;
                     
                            cosTheta2(i,4*z-4+m) = (sqrt(1-((((u0*e0)/(uc*ec))^2)*sind(aux3(i,4*z-4+m))*sind(aux3(i,4*z-4+m)))));
                            angulos3(i,4*z-4+m) = cosd(aux3(i,4*z-4+m));
                     %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                     
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     
                     
                     Nincn2(i,4*z-4+m) = (Ninc/cosTheta2(i,4*z-4+m))*[(N0/angulos3(i,4*z-4+m))+(Ninc/cosTheta2(i,4*z-4+m))*tanh(Yc*d*cosTheta2(i,4*z-4+m))]/((Ninc/cosTheta2(i,4*z-4+m))+(N0/angulos3(i,4*z-4+m))*tanh(Yc*d*cosTheta2(i,4*z-4+m))); %Componente normal da impedância de entrada
                     Nincp2(i,4*z-4+m) = (Ninc*cosTheta2(i,4*z-4+m))*[(N0*angulos3(i,4*z-4+m))+(Ninc*cosTheta2(i,4*z-4+m))*tanh(Yc*d*cosTheta2(i,4*z-4+m))]/[(Ninc*cosTheta2(i,4*z-4+m))+(N0*angulos3(i,4*z-4+m))*tanh(Yc*d*cosTheta2(i,4*z-4+m))]; %Componente paralela da impedância de entrada
                     coefref30(i,4*z-4+m) =  (Nincn2(i,4*z-4+m) - (N0*angulos3(i,4*z-4+m)))/(Nincp2(i,4*z-4+m) + (N0*angulos3(i,4*z-4+m)));
                     
%                      coefref30(i,4*z-4+m) = ([Ninc*cosTheta2(i,4*z-4+m) - N0*angulos3(i,4*z-4+m)]/[Ninc*cosTheta2(i,4*z-4+m) + N0*angulos3(i,4*z-4+m)]);
                      coefref3(i,4*z-4+m)=  abs(coefref30(i,4*z-4+m));
                     
                     Eir3(i,4*z-4+m) = coefref3(i,4*z-4+m)*(Eir2(i,4*z-4+m));  
                     dist4(i,4*z-4+m) = distance(Rx(i,4*z-4+m),Ry(i,4*z-4+m),xTx(z),yTx);
                     pIm2(i,4*z-4+m) = distance(Ix(i,m),Iy(i,m),Rx(i,4*z-4+m),Ry(i,4*z-4+m));
                     Esp3(i,4*z-4+m) = sqrt(pIm2(i,4*z-4+m)/(dist4(i,4*z-4+m)+pIm2(i,4*z-4+m)));
                     Esp3(~isfinite(Esp3))= 0;
                     Efinal(i,4*z-4+m) = Eir3(i,4*z-4+m)*Esp3(i,4*z-4+m)*exp(-im*B*dist4(i,4*z-4+m));
                     Efinal(~isfinite(Efinal))= 0;
                     distTx2(z) = distance(xTi,yTi,xTx(z),yTi);
                     
                     if ( distTx2(z)<= rho || abs(Alfa(i,4*z-4+m)) > 1  || abs(Alfa1(i,4*z-4+m)) > 1 )
                            
                         Efinal(i,4*z-4+m) = 0;
                         
                      end
                    
                      
         end
    end
end
 
% Campo para 3 reflexão
% for z = 1:N
%     for i = 1:Ns   
%         for m = 2:Elementos
%              for q = 2*(m-2)+(m-1) : 2*(m-2)+ m + 1
%                              im = sqrt(-1);
%                              dist3R(i,9*z-9+q) = distance(Rx31(i,9*z-9+q),Ry31(i,9*z-9+q),xTi,yTi); %s
%                            % EspR(i,4*z-4+m) = sqrt(distR(i,4*z-4+m)/(distRI(i,4*z-4+m)+distR(i,4*z-4+m))); %Fator de espalhamento
%                             Esp3R(i,9*z-9+q) = sqrt(rho/dist3R(i,4*z-4+m));
%                             Einput3R(i,9*z-9+q) = E0*Esp3R(i,9*z-9+q)*exp(-im*B*dist3R(i,9*z-9+q));
%                             
%                            %%%%%%%%%%%%%%%%%%%%%%%%%angulos
%                             m31(i,9*z-9+q) = (Ry31(i,9*z-9+q) - yTi)/(Rx31(i,9*z-9+q) - xTi);
%                             m32(i,9*z-9+q) = (Ry31(i,9*z-9+q)- Ry32(i,9*z-9+q))/(Rx31(i,9*z-9+q) - Rx32(i,9*z-9+q));
%                             ANG3R(i,9*z-9+q)= abs((m31(i,9*z-9+q) -  m32(i,9*z-9+q))/(1 + m31(i,9*z-9+q)*m32(i,9*z-9+q)));
%                             ANG3R(isnan(ANG3R))= 0;
%                             TETA3R(i,9*z-9+q)= atand(ANG3R(i,9*z-9+q));
%                             aux3R(i,9*z-9+q) = TETA3R(i,9*z-9+q)/2;
%                             aux3R(isnan(aux2))= 0;
%                             cosTheta3R(i,9*z-9+q) = (sqrt(1-((((u0*e0)/(uc*ec))^2)*sind(aux3R(i,9*z-9+q))*sind(aux3R(i,9*z-9+q)))));
%                             angulos3R(i,9*z-9+q) = cosd(aux3R(i,9*z-9+q));
%                             angulos3R(isnan(angulos3R))= 0;
%                              
%                             %%%%%%%%%%%%%%%%%%%%%%%%%%
%                             coefref3R(i,9*z-9+q) = ([Ninc*cosTheta3R(i,9*z-9+q)- N0*angulos3R(i,9*z-9+q)]/[Ninc*cosTheta3R(i,9*z-9+q) + N0*angulos3R(i,9*z-9+q)]);
%                             coefref30R(i,9*z-9+q)=  abs(coefref3R(i,9*z-9+q));
%                             
%                             Eir3R(i,9*z-9+q) = coefref3R(i,9*z-9+q)*(Einput3R(i,9*z-9+q));
%                             pIm3(i,9*z-9+q) = distance(Ix(i,1),Iy(i,1),Rx31(i,9*z-9+q),Ry31(i,9*z-9+q));
%                             dist3R(i,9*z-9+q) = distance(Rx31(i,9*z-9+q),Ry31(i,9*z-9+q),Rx32(i,9*z-9+q),Ry32(i,9*z-9+q));
%                             Esp2(i,9*z-9+q) = sqrt(pIm3(i,9*z-9+q)/(dist3R(i,9*z-9+q)+ pIm3(i,9*z-9+q))); %Fator de espalhamento
%                             Eir31(i,9*z-9+q) = Eir3R(i,9*z-9+q)*Esp3R(i,9*z-9+q)*exp(-im*B*dist3(i,9*z-9+q));
%                       
%                           
%                     
% 
%                           %%%%%%%%%%%%%%%%%%%angulos
%                             m33(i,9*z-9+q) = (Ry3(i,9*z-9+q)-Ry32(i,9*z-9+q))/(Rx3(i,9*z-9+q)-Rx32(i,9*z-9+q));
%                             m33(~isfinite(m33))= 0;
%                             m33(isnan(m33))= 0;
%                             ANG31R(i,9*z-9+q)= abs((m32(i,9*z-9+q)-m33(i,i,9*z-9+q))/(1 + m32(i,9*z-9+q)*m33(i,9*z-9+q)));
%                             ANG31R(isnan(ANG1))= 0;
%                             TETA31(i,9*z-9+q)= atand(ANG31R(i,9*z-9+q));
%                             TETA31(isnan(TETA31))= 0;
%                             aux31(i,9*z-9+q) = TETA31(i,9*z-9+q)/2;
%                      
%                             cosTheta31R(i,9*z-9+q) = (sqrt(1-((((u0*e0)/(uc*ec))^2)*sind(aux31(i,9*z-9+q))*sind(aux31(i,9*z-9+q)))));
%                             angulos31(i,9*z-9+q) = cosd(aux31(i,9*z-9+q));
%                      %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                      
%                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      
% %                      
% %                      Nincn2(i,9*z-9+q) = (Ninc/cosTheta2(i,(4*z-4+m)))*[(N0/angulos3(i,(4*z-4+m)))+(Ninc/cosTheta2(i,(4*z-4+m)))*tanh(Yc*d*cosTheta2(i,(4*z-4+m)))]/((Ninc/cosTheta2(i,(4*z-4+m)))+(N0/angulos3(i,(4*z-4+m)))*tanh(Yc*d*cosTheta2(i,(4*z-4+m)))); %Componente normal da impedância de entrada
% %                      Nincp2(i,9*z-9+q) = (Ninc*cosTheta2(i,(4*z-4+m)))*[(N0*angulos3(i,(4*z-4+m)))+(Ninc*cosTheta2(i,(4*z-4+m)))*tanh(Yc*d*cosTheta2(i,(4*z-4+m)))]/[(Ninc*cosTheta2(i,(4*z-4+m)))+(N0*angulos3(i,(4*z-4+m)))*tanh(Yc*d*cosTheta2(i,(4*z-4+m)))]; %Componente paralela da impedância de entrada
% %                      coefref30(i,9*z-9+q) =  (Nincn2(i,9*z-9+q) - (N0*angulos3(i,9*z-9+q)))/(Nincp2(i,9*z-9+q) + (N0*angulos3(i,9*z-9+q)));
%                      
%                      coefref31R(i,9*z-9+q) = ([Ninc*cosTheta31R(i,9*z-9+q) - N0*angulos31(i,9*z-9+q)]/[Ninc*cosTheta31R(i,9*z-9+q) + N0*angulos31(i,9*z-9+q)]);
%                      coefref31(i,9*z-9+q)=  abs(coefref31R(i,9*z-9+q));
%                      
%                      Eir32(i,9*z-9+q) = coefref31(i,9*z-9+q)*(Eir31(i,9*z-9+q));  
%                      dist31(i,9*z-9+q) = distance(Rx32(i,9*z-9+q),Ry32(i,9*z-9+q),Rx3(i,9*z-9+q),Ry3(i,9*z-9+q));
%                      pIm31(i,9*z-9+q) = distance(Ix(i,m),Iy(i,m),Rx(i,9*z-9+q),Ry(i,9*z-9+q));
%                      Esp31R(i,9*z-9+q) = sqrt(pIm31(i,9*z-9+q)/(dist31(i,9*z-9+q)+pIm31(i,9*z-9+q)));
%                      Esp31R(~isfinite(Esp31R))= 0;
%                      Efinal3R(i,9*z-9+q) = Eir32(i,9*z-9+q)*Esp31R(i,9*z-9+q)*exp(-im*B*dist31(i,9*z-9+q));
%                      Efinal3R(~isfinite(Efinal3R))= 0;
%                      distTx3(z) = distance(xTi,yTi,xTx(z),yTi);
%                      
%                      if ( distTx3(z)<= rho || abs(Alfa3(i,9*z-9+q)) > 1  || abs(Alfa31(i,9*z-9+q)) > 1 )
%                             
%                          Efinal3R(i,9*z-9+q) = 0;
%              end
%         end
%         end
%     end
% end
    
 S =  Efinal;
 CAMPO = sum(S,1);
 CAMPO(~isfinite(CAMPO))= 0;
 
 for z = 1:N
    for m = 2 : Elementos
    EE(m,z) = CAMPO(1,4*z-4+m);
    end
 end

 EFINAL = EE;
 E2R = sum(EFINAL,1); 
 

D = Ein1;
D(~isfinite(D))= 0;
Ein2 = sum(D,1);

%  for z = 1:N
%     for m = 2 : Elementos
%         for q = (m-2)+(m-1) : (m-2)+ m
%             Campo
                        
    for z = 1:N
        Etotal(z) = 10*log((Ein(z) + Ein2(z) + E2R(z))/N0);
        Etotal(~isfinite(Etotal))= 0;
       

    end
    
 figure
plot(dist, Etotal);
 xlim([0 1.5]);
grid on;
time=toc