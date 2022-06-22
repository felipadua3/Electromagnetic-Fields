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

%%%% célula FDTD 
fx = 3;
fy = 1.75;

%Ponto S2
fx1 = 3;
fy1 = 3.25;

%Ponto S3
fx2 = 4.5;
fy2 = 3.25;

%Ponto S4
fx3 = 4.5;
fy3 = 1.75;

%%Criação de segmentos da área de cálculo
fdtd1 = line([fx fx1], [fy fy1],'Color','r','LineStyle','--'); %Segmento 1
fdtd2 = line([fx1 fx2], [fy1 fy2],'Color','r','LineStyle','--'); %Segmento 2
fdtd3 = line([fx2 fx3], [fy2 fy3],'Color','r','LineStyle','--'); %Segmento 3
fdtd4 = line([fx3 fx], [fy3 fy],'Color','r','LineStyle','--'); %Segmento 4


xlim([0 5]);
ylim([0 5]);

%% Localização da antena receptora e transmissora

%Coordenada ponto de transmissão
xTi = 2.5;
yTi = 2.5;
hold on
plot(xTi,yTi,'k*');
txt1 = '\ Ti';
text(xTi,yTi,txt1,'HorizontalAlignment','left');
hold off

%Coordenada ponto de recepção
xTx = 3;
s0 = distance(xTi,yTi,xTx,1.75);
points = s0/50;
tata = lambda/20;
yTx = [1.75:tata:3.25];
N = length(yTx);
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
             alfa(i,ind) = [(xTi-Ptx(i,aux))*(Ptx(j,aux)-Ptx(i,aux))+(yTi-Pty(i,aux))*(Pty(j,aux)-Pty(i,aux))]/((Ptx(j,aux)-Ptx(i,aux))^2+(Pty(j,aux)-Pty(i,aux))^2);
             Qx(i,ind) = Ptx(i,aux) + alfa(i,aux)*(Ptx(j,aux)-Ptx(i,aux));
             Qy(i,ind) = Pty(i,aux) + alfa(i,aux)*(Pty(j,aux)-Pty(i,aux)); 

%          else
%              alfa(i,ind) = [(xTi-Ptx(j,aux))*(Ptx(i,aux)-Ptx(j,aux))+(yTi-Pty(j,aux))*(Pty(i,aux)-Pty(j,aux))]/[(Ptx(i,aux)-Ptx(j,aux))^2+(Pty(i,aux)-Pty(j,aux))^2];
%              Qx(i,ind) = Ptx(j,aux) + alfa(i,aux)*(Ptx(i,aux)-Ptx(j,aux));
%              Qy(i,ind) = Pty(j,aux) + alfa(i,aux)*(Pty(i,aux)-Pty(j,aux));
%          end
         
         if ((alfa(i,ind)) >= 0 && (alfa(i,ind)) <= 1)
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
            k(~isfinite(k))=0;
            w(i,z) = ((yTx(z)-Iy(i,aux))/(xTx-Ix(i,aux)));
            beta(i,z) = (-k(i,aux)*(xTx+Ix(i,aux))+(yTx(z)+Iy(i,aux)) -2*(Pty(i,aux)-(k(i,aux)*Ptx(i,aux))))/[(k(i,aux)*((xTx-Ix(i,aux))))-(yTx(z)-Iy(i,aux))];
            beta(isnan(beta))= 0;
            beta(~isfinite(beta))= 0;
            ALFA(i,z) = (-w(i,z)*(Ptx(j,aux)+Ptx(i,aux))+(Pty(j,aux)+Pty(i,aux)) -2*(Iy(i,aux)-(w(i,z)*Ix(i,aux))))/((w(i,z)*(Ptx(j,aux)-Ptx(i,aux)))-(Pty(j,aux)-Pty(i,aux)));
            Px(i,z) = 1/2*[(xTx +Ix(i,aux))+ beta(i,z)*(xTx-Ix(i,aux))];
            Px(~isfinite(Px))= 0;
            Py(i,z) = 1/2*[(yTx(z)+Iy(i,aux))+ beta(i,z)*(yTx(z)-Iy(i,aux))];
            Py(~isfinite(Py))= 0;
            Py(isnan(Py))= 0;
            Sx(i,z) = 1/2*((Ptx(j,aux)+Ptx(i,aux))+ ALFA(i,z)*(Ptx(j,aux)-Ptx(i,aux)));
            Sy(i,z) = 1/2*((Pty(j,aux)+Pty(i,aux))+ ALFA(i,z)*(Pty(j,aux)-Pty(i,aux)));
            
%              if (abs(beta(i,z)) >= 1 )
%                 Px(i,z) = NaN;
%                 Py(i,z) = NaN;
%                 Sx(i,z) = NaN;
%                 Sy(i,z) = NaN;
%              end
%                 if (abs(ALFA(i,z)) >= 1  )  
%                 Sx(i,z) = NaN;
%                 Sy(i,z) = NaN;
%                 Px(i,z) = NaN;
%                 Py(i,z) = NaN;
%              end
            if (abs(beta(i,z)) <= 1 &&  abs(ALFA(i,z)) <= 1)
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
       if m == 4
        p = 1
        n = 4;
       end
         %if  (Ptx(n,aux) >= Ptx(p,aux)) && (Pty(n,aux) >= Pty(p,aux))
             alfa(i,m) = ((Ix(i,1)-Ptx(p,aux))*(Ptx(n,aux)-Ptx(p,aux))+(Iy(i,1)-Pty(p,aux))*(Pty(n,aux)-Pty(p,aux)))/((Ptx(n,aux)-Ptx(p,aux))^2 + (Pty(n,aux)-Pty(p,aux))^2);  
             Qx(i,m) = Ptx(p,aux) + alfa(i,m)*(Ptx(n,aux)-Ptx(p,aux));
             Qy(i,m) = Pty(p,aux) + alfa(i,m)*(Pty(n,aux)-Pty(p,aux));
%        else 
%              alfa(i,m) = ((Ix(i,1)-Ptx(n,aux))*(Ptx(p,aux)-Ptx(n,aux))+(Iy(i,1)-Pty(n,aux))*(Pty(p,aux)-Pty(n,aux)))/((Ptx(p,aux)-Ptx(n,aux))^2+(Pty(p,aux)-Pty(n,aux))^2);
%              Qx(i,m) = Ptx(n,aux) + alfa(i,m)*(Ptx(p,aux)-Ptx(n,aux));
%              Qy(i,m) = Pty(n,aux) + alfa(i,m)*(Pty(p,aux)-Pty(n,aux));
%        end
       
       if   ((alfa(i,m)) <= 1 && (alfa(i,m)) >= 0)
           
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
            w(i,4*z-4+m) = ((yTx(z)-Iy(i,m))/(xTx-Ix(i,m)));
            w(isnan(w))= 0;
            Beta(i,4*z-4+m) = ((-k(i,m)*(xTx+Ix(i,m)))+(yTx(z)+Iy(i,m)) -2*(Pty(p,aux)-(k(i,m)*Ptx(p,aux))))/((k(i,m)*((xTx-Ix(i,m))))-(yTx(z)-Iy(i,m)));
            Beta(isnan(Beta))= 0;
            Beta(~isfinite(Beta))= 0;
            Alfa(i,4*z-4+m) = ((-w(i,4*z-4+m)*(Ptx(n,aux)+Ptx(p,aux)))+(Pty(n,aux)+Pty(p,aux)) - 2*(Iy(i,m)-(w(i,4*z-4+m)*Ix(i,m))))/((w(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)))-(Pty(n,aux)-Pty(p,aux)));    
            Pontox(i,4*z-4+m) = 1/2*((xTx + Ix(i,m))+ Beta(i,4*z-4+m)*(xTx-Ix(i,m)));
            Pontox(~isfinite(Pontox))= 0;
            Pontoy(i,4*z-4+m) = 1/2*((yTx(z) + Iy(i,m))+ Beta(i,4*z-4+m)*(yTx(z)-Iy(i,m)));
            Rx(i,4*z-4+m) = 1/2*((Ptx(n,aux)+Ptx(p,aux))+ Alfa(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)));
            Ry(i,4*z-4+m) = 1/2*((Pty(n,aux)+Pty(p,aux))+ Alfa(i,4*z-4+m)*(Pty(n,aux)-Pty(p,aux)));

               
%              if (abs(Alfa(i,4*z-4+m)) > 1 ) 
%                 Rx(i,4*z-4+m) = NaN;
%                 Ry(i,4*z-4+m) = NaN;
%              end
%              if (abs(Beta(i,4*z-4+m)) > 1 )
%                 Pontoy(i,4*z-4+m) = NaN;
%                 Pontox(i,4*z-4+m) = NaN;
%                 Rx(i,4*z-4+m) = NaN;
%                 Ry(i,4*z-4+m) = NaN;
%              end
            
             if (abs(Beta(i,4*z-4+m)) <= 1 &  abs(Alfa(i,4*z-4+m)) <= 1)%(Pontox(i,4*z-4+m) == Rx(i,4*z-4+m) & Pontoy(i,4*z-4+m)== Ry(i,4*z-4+m))
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
            w1(isnan(w1))= 0;
            Beta1(i,4*z-4+m) = ((-k1(i,m)*(Ix(i,1)+Rx(i,4*z-4+m)))+(Iy(i,1)+Ry(i,4*z-4+m)) -2*(Pty(px,aux)-(k1(i,m)*Ptx(px,aux))))/((k1(i,m)*(Ix(i,1)-Rx(i,4*z-4+m)))-(Iy(i,1)-Ry(i,4*z-4+m)));
            Beta1(isnan(Beta1))= 0;
            Beta1(~isfinite(Beta1))= 0;
            Alfa1(i,4*z-4+m) = ((-w1(i,4*z-4+m)*(Ptx(py,aux)+Ptx(px,aux)))+(Pty(py,aux)+Pty(px,aux))-2*(Ry(i,4*z-4+m)-(w1(i,4*z-4+m)*Rx(i,4*z-4+m))))/((w1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)))-(Pty(py,aux)-Pty(px,aux)));
            Pontox1(i,4*z-4+m) = 1/2*((Ix(i,1)+Rx(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Ix(i,1)-Rx(i,4*z-4+m)));
            Pontox1(~isfinite(Pontox1))= 0;
            Pontoy1(i,4*z-4+m) = 1/2*((Iy(i,1)+Ry(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Iy(i,1)-Ry(i,4*z-4+m)));
            Rx1(i,4*z-4+m) = 1/2*((Ptx(py,aux)+Ptx(px,aux))+ Alfa1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)));
            Ry1(i,4*z-4+m) = 1/2*((Pty(py,aux)+Pty(px,aux))+ Alfa1(i,4*z-4+m)*(Pty(py,aux)-Pty(px,aux)));
            
%             if (abs(Alfa1(i,4*z-4+m)) > 1 )
%                 Rx1(i,4*z-4+m) = NaN;
%                 Ry1(i,4*z-4+m) = NaN;
%             end
%             if (abs(Beta1(i,4*z-4+m)) > 1)
%                 Pontoy1(i,4*z-4+m) = NaN;
%                 Pontox1(i,4*z-4+m) = NaN;
%                 Rx1(i,4*z-4+m) = NaN;
%                 Ry1(i,4*z-4+m) = NaN;
%              end
            
             if (abs(Beta1(i,4*z-4+m)) <= 1 &  abs(Alfa1(i,4*z-4+m)) <= 1)%(Pontox1(i,4*z-4+m) == Rx1(i,4*z-4+m) & Pontoy1(i,4*z-4+m) == Ry1(i,4*z-4+m))
                 Qx11(i,4*z-4+m) = Pontox1(i,4*z-4+m);
                 Qy11(i,4*z-4+m) = Pontoy1(i,4*z-4+m);
                 
             else
                 Qx11(i,4*z-4+m) = NaN;
                 Qy11(i,4*z-4+m) = NaN;
             end
             
   hold on
   %plot(Qx11(1,4*z-4+m),Qy11(1,4*z-4+m),'.');
   %plot(Pontox1(1,4*z-4+m),Pontoy1(1,4*z-4+m),'.'); % Ponto de reflexão em cada segmento
   %plot(Rx1(1,4*z-4+m),Ry1(1,4*z-4+m),'.');
   hold off
 
 %----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%%
          
           
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2º SEGMENTO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    elseif i == 2
            p = m + 1;
            n = p + 1;
       if p == 4;
            p = 1;
            n = 4;
       elseif p == 5
           p = 1;
           n = 2;
       end

                %if  (Ptx(n,aux) >= Ptx(p,aux)) && (Pty(n,aux) >= Pty(p,aux))
             alfa(i,m) = ((Ix(i,1)-Ptx(p,aux))*(Ptx(n,aux)-Ptx(p,aux))+(Iy(i,1)-Pty(p,aux))*(Pty(n,aux)-Pty(p,aux)))/((Ptx(n,aux)-Ptx(p,aux))^2+(Pty(n,aux)-Pty(p,aux))^2);  
             Qx(i,m) = Ptx(p,aux) + alfa(i,m)*(Ptx(n,aux)-Ptx(p,aux));
             Qy(i,m) = Pty(p,aux) + alfa(i,m)*(Pty(n,aux)-Pty(p,aux));
%        else 
%              alfa(i,m) = ((Ix(i,1)-Ptx(n,aux))*(Ptx(p,aux)-Ptx(n,aux))+(Iy(i,1)-Pty(n,aux))*(Pty(p,aux)-Pty(n,aux)))/((Ptx(p,aux)-Ptx(n,aux))^2+(Pty(p,aux)-Pty(n,aux))^2);
%              Qx(i,m) = Ptx(n,aux) + alfa(i,m)*(Ptx(p,aux)-Ptx(n,aux));
%              Qy(i,m) = Pty(n,aux) + alfa(i,m)*(Pty(p,aux)-Pty(n,aux));
%        end
%              alfa(i,m) = ((Ix(i,1)-Ptx(p,aux))*(Ptx(n,aux)-Ptx(p,aux))+(Iy(i,1)-Pty(p,aux))*(Pty(n,aux)-Pty(p,aux)))/((Ptx(n,aux)-Ptx(p,aux))^2+(Pty(n,aux)-Pty(p,aux))^2);  
%              Qx(i,m) = Ptx(p,aux) + alfa(i,m)*(Ptx(n,aux)-Ptx(p,aux));
%              Qy(i,m) = Pty(p,aux) + alfa(i,m)*(Pty(n,aux)-Pty(p,aux));

           
       if ((alfa(i,m)) <= 1 && (alfa(i,m)) >= 0)
           
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
            w(i,4*z-4+m) = ((yTx(z)-Iy(i,m))/(xTx-Ix(i,m)));
            w(isnan(w))= 0;
            Beta(i,4*z-4+m) = ((-k(i,m)*(xTx+Ix(i,m)))+(yTx(z)+Iy(i,m)) -2*(Pty(p,aux)-(k(i,m)*Ptx(p,aux))))/((k(i,m)*((xTx-Ix(i,m))))-(yTx(z)-Iy(i,m)));
            Beta(isnan(Beta))= 0;
            Beta(~isfinite(Beta))= 0;
            Alfa(i,4*z-4+m) = ((-w(i,4*z-4+m)*(Ptx(n,aux)+Ptx(p,aux)))+(Pty(n,aux)+Pty(p,aux)) -2*(Iy(i,m)-(w(i,4*z-4+m)*Ix(i,m))))/(w(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux))-(Pty(n,aux)-Pty(p,aux)));   
            Pontox(i,4*z-4+m) = 1/2*[(xTx + Ix(i,m))+ Beta(i,4*z-4+m)*(xTx-Ix(i,m))];
            Pontox(~isfinite(Pontox))= 0;
            Pontoy(i,4*z-4+m) = 1/2*[(yTx(z)+Iy(i,m))+ Beta(i,4*z-4+m)*(yTx(z)-Iy(i,m))];
            Rx(i,4*z-4+m) = 1/2*((Ptx(n,aux)+Ptx(p,aux))+ Alfa(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)));
            Ry(i,4*z-4+m) = 1/2*((Pty(n,aux)+Pty(p,aux))+ Alfa(i,4*z-4+m)*(Pty(n,aux)-Pty(p,aux)));
            
%             if (abs(Alfa(i,4*z-4+m)) > 1 ) 
%                 Rx(i,4*z-4+m) = NaN;
%                 Ry(i,4*z-4+m) = NaN;
%              end
%              if (abs(Beta(i,4*z-4+m)) > 1)
%                 Pontoy(i,4*z-4+m) = NaN;
%                 Pontox(i,4*z-4+m) = NaN;
%                 Rx(i,4*z-4+m) = NaN;
%                 Ry(i,4*z-4+m) = NaN;
%             end
            
            if (abs(Beta(i,4*z-4+m)) <= 1 &&  abs(Alfa(i,4*z-4+m)) <= 1)%(Pontox(i,4*z-4+m) == Rx(i,4*z-4+m) & Pontoy(i,4*z-4+m)== Ry(i,4*z-4+m))
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
            Beta1(i,4*z-4+m) = ((-k1(i,m)*(Ix(i,1)+Rx(i,4*z-4+m)))+(Iy(i,1)+Ry(i,4*z-4+m)) -2*(Pty(px,aux)-(k1(i,m)*Ptx(px,aux))))/((k1(i,m)*(Ix(i,1)-Rx(i,4*z-4+m)))-(Iy(i,1)-Ry(i,4*z-4+m)));
            Beta1(isnan(Beta1))= 0;
            Beta1(~isfinite(Beta1))= 0;
            Alfa1(i,4*z-4+m) = ((-w1(i,4*z-4+m)*(Ptx(py,aux)+Ptx(px,aux)))+(Pty(py,aux)+Pty(px,aux)) -2*(Ry(i,4*z-4+m)-(w1(i,4*z-4+m)*Rx(i,4*z-4+m))))/((w1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)))-(Pty(py,aux)-Pty(px,aux)));
            Pontox1(i,4*z-4+m) = 1/2*[(Ix(i,1)+Rx(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Ix(i,1)-Rx(i,4*z-4+m))];
            Pontox1(~isfinite(Pontox))= 0;
            Pontoy1(i,4*z-4+m) = 1/2*[(Iy(i,1)+Ry(i,4*z-4+m)) + Beta1(i,4*z-4+m)*(Iy(i,1)-Ry(i,4*z-4+m))];
            Rx1(i,4*z-4+m) = 1/2*((Ptx(py,aux)+Ptx(px,aux))+ Alfa1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)));
            Ry1(i,4*z-4+m) = 1/2*((Pty(py,aux)+Pty(px,aux))+ Alfa1(i,4*z-4+m)*(Pty(py,aux)-Pty(px,aux)));
                        
%             if (abs(Alfa1(i,4*z-4+m)) > 1 )
%                 Rx1(i,4*z-4+m) = NaN;
%                 Ry1(i,4*z-4+m) = NaN;
%             end
%             if (abs(Beta1(i,4*z-4+m)) > 1)
%                 Pontoy1(i,4*z-4+m) = NaN;
%                 Pontox1(i,4*z-4+m) = NaN;
%                 Rx1(i,4*z-4+m) = NaN;
%                 Ry1(i,4*z-4+m) = NaN;
%              end
            
             if (abs(Beta1(i,4*z-4+m)) <= 1 &&  abs(Alfa1(i,4*z-4+m)) <= 1)%(Pontox1(i,4*z-4+m) == Rx1(i,4*z-4+m) & Pontoy1(i,4*z-4+m) == Ry1(i,4*z-4+m))
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

                       
        if ((alfa(i,m)) <= 1  && (alfa(i,m)) >= 0)
           
             Ix(i,m) = 2*Qx(i,m) - Ix(i,1);
             Iy(i,m) = 2*Qy(i,m) - Iy(i,1);
       else
              Ix(i,m) = NaN;
              Iy(i,m) = NaN;
        end
             

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculos dos pontos de reflexões DA 3 IMAGEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            k(i,m) = (Pty(n,aux)-Pty(p,aux))/(Ptx(n,aux)-Ptx(p,aux));
            k(~isfinite(k))= 0;
            w(i,4*z-4+m) = ((yTx(z)-Iy(i,m))/(xTx-Ix(i,m)));
            w(isnan(w))= 0;
            Beta(i,4*z-4+m) = ((-k(i,m)*(xTx+Ix(i,m)))+(yTx(z)+Iy(i,m)) -2*(Pty(p,aux)-(k(i,m)*Ptx(p,aux))))/((k(i,m)*((xTx-Ix(i,m))))-(yTx(z)-Iy(i,m)));
            Beta(isnan(Beta))= 0;
            Beta(~isfinite(Beta))= 0;
            Alfa(i,4*z-4+m) = ((-w(i,4*z-4+m)*(Ptx(n,aux)+Ptx(p,aux)))+(Pty(n,aux)+Pty(p,aux)) -2*(Iy(i,m)-(w(i,4*z-4+m)*Ix(i,m))))/((w(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)))-(Pty(n,aux)-Pty(p,aux)));
            Pontox(i,4*z-4+m) = 1/2*[(xTx + Ix(i,m))+ Beta(i,4*z-4+m)*(xTx-Ix(i,m))];
            Pontox(~isfinite(Pontox))= 0;
            Pontoy(i,4*z-4+m) = 1/2*[(yTx(z)+Iy(i,m))+ Beta(i,4*z-4+m)*(yTx(z)-Iy(i,m))];
            Rx(i,4*z-4+m) = 1/2*((Ptx(n,aux)+Ptx(p,aux))+ Alfa(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)));
            Ry(i,4*z-4+m) = 1/2*((Pty(n,aux)+Pty(p,aux))+ Alfa(i,4*z-4+m)*(Pty(n,aux)-Pty(p,aux)));

            
%             if (abs(Alfa(i,4*z-4+m)) > 1) 
%                 Rx(i,4*z-4+m) = NaN;
%                 Ry(i,4*z-4+m) = NaN;
%              end
%              if (abs(Beta(i,4*z-4+m)) > 1)
%                 Pontoy(i,4*z-4+m) = NaN;
%                 Pontox(i,4*z-4+m) = NaN;
%                 Rx(i,4*z-4+m) = NaN;
%                 Ry(i,4*z-4+m) = NaN;
%             end
            
            if (abs(Beta(i,4*z-4+m)) <= 1 &&  abs(Alfa(i,4*z-4+m)) <= 1)%(Pontox(i,4*z-4+m) == Rx(i,4*z-4+m) & Pontoy(i,4*z-4+m)== Ry(i,4*z-4+m))
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
 
 px = 3;
 py = 4;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da imagem de Primeira ordem para ponto de recepção %%%%%%%%%%%%%%%%%%%
 
            k1(i,m) = (Pty(py,aux)-Pty(px,aux))/(Ptx(py,aux)-Ptx(px,aux));
            k1(~isfinite(k1))= 0;
            w1(i,4*z-4+m) = ((Iy(i,1)-Ry(i,4*z-4+m))/(Ix(i,1)-Rx(i,4*z-4+m)));
            Beta1(i,4*z-4+m) = ((-k1(i,m)*(Ix(i,1)+Rx(i,4*z-4+m)))+(Iy(i,1)+Ry(i,4*z-4+m)) -2*(Pty(px,aux)-(k1(i,m)*Ptx(px,aux))))/((k1(i,m)*(Ix(i,1)-Rx(i,4*z-4+m)))-(Iy(i,1)-Ry(i,4*z-4+m)));
            Beta1(isnan(Beta1))= 0;
            Beta1(~isfinite(Beta1))= 0;
            Alfa1(i,4*z-4+m) = ((-w1(i,4*z-4+m)*(Ptx(py,aux)+Ptx(px,aux)))+(Pty(py,aux)+Pty(px,aux)) -2*(Ry(i,4*z-4+m)-(w1(i,4*z-4+m)*Rx(i,4*z-4+m))))/((w1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)))-(Pty(py,aux)-Pty(px,aux)));
            Pontox1(i,4*z-4+m) = 1/2*[(Ix(i,1)+Rx(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Ix(i,1)-Rx(i,4*z-4+m))];
            Pontox1(~isfinite(Pontox1))= 0;
            Pontoy1(i,4*z-4+m) = 1/2*[(Iy(i,1)+Ry(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Iy(i,1)-Ry(i,4*z-4+m))];
            Rx1(i,4*z-4+m) = 1/2*((Ptx(py,aux)+Ptx(px,aux))+ Alfa1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)));
            Ry1(i,4*z-4+m) = 1/2*((Pty(py,aux)+Pty(px,aux))+ Alfa1(i,4*z-4+m)*(Pty(py,aux)-Pty(px,aux)));
  
%  
%                         if (abs(Alfa1(i,4*z-4+m)) > 1)
%                             Rx1(i,4*z-4+m) = NaN;
%                             Ry1(i,4*z-4+m) = NaN;
%                             
%                         end
%                         if (abs(Beta1(i,4*z-4+m)) > 1)
%                             Pontoy1(i,4*z-4+m) = NaN;
%                             Pontox1(i,4*z-4+m) = NaN;
%                             Rx1(i,4*z-4+m) = NaN;
%                             Ry1(i,4*z-4+m) = NaN;
%                          end

                         if (abs(Beta1(i,4*z-4+m)) <= 1 &  abs(Alfa1(i,4*z-4+m)) <= 1)%(Pontox1(i,4*z-4+m) == Rx1(i,4*z-4+m) & Pontoy1(i,4*z-4+m) == Ry1(i,4*z-4+m))
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
       if m == 2
        p = 1;
        n = 2;
       elseif m == 3
        p = 2;
        n = 3;
       elseif m == 4
        p = 3;
        n = 4;
       end
           
%        if  (Ptx(n,aux) >= Ptx(p,aux)) && (Pty(n,aux) >= Pty(p,aux))
             alfa(i,m) = ((Ix(i,1)-Ptx(p,aux))*(Ptx(n,aux)-Ptx(p,aux))+(Iy(i,1)-Pty(p,aux))*(Pty(n,aux)-Pty(p,aux)))/((Ptx(n,aux)-Ptx(p,aux))^2+(Pty(n,aux)-Pty(p,aux))^2);  
             Qx(i,m) = Ptx(p,aux) + alfa(i,m)*(Ptx(n,aux)-Ptx(p,aux));
             Qy(i,m) = Pty(p,aux) + alfa(i,m)*(Pty(n,aux)-Pty(p,aux));
%        else 
%              alfa(i,m) = ((Ix(i,1)-Ptx(n,aux))*(Ptx(p,aux)-Ptx(n,aux))+(Iy(i,1)-Pty(n,aux))*(Pty(p,aux)-Pty(n,aux)))/((Ptx(p,aux)-Ptx(n,aux))^2+(Pty(p,aux)-Pty(n,aux))^2);
%              Qx(i,m) = Ptx(n,aux) + alfa(i,m)*(Ptx(p,aux)-Ptx(n,aux));
%              Qy(i,m) = Pty(n,aux) + alfa(i,m)*(Pty(p,aux)-Pty(n,aux));
%        end

        
             
             
         if ((alfa(i,m)) <= 1  && (alfa(i,m)) >= 0)
           
             Ix(i,m) = 2*Qx(i,m) - Ix(i,1);
             Iy(i,m) = 2*Qy(i,m) - Iy(i,1);
         else
              Ix(i,m) = NaN;
              Iy(i,m) = NaN;
        end
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% da imagem de segunda ordem para ponto de recepção %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            k(i,m) = (Pty(n,aux)-Pty(p,aux))/(Ptx(n,aux)-Ptx(p,aux));
            k(~isfinite(k))= 0;
            w(i,4*z-4+m) = ((yTx(z)-Iy(i,m))/(xTx-Ix(i,m)));
            w(isnan(w))= 0;
            Beta(i,4*z-4+m) = ((-k(i,m)*(xTx+Ix(i,m)))+(yTx(z)+Iy(i,m)) -2*(Pty(p,aux)-(k(i,m)*Ptx(p,aux))))/((k(i,m)*((xTx-Ix(i,m))))-(yTx(z)-Iy(i,m)));
            Beta(isnan(Beta))= 0;
            Beta(~isfinite(Beta))= 0;
            Alfa(i,4*z-4+m) = ((-w(i,4*z-4+m)*(Ptx(n,aux)+Ptx(p,aux)))+(Pty(n,aux)+Pty(p,aux)) -2*(Iy(i,m)-(w(i,4*z-4+m)*Ix(i,m))))/(w(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux))-(Pty(n,aux)-Pty(p,aux)));
              
            Pontox(i,4*z-4+m) = 1/2*[(xTx + Ix(i,m))+ Beta(i,4*z-4+m)*(xTx-Ix(i,m))];
            Pontox(~isfinite(Pontox))= 0;
            Pontoy(i,4*z-4+m) = 1/2*[(yTx(z)+Iy(i,m))+ Beta(i,4*z-4+m)*(yTx(z)-Iy(i,m))];
            Rx(i,4*z-4+m) = 1/2*((Ptx(n,aux)+Ptx(p,aux))+ Alfa(i,4*z-4+m)*(Ptx(n,aux)-Ptx(p,aux)));
            Ry(i,4*z-4+m) = 1/2*((Pty(n,aux)+Pty(p,aux))+ Alfa(i,4*z-4+m)*(Pty(n,aux)-Pty(p,aux)));
               
            
%             if (abs(Alfa(i,4*z-4+m)) > 1 ) 
%                 Rx(i,4*z-4+m) = NaN;
%                 Ry(i,4*z-4+m) = NaN;
%              end
%              if (abs(Beta(i,4*z-4+m)) > 1)
%                 Pontoy(i,4*z-4+m) = NaN;
%                 Pontox(i,4*z-4+m) = NaN;
%                 Rx(i,4*z-4+m) = NaN;
%                 Ry(i,4*z-4+m) = NaN;
%             end
            
            if (abs(Beta(i,4*z-4+m)) <= 1 &  abs(Alfa(i,4*z-4+m)) <= 1)%(Pontox(i,4*z-4+m) == Rx(i,4*z-4+m) & Pontoy(i,4*z-4+m)== Ry(i,4*z-4+m))
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
            Beta1(i,4*z-4+m) = ((-k1(i,m)*(Ix(i,1)+Rx(i,4*z-4+m)))+(Iy(i,1)+Ry(i,4*z-4+m)) -2*(Pty(px,aux)-(k1(i,m)*Ptx(px,aux))))/((k1(i,m)*(Ix(i,1)-Rx(i,4*z-4+m)))-(Iy(i,1)-Ry(i,4*z-4+m)));
            Beta1(isnan(Beta1))= 0;
            Beta1(~isfinite(Beta1))= 0;
            Alfa1(i,4*z-4+m) = ((-w1(i,4*z-4+m)*(Ptx(py,aux)+Ptx(px,aux)))+(Pty(py,aux)+Pty(px,aux)) -2*(Ry(i,4*z-4+m)-(w1(i,4*z-4+m)*Rx(i,4*z-4+m))))/((w1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)))-(Pty(py,aux)-Pty(px,aux)));
            Pontox1(i,4*z-4+m) = 1/2*[(Ix(i,1)+Rx(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Ix(i,1)-Rx(i,4*z-4+m))];
            Pontox1(~isfinite(Pontox))= 0;
            Pontoy1(i,4*z-4+m) = 1/2*[(Iy(i,1)+Ry(i,4*z-4+m))+ Beta1(i,4*z-4+m)*(Iy(i,1)-Ry(i,4*z-4+m))];
            Rx1(i,4*z-4+m) = 1/2*((Ptx(py,aux)+Ptx(px,aux))+ Alfa1(i,4*z-4+m)*(Ptx(py,aux)-Ptx(px,aux)));
            Ry1(i,4*z-4+m) = 1/2*((Pty(py,aux)+Pty(px,aux))+ Alfa1(i,4*z-4+m)*(Pty(py,aux)-Pty(px,aux)));
  
 

%                         if (abs(Alfa1(i,4*z-4+m)) > 1 )
%                             Rx1(i,4*z-4+m) = NaN;
%                             Ry1(i,4*z-4+m) = NaN;
%                         end
%                         if (abs(Beta1(i,4*z-4+m)) > 1)
%                             Pontoy1(i,4*z-4+m) = NaN;
%                             Pontox1(i,4*z-4+m) = NaN;
%                             Rx1(i,4*z-4+m) = NaN;
%                             Ry1(i,4*z-4+m) = NaN;
%                          end

                         if (abs(Beta1(i,4*z-4+m)) <= 1 &&  abs(Alfa1(i,4*z-4+m)) <= 1)%(Pontox1(i,4*z-4+m) == Rx1(i,4*z-4+m) & Pontoy1(i,4*z-4+m) == Ry1(i,4*z-4+m))
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

        if (alfa3(i,q) <= 1  && alfa3(i,q) >= 0)
           
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
            w3(i,9*z-9+q) = ((Iy3(i,q)-yTx(z))/(Ix3(i,q)-xTx));
            Beta3(i,9*z-9+q) = ((-k3(i,q)*(xTx+Ix3(i,q)))+(yTx(z)+Iy3(i,q)) -(2*(Pty(p1(i,q),aux)-(k3(i,q)*Ptx(p1(i,q),aux)))))/((k3(i,q)*((Ix3(i,q)-xTx)))-(Iy3(i,q)-yTx(z)));
            Beta3(isnan(Beta3))= 0;
            Beta3(~isfinite(Beta3))= 0;
            Alfa3(i,9*z-9+q) = (-w3(i,9*z-9+q)*(Ptx(p2(i,q),aux)+Ptx(p1(i,q),aux))+(Pty(p2(i,q),aux)+Pty(p1(i,q),aux)) - (2*(yTx(z)-(w3(i,9*z-9+q)*xTx))))/((w3(i,9*z-9+q)*(Ptx(p2(i,q),aux)-Ptx(p1(i,q),aux)))-(Pty(p2(i,q),aux)-Pty(p1(i,q),aux)));
            Pontox3(i,9*z-9+q) = 1/2*[(xTx + Ix3(i,q))+ Beta3(i,9*z-9+q)*(Ix3(i,q)-xTx)];
            Pontox3(~isfinite(Pontox3))= 0;
            Pontoy3(i,9*z-9+q) = 1/2*[(yTx(z)+Iy3(i,q))+ Beta3(i,9*z-9+q)*(Iy3(i,q)-yTx(z))];
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

         if (alfa4(i,c) <= 1  && alfa4(i,c) >= 0)
           
             Ix4(i,c) = 2*Qx4(i,c) - Ix3(i,q);
             Iy4(i,c) = 2*Qy4(i,c) - Iy3(i,q);
        else
              Ix4(i,c) = NaN;
              Iy4(i,c) = NaN;
             
        end
        
        
            k4(i,c) = (Pty(ptx1(i,c),aux)-Pty(ptx0(i,c),aux))/(Ptx(ptx1(i,c),aux)-Ptx(ptx0(i,c),aux));
            k4(~isfinite(k4))= 0;
            w4(i,27*z-27+c) = ((Iy4(i,c)-yTx(z))/(Ix4(i,c)-xTx));
            Beta4(i,27*z-27+c) = ((-k4(i,c)*(xTx+Ix4(i,c)))+(yTx(z)+Iy4(i,c)) -(2*(Pty(ptx0(i,c),aux)-(k4(i,c)*Ptx(ptx0(i,c),aux)))))/((k4(i,c)*((Ix4(i,c)-xTx)))-(Iy4(i,c)-yTx(z)));
            Beta4(isnan(Beta4))= 0;
            Beta4(~isfinite(Beta4))= 0;
            Alfa4(i,27*z-27+c) = (-w4(i,27*z-27+c)*(Ptx(ptx1(i,c),aux)+Ptx(ptx0(i,c),aux))+(Pty(ptx1(i,c),aux)+Pty(ptx0(i,c),aux)) - (2*(yTx(z)-(w4(i,27*z-27+c)*xTx))))/((w4(i,27*z-27+c)*(Ptx(ptx1(i,c),aux)-Ptx(ptx0(i,c),aux)))-(Pty(ptx1(i,c),aux)-Pty(ptx0(i,c),aux)));
            Pontox4(i,27*z-27+c) = 1/2*[(xTx + Ix4(i,c))+ Beta4(i,27*z-27+c)*(Ix4(i,c)-xTx)];
            Pontox4(~isfinite(Pontox4))= 0;
            Pontoy4(i,27*z-27+c) = 1/2*[(yTx(z)+Iy4(i,c))+ Beta4(i,27*z-27+c)*(Iy4(i,c)-yTx(z))];
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
         plot(Rx44(i,27*z-27+c),Ry44(i,27*z-27+c),'.');
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
for Tst = 1:N
% Tst = N;
 line ([xTi,xTx],[yTi,yTx(Tst)]);
% 
for i = 1:Ns
     if (abs(ALFA(i,Tst)) < 1 && abs(beta(i,Tst)) < 1 )
     line ([xTi,Sx(i,Tst),xTx],[yTi,Sy(i,Tst),yTx(Tst)]);
     end
    if Nr > 1
        for m = 2:Elementos
            if ((abs(Alfa(i,4*Tst-4+m)) < 1   && abs(Alfa1(i,4*Tst-4+m)) < 1) && (abs(Beta(i,4*Tst-4+m)) < 1 && abs(Beta1(i,4*Tst-4+m)) < 1)  )
             line ([xTi,Rx1(i,4*Tst-4+m),Rx(i,4*Tst-4+m),xTx],[yTi,Ry1(i,4*Tst-4+m),Ry(i,4*Tst-4+m),yTx(Tst)]);
            end
                    if Nr > 2
                     for q = 2*(m-2)+(m-1) : 2*(m-2)+ m + 1 
                             if (abs(Alfa31(i,9*Tst-9+q))< 1  && abs(Alfa32(i,9*Tst-9+q)) < 1  && abs(Alfa3(i,9*Tst-9+q))< 1 )
                              line ([xTi,Rx31(i,9*Tst-9+q),Rx32(i,9*Tst-9+q),Rx3(i,9*Tst-9+q),xTx],[yTi,Ry31(i,9*Tst-9+q),Ry32(i,9*Tst-9+q),Ry3(i,9*Tst-9+q),yTx(Tst)]);
                             end
                             if Nr > 3
                                        for c = 2*(q-1) + q : 2*(q-1) + q + 2
                                            if (abs(Alfa4(i,27*Tst-27+c))< 1  && abs(Alfa43(i,27*Tst-27+c)) < 1  && abs(Alfa40(i,27*Tst-27+c))< 1 && abs(Alfa44(i,27*Tst-27+c))< 1 )
                                            line ([xTi,Rx40(i,27*Tst-27+c),Rx44(i,27*Tst-27+c),Rx43(i,27*Tst-27+c),Rx4(i,27*Tst-27+c),xTx],[yTi,Ry40(i,27*Tst-27+c),Ry44(i,27*Tst-27+c),Ry43(i,27*Tst-27+c),Ry4(i,27*Tst-27+c),yTx(Tst)])
                                            end
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
rho = 0.01;
Pt = 1;
Gt = 1;
E0 = 100;%abs(sqrt((N0*Pt*Gt)/(2*pi*rho*rho))); % Ponto de Excitação
H0 = E0/N0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Raio direto
 for z = 1:N
    im = sqrt(-1);
    dist(z) = distance(xTi,yTi,xTx,yTx(z));
    A(z) = sqrt(rho/dist(z)); %Fator de espalhamento
    
%     if dist(z)<= rho 
%         Ein(z) = 0;
%     else   
    Ein(z) = E0*A(z)*(exp(-im*B*dist(z))); %Campo elétrico referente na antena receptora
    Ein(~isfinite(Ein))= 0;

%     end
 end
 
% Campo para primeira reflexão

 for z = 1:N
  for i = 1:Ns
  if ( abs(beta(i,z)) < 1 && abs(ALFA(i,z)) < 1   )   
     im = sqrt(-1);
     dist1(i,z) = distance(xTi,yTi,Sx(i,z),Sy(i,z)); %s
     %dist1(i,z) = distance(Ix(i,1),Iy(i,1),Sx(i,z),Sy(i,z)); %s
     %distI(i,z) = distance(Sx(i,z),Sy(i,z),xTx,yTx(z)); %s'
     Esp(i,z) = sqrt(rho/dist1(i,z)); %Fator de espalhamento
     Ein0(i,z) = E0*Esp(i,z)*(exp(-im*B*dist1(i,z)));
     
%      %%%%%%%%%%%%%%%%%%%angulos
     
     m1(i,z) = (Sy(i,z) - yTi)/(Sx(i,z) - xTi);
     m2(i,z) = (Sy(i,z)-yTx(z))/(Sx(i,z) - xTx);
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
     
     coefref1(i,z) = (((Ninc*cosTheta(i,z)- N0*an2(i,z))/(Ninc*cosTheta(i,z) + N0*an2(i,z)))); %% COmponente paralela.
     coefperp(i,z) = (((N0*cosTheta(i,z)- Ninc*an2(i,z))/(N0*cosTheta(i,z) + Ninc*an2(i,z)))); %% COmponente perpendicular.
     
%      coefrefl(i,z) = ((N0*cosTheta(i,z)- Ninc*an2(i,z))/(N0*cosTheta(i,z) + Ninc*an2(i,z)));
     coefref10(i,z) =  abs(coefref1(i,z));
     coefperp1(i,z) = abs(coefperp(i,z)); %componente perpendicular
     
     Eir(i,z) = coefref10(i,z)*(Ein0(i,z));
     Eip(i,z) = coefperp1(i,z)*(Ein0(i,z)); %componente perpendicular
     Eir(isnan(Eir))= 0;
     Eip(isnan(Eip))= 0;
     dist2(i,z) = distance(Sx(i,z),Sy(i,z),xTx,yTx(z));
     pIm(i,z) = distance(Ix(i,1),Iy(i,1),Sx(i,z),Sy(i,z));
     Esp1(i,z) = sqrt(pIm(i,z)/(dist2(i,z) + pIm(i,z))); %Fator de espalhamento
     Esp1(~isfinite(Esp1))= 0;
     
     Ein1(i,z) = (Eir(i,z))*Esp1(i,z)*exp(-im*B*dist2(i,z));
    % Ein1(i,z) =  (Ein0(i,z))*(coefref10(i,z) + coefperp1(i,z))*Esp1(i,z)*exp(-im*B*dist2(i,z)); % levando em consideração 2 componentes
     Ein1(~isfinite(Ein1))= 0;
     distTx(z) = distance(xTi,yTi,xTx,yTx(z));
     Ein1(~isfinite(Ein1))= 0;
     
%       if (distTx(z)<= rho || abs(ALFA(i,z)) > 1 || abs(beta(i,z)) > 1 )
  else     
        Ein1(i,z) = 0;
     
        
      end
  end
 end
 Ein1(3,:) = 0;
 
%  CAMPO PARA 2 REFLEXÕES     
for z = 1:N
    for i = 1:Ns   
        for m = 2:Elementos                    
               if (abs(Alfa(i,4*z-4+m)) < 1   && abs(Alfa1(i,4*z-4+m)) < 1 && abs(Beta(i,4*z-4+m)) < 1 && abs(Beta1(i,4*z-4+m)) < 1  )
                            im = sqrt(-1);
                            distR(i,4*z-4+m) = distance(xTi,yTi,Rx1(i,4*z-4+m),Ry1(i,4*z-4+m)); %s
                            distRI(i,4*z-4+m) = distance(xTi,yTi,xTx,yTx(z)); %s'
                           % EspR(i,4*z-4+m) = sqrt(distR(i,4*z-4+m)/(distRI(i,4*z-4+m)+distR(i,4*z-4+m))); %Fator de espalhamento
                            EspR(i,4*z-4+m) = sqrt(rho/distR(i,4*z-4+m));
                            EspR(isnan(EspR))= 0;
                            Einput(i,4*z-4+m) = E0*EspR(i,4*z-4+m)*(exp(-im*B*distR(i,4*z-4+m)));
                            Einput(isnan(Einput))= 0;
                            
                     %%%%%%%%%%%%%%%%%%%%%%%%%angulos
                            m11(i,4*z-4+m) = (Ry1(i,4*z-4+m) - yTi)/(Rx1(i,4*z-4+m) - xTi);
                            m11(isnan(m11))= 0;
                            m3(i,4*z-4+m) = (Ry(i,4*z-4+m) - Ry1(i,4*z-4+m))/(Rx(i,4*z-4+m) - Rx1(i,4*z-4+m));
                            m3(isnan(m3))= 0;
                            ANG0(i,4*z-4+m)= abs((m3(i,4*z-4+m)-  m11(i,4*z-4+m))/(1 + m11(i,4*z-4+m)*m3(i,4*z-4+m)));
                            ANG0(isnan(ANG0))= 0;
                            TETA0(i,4*z-4+m)= atand(ANG0(i,4*z-4+m));
                            aux2(i,4*z-4+m) = TETA0(i,4*z-4+m)/2;
                            aux2(isnan(aux2))= 0;
                            cosTheta1(i,4*z-4+m) = (sqrt(1-((((u0*e0)/(uc*ec))^2)*sind(aux2(i,4*z-4+m))*sind(aux2(i,4*z-4+m)))));
                            angulos1(i,4*z-4+m) = cosd(aux2(i,4*z-4+m));
                            angulos1(isnan(angulos1))= 0;
                             
                     %%%%%%%%%%%%%%%%%%%%%%%%%%
                     
                  %  Nincn1(i,4*z-4+m) = (Ninc/cosTheta1(i,4*z-4+m))*[(N0/angulos1(i,4*z-4+m))+(Ninc/cosTheta1(i,4*z-4+m))*tanh(Yc*d*cosTheta(i,z))]/((Ninc/cosTheta1(i,4*z-4+m))+(N0/angulos1(i,4*z-4+m))*tanh(Yc*d*cosTheta1(i,4*z-4+m))); %Componente normal da impedância de entrada
                 %   Nincp1(i,4*z-4+m) = (Ninc*cosTheta1(i,4*z-4+m))*[(N0*angulos1(i,4*z-4+m))+(Ninc*cosTheta1(i,4*z-4+m))*tanh(Yc*d*cosTheta(i,z))]/[(Ninc*cosTheta1(i,4*z-4+m))+(N0*angulos1(i,4*z-4+m))*tanh(Yc*d*cosTheta1(i,4*z-4+m))]; %Componente paralela da impedância de entrada
                   % coefref20(i,4*z-4+m) =  (Nincn1(i,4*z-4+m) - (N0*angulos1(i,4*z-4+m)))/(Nincp1(i,4*z-4+m) + (N0*angulos1(i,4*z-4+m)));

                    coefref20(i,4*z-4+m) = ([Ninc*cosTheta1(i,4*z-4+m)- N0*angulos1(i,4*z-4+m)]/[Ninc*cosTheta1(i,4*z-4+m) + N0*angulos1(i,4*z-4+m)]);
                     coefref2(i,4*z-4+m)=  abs(coefref20(i,4*z-4+m));
                     
                     %%%%%%%%%Componente perpendicular
                     coefp20(i,4*z-4+m) = ([N0*cosTheta1(i,4*z-4+m)- Ninc*angulos1(i,4*z-4+m)]/[N0*cosTheta1(i,4*z-4+m) + Ninc*angulos1(i,4*z-4+m)]);
                     coefp2(i,4*z-4+m)=  abs(coefp20(i,4*z-4+m));
                     
                     Eirp(i,4*z-4+m) = coefp2(i,4*z-4+m)*(Einput(i,4*z-4+m)); % componente perpendicular
                     Eir1(i,4*z-4+m) = coefref2(i,4*z-4+m)*(Einput(i,4*z-4+m));
                     Eir1(isnan(Eir1))= 0;
                     pIm1(i,4*z-4+m) = distance(Ix(i,1),Iy(i,1),Rx1(i,4*z-4+m),Ry1(i,4*z-4+m));
                     dist3(i,4*z-4+m) = distance(Rx1(i,4*z-4+m),Ry1(i,4*z-4+m),Rx(i,4*z-4+m),Ry(i,4*z-4+m));
                     Esp2(i,4*z-4+m) = sqrt(pIm1(i,4*z-4+m)/(dist3(i,4*z-4+m) + pIm1(i,4*z-4+m))); %Fator de espalhamento
                     Esp2(isnan(Esp2))= 0;
                     Eir2(i,4*z-4+m) = Eir1(i,4*z-4+m)*Esp2(i,4*z-4+m)*exp(-im*B*dist3(i,4*z-4+m));
                     Eir2(isnan(Eir2))= 0;
                     %Eir2(i,4*z-4+m) = (Einput(i,4*z-4+m))*(coefref2(i,4*z-4+m) + coefp2(i,4*z-4+m))*Esp2(i,4*z-4+m)*exp(-im*B*dist3(i,4*z-4+m));
                      
                          
                    

                          %%%%%%%%%%%%%%%%%%%angulos
                            m4(i,4*z-4+m) = (yTx(z)-Ry(i,4*z-4+m))/(xTx-Rx(i,4*z-4+m));
                            m4(~isfinite(m4))= 0;
                            m4(isnan(m4))= 0;
                            ANG1(i,4*z-4+m)= abs((m4(i,4*z-4+m)-m3(i,4*z-4+m))/(1 + m4(i,4*z-4+m)*m3(i,4*z-4+m)));
                            ANG1(isnan(ANG1))= 0;
                            TETA1(i,4*z-4+m)= atand(ANG1(i,4*z-4+m));
                            TETA1(isnan(TETA1))= 0;
                            aux3(i,4*z-4+m) = TETA1(i,4*z-4+m)/2;
                     
                            cosTheta2(i,4*z-4+m) = sqrt(1-((((u0*e0)/(uc*ec))^2)*sind(aux3(i,4*z-4+m))*sind(aux3(i,4*z-4+m))));
                            angulos3(i,4*z-4+m) = cosd(aux3(i,4*z-4+m));
                     %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                     
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     
%                      
                   %  Nincn2(i,4*z-4+m) = (Ninc/cosTheta2(i,4*z-4+m))*[(N0/angulos3(i,4*z-4+m))+(Ninc/cosTheta2(i,4*z-4+m))*tanh(Yc*d*cosTheta2(i,4*z-4+m))]/((Ninc/cosTheta2(i,4*z-4+m))+(N0/angulos3(i,4*z-4+m))*tanh(Yc*d*cosTheta2(i,4*z-4+m))); %Componente normal da impedância de entrada
                   %  Nincp2(i,4*z-4+m) = (Ninc*cosTheta2(i,4*z-4+m))*[(N0*angulos3(i,4*z-4+m))+(Ninc*cosTheta2(i,4*z-4+m))*tanh(Yc*d*cosTheta2(i,4*z-4+m))]/[(Ninc*cosTheta2(i,4*z-4+m))+(N0*angulos3(i,4*z-4+m))*tanh(Yc*d*cosTheta2(i,4*z-4+m))]; %Componente paralela da impedância de entrada
                    % coefref30(i,4*z-4+m) =  (Nincn2(i,4*z-4+m) - (N0*angulos3(i,4*z-4+m)))/(Nincp2(i,4*z-4+m) + (N0*angulos3(i,4*z-4+m)));
                     
                     coefref30(i,4*z-4+m) = ([Ninc*cosTheta2(i,4*z-4+m) - N0*angulos3(i,4*z-4+m)]/[Ninc*cosTheta2(i,4*z-4+m) + N0*angulos3(i,4*z-4+m)]);
                     coefref3(i,4*z-4+m)=  abs(coefref30(i,4*z-4+m));
                     
                     %%%% COmpornente perpendicular
                     coefp30(i,4*z-4+m) = ([N0*cosTheta2(i,4*z-4+m) - Ninc*angulos3(i,4*z-4+m)]/[N0*cosTheta2(i,4*z-4+m) + Ninc*angulos3(i,4*z-4+m)]);
                     coefp3(i,4*z-4+m)=  (coefp30(i,4*z-4+m));
                     Eip3(i,4*z-4+m) = coefp3(i,4*z-4+m)*(Eir2(i,4*z-4+m));  
                     
                     Eir3(i,4*z-4+m) = coefref3(i,4*z-4+m)*(Eir2(i,4*z-4+m));
                     Eir3(isnan(Eir3))= 0;
                     dist4(i,4*z-4+m) = distance(Rx(i,4*z-4+m),Ry(i,4*z-4+m),xTx,yTx(z));
                     pIm2(i,4*z-4+m) = distance(Ix(i,m),Iy(i,m),Rx(i,4*z-4+m),Ry(i,4*z-4+m));
                     Esp3(i,4*z-4+m) = sqrt(pIm2(i,4*z-4+m)/(dist4(i,4*z-4+m)+pIm2(i,4*z-4+m)));
                     Esp3(~isfinite(Esp3))= 0;
                     Esp3(isnan(Esp3))= 0;
                     Efinal(i,4*z-4+m) = Eir3(i,4*z-4+m)*Esp3(i,4*z-4+m)*exp(-im*B*dist4(i,4*z-4+m));
                     %Efinal(i,4*z-4+m) = (Eir2(i,4*z-4+m))*(-coefref3(i,4*z-4+m) + coefp3(i,4*z-4+m))*Esp3(i,4*z-4+m)*exp(-im*B*dist4(i,4*z-4+m));
                     Efinal(~isfinite(Efinal))= 0;
                     distTx2(z) = distance(xTi,yTi,xTx,yTx(z));
                    
                     
%                      if (distTx2(z)<= rho || abs(Alfa(i,4*z-4+m)) > 1  || abs(Alfa1(i,4*z-4+m)) > 1  || abs(Beta1(i,4*z-4+m)) > 1 || abs(Beta(i,4*z-4+m)) > 1)
                        else      
                         Efinal(i,4*z-4+m) = 0;
                         
                      end
                    
                      
         end
    end
end
   

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
        H(z) = (Ein(z) + Ein2(z) + E2R(z));
        
    end
 Etotal = 20*log10(H/max(H));
Etotal(~isfinite(Etotal))= 0;    
 figure
plot(dist, Etotal);
xlim([0 1.5]);
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FDTD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
%% Pontos Inicialização
%%%% célula FDTD 
fx = 3;
fy = 1.75;

%Ponto S2
fx1 = 3;
fy1 = 3.25;

%Ponto S3
fx2 = 4.5;
fy2 = 3.25;

%Ponto S4
fx3 = 4.5;
fy3 = 1.75;

%%Criação de segmentos da área de cálculo
fdtd1 = line([fx fx1], [fy fy1],'Color','r','LineStyle','--'); %Segmento 1
fdtd2 = line([fx1 fx2], [fy1 fy2],'Color','r','LineStyle','--'); %Segmento 2
fdtd3 = line([fx2 fx3], [fy2 fy3],'Color','r','LineStyle','--'); %Segmento 3
fdtd4 = line([fx3 fx], [fy3 fy],'Color','r','LineStyle','--'); %Segmento 4
%% Dimensão da celula em X

dx = lambda/20; % Dimensão em x de uma célula unitária.
xmax = 2; % Área de cálculo (dimensão x)
xmin = 0;
%axis celula
zx = xmin : dx :xmax ;
Nx = length(zx); %passo espacial em x


%% dimensao da celula em Y

dy = lambda/20; % Dimensão em y de uma célula unitária.
ymax  = 2; % Área de calculo (dimensao y); 
ymin = 0;
zy = ymin : dy : ymax;
Ny = length(zy); %passo espacial em y
posx = floor(0.5*Nx/xmax);
posy0 = floor(1.75*Ny/ymax);
posy1 = floor(0.25*Ny/ymax);
%% Passo de tempo 
%dt = 0.8/(c0*sqrt((1/(dx^2))+1/(dy^2))); %Passo temporal pelo critério de Courant;
dt = 0.9*dx/(c0*sqrt(2));
STEP = 3000;

%% position of source

 xsource = floor(1*Nx/xmax); % position´[3,2.5]
 ysource = floor(0.5*Ny/ymax);

x1 = 0.5;
y1 = 1;
x2 = 1.5;
y2 = 1;

xmed = floor(1*Nx/xmax); % position´[3,2.5]
ymed = floor(0.5*Ny/ymax);

gf = 0 : dx : 1 ;
Nz = length(gf);

% %Coordenada ponto de recepção
% s0 = distance(x1,y1,x2,y2);
% Npoints = Nx/2;
% points = s0/Npoints;
% xTx = [x1+points:points:x2];
% yTx = 1;

%% Campos iniciais 
Hz = zeros(Nx,Ny);
Ex = zeros(Nx+1,Ny);
Ey = zeros(Nx,Ny+1);
Htotal = zeros(Nz,STEP);

%% Coeficientes
mEz =  (dt)/(e0*dx);
mEz0 = (dt)/(e0*dy);
mHx = (dt)/(u0*dy);
mHy = (dt)/(u0*dx);



A = ((2*e0*er) - (cond*dt))/((2*e0*er) + (cond*dt));
%B = 1/(er*( 1 + ((cond*dt)/(2*e0*er))));
B = (2*dt)/(((2*e0*er) + (cond*dt))*dx);
C = (dt)/(ur*u0*dy);

%% Abrir janela de Figura 
figure('Color','w');
dist = 4; % Área presente do primeiro meio
dist1 = 1;
N0x = dist/dx ;
N0x1 = dist1/dx ;
N0y = dist/dy ;
N0y1 = dist1/dy ;

%% Parametros da PML
pmlx = Nx/3.5; % tamanho da PML ao longo de x
pmly = Ny/3; % tamanho da PML ao longo de y

xn = zeros(1,Nx);
yn = zeros(1,Ny);
fi1 = zeros(1,Nx);
fi2 = ones(1,Nx);
fi3 = ones(1,Nx);
gi2 = ones(1,Nx);
gi3 = ones(1,Nx);
fj1 = zeros(1,Ny);
fj2 = ones(1,Ny);
fj3 = ones(1,Ny);
gj2 = ones(1,Ny);
gj3 = ones(1,Ny);

gaz = ones(Nx,Ny);
gbz = zeros(Nx,Ny);

%% Determinação coeficientes da PML
for i = 1:pmlx
    xnum = pmlx - i ; 
    xn(i)=0.333*((xnum/pmlx)^3 );
    fi1(i)=xn(i); % para dimensão x 
    fi1(Nx-i)=xn(i);  % para segundo lado da pml em x
    fi2(i)=1/(1+xn(i));
    fi2(Nx-i)=1/(1+xn(i));
    fi3(i)=(1-xn(i))/(1+xn(i));
    fi3(Nx-i)=(1-xn(i))/(1+xn(i));
    gi2(i)=1/(1+xn(i));
    gi2(Nx-i)=1/(1+xn(i));
    gi3(i)=(1-xn(i))/(1+xn(i));
    gi3(Nx-i)=(1-xn(i))/(1+xn(i));
    gaz(i)=1/(1+2*xn(i));
    gaz(Nx-i)=1/(1+2*xn(i));

end

for j = 1:pmly
    ynum = pmly-j ; % tornar a condutividade 0 no limite do contorno da pml
    yn(j)=0.333*((ynum/pmly)^3 );
    fj1(j)=yn(j); 
    fj1(Ny-j)=yn(j);
    fj2(j)=1/(1+yn(j));
    fj2(Ny-j)=1/(1+yn(j));
    fj3(j)=(1-yn(j))/(1+yn(j));
    fj3(Ny-j)=(1-yn(j))/(1+yn(j));
    gj2(j)=1/(1+yn(j));
    gj2(Ny-j)=1/(1+yn(j));
    gj3(j)=(1-yn(j))/(1+yn(j));
    gj3(Ny-j)=(1-yn(j))/(1+yn(j));
    gaz(:,j)=1/(1+2*yn(j));
    gaz(:,Ny-j)=1/(1+2*yn(j));
end
%1refl
% H1rl = [-2.67743217772210 + 5.21410836636886i -2.30752765325678 + 5.37541395337246i -1.91733802257967 + 5.51403057120064i -1.50816795620705 + 5.62707192247315i -1.08168127196818 + 5.71166877447814i -0.639917761412393 + 5.76501033959695i -0.185305799288049 + 5.78439009525459i 0.279330101705455 + 5.76725575758040i 0.750767643366925 + 5.71126298683232i 1.22539184967896 + 5.61433225511323i 1.69920775109727 + 5.47470815013275i 2.16786141490690 + 5.29102022554778i 2.62667001119303 + 5.06234434224950i 3.07066149365684 + 4.78826328006110i 3.49462433190895 + 4.46892524056855i 3.89316755579707 + 4.10509871478695i 4.26079116241259 + 3.69822206021341i 4.59196669331684 + 3.25044602717044i 4.88122751503784 + 2.76466740117882i 5.12326803311383 + 2.24455189356773i 5.31305074342832 + 1.69454442373597i 5.44591968032947 + 1.11986500022182i 5.51771846665431 + 0.526488530234396i 5.52491081544661 - 0.0788929261453348i 5.46470098753253 - 0.688926686210186i 5.33515138525043 - 1.29567566317450i 5.13529417374985 - 1.89071478052408i 4.86523358152017 - 2.46524762399021i 4.52623535586571 - 3.01024284236182i 4.12079975170871 - 3.51658922004433i 3.65271442775297 - 3.97526770859511i 3.12708372606506 - 4.37753803489303i 2.55033103126064 - 4.71513681680451i 1.93017125311418 - 4.98048343293593i 1.27555095789452 - 5.16688923385718i 0.596554291661344 - 5.26876507301026i -0.0957264086831759 - 5.28182160317749i -0.789355540444054 - 5.20325635679007i -1.47175392434930 - 5.03192133364532i -2.12993139296062 - 4.76846468509479i -2.75074924724134 - 4.41544013482077i -3.32120834132701 - 3.97737803505779i -3.82875733959799 - 3.46081244112034i -4.26161450808589 - 2.87425930807322i -4.60909529914098 - 2.22814187588167i -4.86193701225390 - 1.53466050962686i -5.01261101459893 - 0.807605686269532i -5.05561243168043 - 0.0621144457927882i -4.98771691863428 + 0.685627581350131i -4.80819413910616 + 1.41873454102829i -4.51896794735195 + 2.12001310055476i -4.12471401742727 + 2.77239156471901i -3.63288680691986 + 3.35937543058468i -3.05366928412190 + 3.86551766889375i -2.39984077417565 + 4.27689046617404i -1.68656056233350 + 4.58154391652825i -0.931067484453189 + 4.76993629536741i -0.152298572008795 + 4.83532015267284i 0.629567179992534 + 4.77406859436560i 1.39363078564829 + 4.58592682681753i 2.11884881997912 + 4.27417535365280i 2.78464752775954 + 3.84569314694892i 3.37155789757319 + 3.31091165401468i 3.86185195971917 + 2.68365360725125i 4.24015896229801 + 1.98085421193231i 4.49403913659442 + 1.22216630122220i 4.61449258108478 + 0.429455348936748i];
% H1rl = [0.00000000000000 + 0.00000000000000i -4.46850808598830 + 10.4094443601913i -3.71290911525067 + 10.6778742863693i -2.92055468883679 + 10.8967779218249i -2.09466677611861 + 11.0605990213034i -1.23919541644627 + 11.1638945179863i -0.358843137299135 + 11.2014232534441i 0.540920415999826 + 11.1682427512139i 1.45385528981574 + 11.0598132169062i 2.37296111319245 + 10.8721076620643i 3.29050166087314 + 10.6017267454051i 4.19804557841003 + 10.2460166089966i 5.08652460466652 + 9.80318766514982i 5.94631041338666 + 9.27243197046657i 6.76731091936729 + 8.65403651608399i 7.53908655379854 + 7.94948947845367i 8.25098660684585 + 7.16157622489170i 8.89230526477764 + 6.29446166548859i 9.45245643739549 + 5.35375540128405i 9.92116588520334 + 4.34655605174307i 10.2886785234264 + 3.28147116626903i 10.5459781133496 + 2.16860924792983i 10.6850158652066 + 1.01954065478352i 10.6989437887530 - 0.152775494547226i 10.5823479583584 - 1.33410079122496i 10.3314762322852 - 2.50906513278054i 9.94445440640458 - 3.66135342869900i 9.42148431831834 - 4.77393149599574i 8.76501707714156 - 5.82931019825371i 7.97989440571174 - 6.80984574234212i 7.07345107356941 - 7.69807281605248i 6.05557159662310 - 8.47706595326511i 4.93869480569423 - 9.13082318328667i 3.73776055931976 - 9.64466469591140i 2.47009380858220 - 10.0056379772136i 1.15522241835647 - 10.2029196914240i -0.185373393311093 - 10.2282035533725i -1.52858043118819 - 10.0760626079479i -2.85003921936012 - 9.74427376231707i -4.12459440675042 - 9.23409215617515i -5.32680301207914 - 8.55046305416834i -6.43149229762233 - 7.70215944567499i -7.41435670646621 - 6.70183447440326i -8.25258100888628 - 5.56597921642463i -8.92547466732796 - 4.31477818913852i -9.41510053938440 - 2.97185729793350i -9.70687949027876 - 1.56392168661215i -9.79015137638875 - 0.120284107057186i -9.65867228073378 + 1.32771210213476i -9.31102791304513 + 2.74736762504590i -8.75094380106961 + 4.10538771609508i -7.98747434868205 + 5.36871317966131i -7.03505505092863 + 6.50540254815642i -5.91340514111462 + 7.48554277804800i -4.64727167596590 + 8.28216277452266i -3.26601048514393 + 8.87212164431590i -1.80300443074122 + 9.23694191721414i -0.294924916520616 + 9.36355721244151i 1.21915160171492 + 9.24494407981491i 2.69875441178889 + 8.88060911359507i 4.10313309645518 + 8.27690498173376i 5.39244675003517 + 7.44715275637858i 6.52899379403754 + 6.41155284838145i 7.47844416275541 + 5.19687286393678i 8.21103252312615 + 3.83590768703932i 8.70266936661727 + 2.36671486546475i 8.93592645886952 + 0.831636707185859i];
H1rl = [0.00000000000000 + 0.00000000000000i -1.15376382662839 + 2.68770697668623i -0.958669011289837 + 2.75701528560032i -0.754083978103524 + 2.81353596123658i -0.540840635984090 + 2.85583438723907i -0.319958880706197 + 2.88250516979847i -0.0926528996440246 + 2.89219504762730i 0.139665050852728 + 2.88362787879020i 0.375383821683462 + 2.85563149341616i 0.612695924839478 + 2.80716612755661i 0.849603875548635 + 2.73735407506638i 1.08393070745345 + 2.64551011277389i 1.31333500559651 + 2.53117217112475i 1.53533074682842 + 2.39413164003055i 1.74731216595448 + 2.23446262028428i 1.94658377789853 + 2.05254935739348i 2.13039558120629 + 1.84911103010670i 2.29598334665842 + 1.62522301358522i 2.44061375751892 + 1.38233370058941i 2.56163401655691 + 1.12227594678387i 2.65652537171416 + 0.847272211867986i 2.72295984016473 + 0.559932500110907i 2.75885923332716 + 0.263244265117198i 2.76245540772331 - 0.0394464630726674i 2.73235049376626 - 0.344463343105093i 2.66757569262521 - 0.647837831587249i 2.56764708687493 - 0.945357390262042i 2.43261679076008 - 1.23262381199511i 2.26311767793286 - 1.50512142118091i 2.06039987585435 - 1.75829461002217i 1.82635721387649 - 1.98763385429755i 1.56354186303253 - 2.18876901744651i 1.27516551563032 - 2.35756840840226i 0.965085626557091 - 2.49024171646797i 0.637775478947262 - 2.58344461692859i 0.298277145830672 - 2.63438253650513i -0.0478632043415880 - 2.64091080158874i -0.394677770222027 - 2.60162817839504i -0.735876962174649 - 2.51596066682266i -1.06496569648031 - 2.38423234254739i -1.37537462362067 - 2.20772006741039i -1.66060417066350 - 1.98868901752889i -1.91437866979899 - 1.73040622056017i -2.13080725404294 - 1.43712965403661i -2.30454764957049 - 1.11407093794083i -2.43096850612695 - 0.767330254813428i -2.50630550729947 - 0.403802843134766i -2.52780621584022 - 0.0310572228963939i -2.49385845931714 + 0.342813790675065i -2.40409706955308 + 0.709367270514147i -2.25948397367598 + 1.06000655027738i -2.06235700871363 + 1.38619578235950i -1.81644340345993 + 1.67968771529234i -1.52683464206095 + 1.93275883444688i -1.19992038708783 + 2.13844523308702i -0.843280281166751 + 2.29077195826412i -0.465533742226595 + 2.38496814768370i -0.0761492860043975 + 2.41766007633642i 0.314783589996267 + 2.38703429718280i 0.696815392824145 + 2.29296341340877i 1.05942440998956 + 2.13708767682640i 1.39232376387977 + 1.92284657347446i 1.68577894878660 + 1.65545582700734i 1.93092597985959 + 1.34182680362563i 2.12007948114900 + 0.990427105966155i 2.24701956829721 + 0.611083150611099i 2.30724629054239 + 0.214727674468374i];
%total
%H1rl = [11.3803442644502 + 8.55545285006876i 11.4224125127457 + 5.28549525548890i 10.8401322706352 + 2.08148275067413i 9.63896299787229 - 0.933036303671895i 7.85048866449493 - 3.63257116526907i 5.53363329828889 - 5.89543859735748i 2.77457787232669 - 7.60998629086175i -0.314781848153859 - 8.68099060682719i -3.59992693956733 - 9.03587838889621i -6.92836117684720 - 8.63038805748497i -10.1357183111186 - 7.45327195207255i -13.0531339352386 - 5.52965380386501i -15.5155811798352 - 2.92269480053405i -17.3707731366276 + 0.266710302513266i -18.4881540181021 + 3.90239046249990i -18.7674433502191 + 7.81757375079093i -18.1461698820571 + 11.8217090288409i -16.6056398689534 + 15.7092153860012i -14.1748313834786 + 19.2697935274047i -10.9317932118900 + 22.2997959546489i -7.00225165947404 + 24.6140364031303i -2.55528608329139 + 26.0573342057102i 2.20388396064662 + 26.5150453860002i 7.04376128390337 + 25.9218365417076i 11.7186481547077 + 24.2680142181246i 15.9823197620379 + 21.6028321841881i 19.6025211276785 + 18.0343585253835i 22.3754017451232 + 13.7256863819163i 24.1389432534084 + 8.88750521717532i 24.7844390609257 + 3.76729907548846i 24.2651557266215 - 1.36431274125821i 22.6014440253629 - 6.22934725863492i 19.8817677528543 - 10.5594040329853i 16.2593700800641 - 14.1127106032200i 11.9445855128735 - 16.6901228047513i 7.19311119782636 - 18.1489451033554i 2.29085268928369 - 18.4135286990972i -2.46376656007620 - 17.4817791683305i -6.77891814880921 - 15.4269528180513i -10.3886352296389 - 12.3944274673629i -13.0708939636798 - 8.59347929335205i -14.6631970781437 - 4.28445865057496i -15.0744711931801 + 0.237892269873472i -14.2923086136620 + 4.66392575926204i -12.3848859410584 + 8.69007786526645i -9.49725712315697 + 12.0398994198655i -5.84212238368281 + 14.4836106483565i -1.68558748486018 + 15.8547166823159i 2.67118160013059 + 16.0623900976424i 6.91316662968730 + 15.0986006020579i 10.7340710608826 + 13.0393343720008i 13.8586144423813 + 10.0396716490276i 16.0628741788892 + 6.32295022946235i 17.1911026344174 + 2.16469962733342i 17.1676739661497 - 2.12755015023566i 16.0031536405079 - 6.23713464738133i 13.7939114580452 - 9.86210124459091i 10.7151851951031 - 12.7382226922964i 7.00800986024712 - 14.6593857833684i 2.96091716149993 - 15.4937261519772i -1.11225875756727 - 15.1941769758903i -4.89679643087009 - 13.8025101460163i -8.10183798850836 - 11.4464435930663i -10.4835241693101 - 8.32993039814296i -11.8646501288793 - 4.71729149310514i -12.1492194055909 - 0.912359926263790i -11.3306508777835 + 2.76577148501728i];

% baixo
%total
%H2rl = [11.3803442644502 + 8.55545285006895i 11.4224125127458 + 5.28549525548906i 10.8401322706352 + 2.08148275067427i 9.63896299787245 - 0.933036303671689i 7.85048866449511 - 3.63257116526891i 5.53363329828903 - 5.89543859735740i 2.77457787232682 - 7.60998629086170i -0.314781848153692 - 8.68099060682716i -3.59992693956713 - 9.03587838889619i -6.92836117684697 - 8.63038805748501i -10.1357183111185 - 7.45327195207259i -13.0531339352386 - 5.52965380386503i -15.5155811798351 - 2.92269480053406i -17.3707731366275 + 0.266710302513203i -18.4881540181020 + 3.90239046249983i -18.7674433502191 + 7.81757375079073i -18.1461698820571 + 11.8217090288408i -16.6056398689534 + 15.7092153860011i -14.1748313834786 + 19.2697935274047i -10.9317932118900 + 22.2997959546489i -7.00225165947419 + 24.6140364031302i -2.55528608329140 + 26.0573342057101i 2.20388396064664 + 26.5150453860001i 7.04376128390332 + 25.9218365417076i 11.7186481547076 + 24.2680142181247i 15.9823197620378 + 21.6028321841881i 19.6025211276785 + 18.0343585253835i 22.3754017451231 + 13.7256863819165i 24.1389432534084 + 8.88750521717540i 24.7844390609257 + 3.76729907548852i 24.2651557266216 - 1.36431274125811i 22.6014440253630 - 6.22934725863484i 19.8817677528543 - 10.5594040329853i 16.2593700800642 - 14.1127106032200i 11.9445855128737 - 16.6901228047512i 7.19311119782634 - 18.1489451033554i 2.29085268928371 - 18.4135286990972i -2.46376656007616 - 17.4817791683306i -6.77891814880918 - 15.4269528180513i -10.3886352296388 - 12.3944274673631i -13.0708939636798 - 8.59347929335213i -14.6631970781436 - 4.28445865057510i -15.0744711931801 + 0.237892269873292i -14.2923086136620 + 4.66392575926197i -12.3848859410585 + 8.69007786526642i -9.49725712315714 + 12.0398994198653i -5.84212238368287 + 14.4836106483565i -1.68558748486018 + 15.8547166823159i 2.67118160013051 + 16.0623900976424i 6.91316662968719 + 15.0986006020579i 10.7340710608825 + 13.0393343720008i 13.8586144423813 + 10.0396716490277i 16.0628741788892 + 6.32295022946236i 17.1911026344174 + 2.16469962733342i 17.1676739661497 - 2.12755015023556i 16.0031536405079 - 6.23713464738130i 13.7939114580453 - 9.86210124459082i 10.7151851951032 - 12.7382226922964i 7.00800986024715 - 14.6593857833684i 2.96091716149991 - 15.4937261519772i -1.11225875756704 - 15.1941769758904i -4.89679643087000 - 13.8025101460163i -8.10183798850828 - 11.4464435930664i -10.4835241693101 - 8.32993039814307i -11.8646501288794 - 4.71729149310512i -12.1492194055909 - 0.912359926263779i -11.3306508777835 + 2.76577148501715i]; 
%1 refl
% H2rl = [-2.67743217772212 + 5.21410836636885i -2.30752765325676 + 5.37541395337247i -1.91733802257966 + 5.51403057120065i -1.50816795620701 + 5.62707192247316i -1.08168127196816 + 5.71166877447814i -0.639917761412374 + 5.76501033959695i -0.185305799288049 + 5.78439009525459i 0.279330101705475 + 5.76725575758040i 0.750767643366944 + 5.71126298683232i 1.22539184967898 + 5.61433225511322i 1.69920775109723 + 5.47470815013276i 2.16786141490688 + 5.29102022554779i 2.62667001119306 + 5.06234434224948i 3.07066149365685 + 4.78826328006109i 3.49462433190898 + 4.46892524056853i 3.89316755579705 + 4.10509871478697i 4.26079116241259 + 3.69822206021341i 4.59196669331684 + 3.25044602717044i 4.88122751503784 + 2.76466740117882i 5.12326803311384 + 2.24455189356771i 5.31305074342832 + 1.69454442373597i 5.44591968032947 + 1.11986500022182i 5.51771846665431 + 0.526488530234377i 5.52491081544661 - 0.0788929261453548i 5.46470098753253 - 0.688926686210166i 5.33515138525042 - 1.29567566317454i 5.13529417374985 - 1.89071478052410i 4.86523358152017 - 2.46524762399021i 4.52623535586573 - 3.01024284236181i 4.12079975170871 - 3.51658922004433i 3.65271442775294 - 3.97526770859513i 3.12708372606506 - 4.37753803489303i 2.55033103126064 - 4.71513681680451i 1.93017125311416 - 4.98048343293594i 1.27555095789452 - 5.16688923385718i 0.596554291661307 - 5.26876507301026i -0.0957264086831759 - 5.28182160317749i -0.789355540444073 - 5.20325635679007i -1.47175392434932 - 5.03192133364532i -2.12993139296063 - 4.76846468509478i -2.75074924724132 - 4.41544013482078i -3.32120834132699 - 3.97737803505780i -3.82875733959797 - 3.46081244112035i -4.26161450808588 - 2.87425930807324i -4.60909529914099 - 2.22814187588165i -4.86193701225390 - 1.53466050962686i -5.01261101459893 - 0.807605686269514i -5.05561243168043 - 0.0621144457927709i -4.98771691863428 + 0.685627581350131i -4.80819413910617 + 1.41873454102828i -4.51896794735197 + 2.12001310055473i -4.12471401742727 + 2.77239156471901i -3.63288680691985 + 3.35937543058469i -3.05366928412189 + 3.86551766889376i -2.39984077417568 + 4.27689046617402i -1.68656056233348 + 4.58154391652825i -0.931067484453172 + 4.76993629536741i -0.152298572008777 + 4.83532015267284i 0.629567179992568 + 4.77406859436559i 1.39363078564831 + 4.58592682681753i 2.11884881997913 + 4.27417535365280i 2.78464752775954 + 3.84569314694892i 3.37155789757321 + 3.31091165401467i 3.86185195971919 + 2.68365360725123i 4.24015896229801 + 1.98085421193229i 4.49403913659442 + 1.22216630122218i 4.61449258108478 + 0.429455348936715i];
% H2rl = [0.00000000000000 + 0.00000000000000i -4.46850808598826 + 10.4094443601913i -3.71290911525063 + 10.6778742863693i -2.92055468883671 + 10.8967779218250i -2.09466677611857 + 11.0605990213034i -1.23919541644623 + 11.1638945179863i -0.358843137299134 + 11.2014232534441i 0.540920415999864 + 11.1682427512139i 1.45385528981578 + 11.0598132169062i 2.37296111319249 + 10.8721076620643i 3.29050166087307 + 10.6017267454051i 4.19804557840999 + 10.2460166089966i 5.08652460466659 + 9.80318766514978i 5.94631041338669 + 9.27243197046655i 6.76731091936736 + 8.65403651608394i 7.53908655379851 + 7.94948947845370i 8.25098660684585 + 7.16157622489170i 8.89230526477764 + 6.29446166548859i 9.45245643739549 + 5.35375540128405i 9.92116588520336 + 4.34655605174304i 10.2886785234264 + 3.28147116626903i 10.5459781133496 + 2.16860924792983i 10.6850158652066 + 1.01954065478348i 10.6989437887530 - 0.152775494547265i 10.5823479583584 - 1.33410079122492i 10.3314762322851 - 2.50906513278061i 9.94445440640457 - 3.66135342869903i 9.42148431831834 - 4.77393149599574i 8.76501707714158 - 5.82931019825368i 7.97989440571174 - 6.80984574234212i 7.07345107356935 - 7.69807281605254i 6.05557159662310 - 8.47706595326511i 4.93869480569423 - 9.13082318328667i 3.73776055931973 - 9.64466469591141i 2.47009380858220 - 10.0056379772136i 1.15522241835640 - 10.2029196914240i -0.185373393311093 - 10.2282035533725i -1.52858043118822 - 10.0760626079479i -2.85003921936015 - 9.74427376231706i -4.12459440675045 - 9.23409215617514i -5.32680301207911 - 8.55046305416836i -6.43149229762230 - 7.70215944567502i -7.41435670646618 - 6.70183447440328i -8.25258100888626 - 5.56597921642466i -8.92547466732797 - 4.31477818913849i -9.41510053938440 - 2.97185729793350i -9.70687949027877 - 1.56392168661212i -9.79015137638876 - 0.120284107057152i -9.65867228073378 + 1.32771210213476i -9.31102791304515 + 2.74736762504587i -8.75094380106964 + 4.10538771609502i -7.98747434868205 + 5.36871317966131i -7.03505505092861 + 6.50540254815645i -5.91340514111460 + 7.48554277804803i -4.64727167596596 + 8.28216277452263i -3.26601048514390 + 8.87212164431591i -1.80300443074119 + 9.23694191721414i -0.294924916520582 + 9.36355721244151i 1.21915160171499 + 9.24494407981490i 2.69875441178892 + 8.88060911359506i 4.10313309645521 + 8.27690498173374i 5.39244675003518 + 7.44715275637858i 6.52899379403756 + 6.41155284838143i 7.47844416275545 + 5.19687286393673i 8.21103252312617 + 3.83590768703929i 8.70266936661728 + 2.36671486546472i 8.93592645886952 + 0.831636707185796i];
% H2r2 = [0.00000000000000 + 0.00000000000000i -1.53915187852381 + 0.351443424978461i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.53023157042088 + 0.386096704428652i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.52022471540133 + 0.421570758331144i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.50906447775834 + 0.457824017342179i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.49668354798198 + 0.494810923626804i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.48301437267739 + 0.532481743266720i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.46798940319173 + 0.570782384330001i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.45154136348489 + 0.609654222123366i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.43360353770094 + 0.649033933248001i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.41411007780678 + 0.688853340182840i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.39299633156229 + 0.729039268219747i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.37019919097182 + 0.769513416672889i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.34565746123835 + 0.810192246381250i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.31931225010093 + 0.850986885613115i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -1.29110737727979 + 0.891803056568104i 0.00000000000000 + 0.00000000000000i]
H2rl = [0.00000000000000 + 0.00000000000000i -1.15376382662839 + 2.68770697668623i -0.958669011289837 + 2.75701528560032i -0.754083978103524 + 2.81353596123658i -0.540840635984090 + 2.85583438723907i -0.319958880706197 + 2.88250516979847i -0.0926528996440246 + 2.89219504762730i 0.139665050852728 + 2.88362787879020i 0.375383821683462 + 2.85563149341616i 0.612695924839478 + 2.80716612755661i 0.849603875548635 + 2.73735407506638i 1.08393070745345 + 2.64551011277389i 1.31333500559651 + 2.53117217112475i 1.53533074682842 + 2.39413164003055i 1.74731216595448 + 2.23446262028428i 1.94658377789853 + 2.05254935739348i 2.13039558120629 + 1.84911103010670i 2.29598334665842 + 1.62522301358522i 2.44061375751892 + 1.38233370058941i 2.56163401655691 + 1.12227594678387i 2.65652537171416 + 0.847272211867986i 2.72295984016473 + 0.559932500110907i 2.75885923332716 + 0.263244265117198i 2.76245540772331 - 0.0394464630726674i 2.73235049376626 - 0.344463343105093i 2.66757569262521 - 0.647837831587249i 2.56764708687493 - 0.945357390262042i 2.43261679076008 - 1.23262381199511i 2.26311767793286 - 1.50512142118091i 2.06039987585435 - 1.75829461002217i 1.82635721387649 - 1.98763385429755i 1.56354186303253 - 2.18876901744651i 1.27516551563032 - 2.35756840840226i 0.965085626557091 - 2.49024171646797i 0.637775478947262 - 2.58344461692859i 0.298277145830672 - 2.63438253650513i -0.0478632043415880 - 2.64091080158874i -0.394677770222027 - 2.60162817839504i -0.735876962174649 - 2.51596066682266i -1.06496569648031 - 2.38423234254739i -1.37537462362067 - 2.20772006741039i -1.66060417066350 - 1.98868901752889i -1.91437866979899 - 1.73040622056017i -2.13080725404294 - 1.43712965403661i -2.30454764957049 - 1.11407093794083i -2.43096850612695 - 0.767330254813428i -2.50630550729947 - 0.403802843134766i -2.52780621584022 - 0.0310572228963939i -2.49385845931714 + 0.342813790675065i -2.40409706955308 + 0.709367270514147i -2.25948397367598 + 1.06000655027738i -2.06235700871363 + 1.38619578235950i -1.81644340345993 + 1.67968771529234i -1.52683464206095 + 1.93275883444688i -1.19992038708783 + 2.13844523308702i -0.843280281166751 + 2.29077195826412i -0.465533742226595 + 2.38496814768370i -0.0761492860043975 + 2.41766007633642i 0.314783589996267 + 2.38703429718280i 0.696815392824145 + 2.29296341340877i 1.05942440998956 + 2.13708767682640i 1.39232376387977 + 1.92284657347446i 1.68577894878660 + 1.65545582700734i 1.93092597985959 + 1.34182680362563i 2.12007948114900 + 0.990427105966155i 2.24701956829721 + 0.611083150611099i 2.30724629054239 + 0.214727674468374i];
H2r2 = [0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.830851759340075 + 0.323325858946346i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.820438318418481 + 0.347524808097901i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.808978189438900 + 0.372099452516140i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.796418714505258 + 0.396995183763290i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.782708296672767 + 0.422152901766491i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.767796757478086 + 0.447508885021370i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.751635717859946 + 0.472994679894661i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.734179002257626 + 0.498537011787959i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.715383065440499 + 0.524057721043868i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.695207441363838 + 0.549473726582840i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.673615213064930 + 0.574697020349937i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.650573502308978 + 0.599634695722950i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i 0.00000000000000 + 0.00000000000000i -0.626053977366226 + 0.624189013084309i];

reflx2 = find(H2r2~=0);
tamanho = length(reflx2);
Hrex = H2r2(reflx2);
for xx = 1 : tamanho
    H2rl(xx) = H2rl(xx) + Hrex(xx);
end
%Incialização dos campos
Dz=zeros(Nx,Nx);
Ez=zeros(Nx,Ny);
Ey=zeros(Nx,Ny);
Ex=zeros(Nx,Ny);
Hz=zeros(Nx,Ny);
Hx1=zeros(Nx,Ny);
Hy1=zeros(Nx,Ny);
IEx = zeros(Nx,Ny);
IEy = zeros(Nx,Ny);
IHx1 = zeros(Nx,Ny);
IHy1 = zeros(Nx,Ny);
Iz = zeros(Nx,Ny);
im = sqrt(-1);

%% Loop Principal  
for T = 1:STEP
   
 % Calculating Hx
for ny = 1: Ny
    for nx = 1:Nx-1
        if (ny <= floor(1.5*Ny/ymax))
        Chx = Hz(nx+1,ny) - Hz(nx,ny);
        IEx(nx,ny) = IEx(nx,ny) + fi1(nx)*Chx; 
        %Hx(nx,ny) = fj3(ny)*Hx(nx,ny)-fj2(ny)*mHx*(Cex+IHx(nx,ny));
        Ex(nx,ny) = fi3(nx)*Ex(nx,ny) + fi2(nx)*mEz*(Chx + IEx(nx,ny)) ;
        else
        Chx = Hz(nx+1,ny) - Hz(nx,ny);
        IEx(nx,ny) = IEx(nx,ny) + fi1(nx)*Chx; 
        %Hx(nx,ny) = fj3(ny)*Hx(nx,ny)-fj2(ny)*mHx*(Cex+IHx(nx,ny));
        Ex(nx,ny) = A*fi3(nx)*Ex(nx,ny) + fi2(nx)*B*(Chx + IEx(nx,ny)) ;
        end
        
       
    end
end

 % Calculating Hy
% Calculating Hy
for nx = 1: Nx
    for ny = 1:Ny-1
        if (ny <= floor(1.5*Ny/ymax))
        Chy = Hz(nx,ny+1)- Hz(nx,ny);
        IEy(nx,ny) = IEy(nx,ny)+ fj1(ny)*Chy;
        %Hy(nx,ny) = fi3(nx)*Hy(nx,ny)+fi2(nx)*mHy*(Cey+IHy(nx,ny));
        Ey(nx,ny) = fj3(ny)*Ey(nx,ny) - fj2(ny)*mEz*(Chy + IEy(nx,ny));
        else
        Chy = Hz(nx,ny+1)- Hz(nx,ny);
        IEy(nx,ny) = IEy(nx,ny)+ fj1(ny)*Chy;
        %Hy(nx,ny) = fi3(nx)*Hy(nx,ny)+fi2(nx)*mHy*(Cey+IHy(nx,ny));
        Ey(nx,ny) = A*fj3(ny)*Ey(nx,ny) - fj2(ny)*B*(Chy + IEy(nx,ny));
        end
    end
end


 % Calculating Dz
for nx = 2: Nx
    for ny = 2:Ny
            %if ((nx <= floor(1.5*Nx/xdim) && ny >= floor(Nx/xdim)) || (nx <= floor(1.75*Nx/xdim) && ny >= floor(1.25*Nx/xdim)) || (nx <= floor(2*Nx/xdim) && ny >= floor(1.5*Nx/xdim)));%(ny >= N0y && nx <= N0x )
                Hz(nx,ny) = gi3(nx)*gj3(ny)*Hz(nx,ny) + gi2(nx)*gj2(ny)*C*((Ex(nx,ny)-Ex(nx-1,ny))-(Ey(nx,ny)-Ey(nx,ny-1))); % can make separte Rx,Ry
            %else
                
       
    end
end

% Calculating Ez
% Ez = gaz.*(Dz-Iz);
% Iz = Iz + gbz.*Ez;

%% creating sine source
Hi = abs(H);
Hi1 = abs(H1rl);
Hi2 = abs(H2rl);

ang = angle(H);
ang1 = angle(H1rl);
ang2 = angle(H2rl);
%% Assigning source
    for aux = 1:N
    Hz(posy0-aux,posx) = Hi(aux)*sin(2*pi*f*T*dt + ang(aux)); % Hard source acts as Metal wall
    end
 H1r = [3:tata:4]
 Nh1r = length(H1r);
    for aux2 = 1:Nh1r
         Hz(posy1,posx+aux2) = Hi2(aux2)*sin(2*pi*f*T*dt + ang2(aux2)); 
         Hz(posy0,posx+aux2) = Hi2(aux2)*sin(2*pi*f*T*dt + ang2(aux2)); 
    end
    
if T >= 2500
    for r = 1:Nz
        hx(r,T) = Hz(xmed,ymed+r);
    end
end 
    % mostrar resultados
        
    if ~mod(T,50)
       % mostrar campos
        imagesc(zx(1:Nx),zy(1:Ny),Hz(1:Nx,1:Ny));
        %mesh(zx,zy,Ez,'linewidth',2); % Tamanho da espessura da linha
        %surf(zx,zy,Ez);
        %contour(zx,zy,Ez);
        %zlim([-3 1]); %Limite dos eixos
        xlabel('zx');
        ylabel('zy');
        %axis([0 Nx 0 Ny -1 1]);
        title(['Field at step ' num2str(T) ' of ' num2str(STEP)]);
        hold on
         line([1.5 1.5],[0 2],'Color','b');   
        plot(xTx,yTx,'--');
        hold off
      
        getframe;
         
    end
  
end

campo = max(hx,[],2);
F = campo';
% rho = 0.15;
% Np = rho/dx;
% F(1:Np) = 0;
Vlr = 20*log10(F/max(F));
Vlr(~isfinite(Vlr))= 0;
plot(gf,Vlr);

time = toc