close all;
clc;
clear all;
tic
%% Pontos Inicialização

c0 = 3*10^8;
fmax = 1*10^9; %Frequencia arbirtrária
e0 = 8.8541878*10^-12;
u0 = 4*pi*10^-7;
er = 7;
ur = 1;
w = 2*pi*fmax;
lambda = c0/fmax;
constB = 2*pi/lambda;
N = sqrt(u0/e0);
cond = 4.73*10^-2;
cmat = 1/(sqrt(er*e0*ur*u0));
Beta = w*sqrt(u0*e0);
%% Dimensão da celula em X

dx = lambda/30; % Dimensão em x de uma célula unitária.
xmax = 5; % Área de cálculo (dimensão x)
xmin = 0;
%axis celula
zx = xmin : dx :xmax ;
Nx = length(zx); %passo espacial em x

%% dimensao da celula em Y

dy = lambda/30; % Dimensão em y de uma célula unitária.
ymax  = 5; % Área de calculo (dimensao y); 
ymin = 0;
zy = ymin : dy : ymax;
Ny = length(zy); %passo espacial em y

%% Passo de tempo 
%dt = 0.8/(c0*sqrt((1/(dx^2))+1/(dy^2))); %Passo temporal pelo critério de Courant;
dt = 0.9*dx/(c0*sqrt(2));
STEP = 2000;

%% position of source
xsource=floor(Nx/2); 
ysource=floor(Ny/2); 
x1 = 2.5;
y1 = 2.5;
x2 = 4;
y2 = 2.5;


gf = 0 : dx : 1.5 ;
Nz = length(gf);


%Coordenada ponto de recepção
s0 = distance(x1,y1,x2,y2);
Npoints = Nx/2;
points = s0/Npoints;
xTx = [x1+points:points:x2];
yTx = 2.5;

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

%% parametros de alimentação
% tau = 100;
% t0 = 10;



%% Source
  tau = 13;
  t0 = 6*tau;
  cc =  32;
  alfa = (4/(cc*dt))^2;

  
%% Main Loop
 im = sqrt(-1);
t = 0;  
for T = 1:STEP
    t = t+1;
   H0 = 325;
   
 Hz(xsource,ysource) = H0*sin(w*dt*T);

  
  for nx = 1:Nx
    for ny = 1:Ny-1
        if nx <= N0x && nx >= N0x1 && ny <= N0y && ny >= N0y1
            Ex(nx,ny) = Ex(nx,ny) + mEz*(Hz(nx,ny+1) - Hz(nx,ny)); 
        else
            Ex(nx,ny) = A*Ex(nx,ny) + B*(Hz(nx,ny+1) - Hz(nx,ny)); 
        end
         end
  end
  
 
    
    for nx = 1:Nx-1
        for ny = 1:Ny
        if nx <= N0x && nx >= N0x1 && ny <= N0y && ny >= N0y1  
            Ey(nx,ny) = Ey(nx,ny) - mEz*(Hz(nx,ny) - Hz(nx+1,ny)); 
        else
            Ey(nx,ny) = A*Ey(nx,ny) - B*(Hz(nx,ny) - Hz(nx+1,ny)); 
        end
        end
    end


    for ny = 2:Ny
        for nx = 2:Nx
            if nx <= N0x && nx >= N0x1 && ny <= N0y && ny >= N0y1
            Hz(nx,ny) = Hz(nx,ny)+ C*(Ey(nx,ny) - Ey(nx-1,ny) + Ex(nx,ny) - Ex(nx,ny-1));
            else
            Hz(nx,ny) = Hz(nx,ny)+ C*(Ey(nx,ny) - Ey(nx-1,ny) + Ex(nx,ny) - Ex(nx,ny-1));    
            end
            
        end
    end
  


 % Ez(xsource,ysource) =  exp(-(((t - t0)^2)/(tau^2)));
if T >= 2500
    for z = 1:Nz
        hx(z,T) = Hz(xsource + z,ysource);
    end
end 
    % mostrar resultados
        
    if ~mod(T,15)
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
         line([4 4],[4 1],'Color','b');
         line([1 4],[4 4],'Color','b');
         line([1 1],[1 4],'Color','b');
         line([4 1],[1 1],'Color','b');     
        plot(xTx,yTx,'--');
        hold off
      
        getframe;
         
    end
  
end

campo = max(hx,[],2);
F = campo';
rho = 0.15;
Np = rho/dx;
F(1:Np) = 0;
Vlr = 20*log10(F/max(F));
Vlr(~isfinite(Vlr))= 0;

plot(gf,Vlr);
time=toc