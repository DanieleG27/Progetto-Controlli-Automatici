
close all; clear all; 

%% parametri fisici del sistema ed equilibrio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=0.2;
J_0=2;
J=[0.5, 0.01, 0.01, 0.08];
psi=[-1, 2, 1.5, -2];
k=20;
teta_e=17/36*pi;

%%% - coppia d'equilibrio
x_e = [teta_e, 0];
u_e = 20*teta_e;

%%% - calcolo momento d'inerzia in x_e (per usarlo nelle matrici del sist. lin.)
J_x1 = momentoDiInerzia(x_e(1),J_0,J,psi);
J_punto_x1 = derivata(x_e(1),J,psi);

%%% - solo per visualizzione, pulsazione minima e massima
omega_plot_min = 1e-2;
omega_plot_max = 1e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNZIONE DI TRASFERIMENTO (definizione matrici e calcolo della FdT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[0 1; (-k*J_x1 + k*x_e(1)*J_punto_x1)/(J_x1*J_x1) -beta/J_x1];
B=[0; 1/J_x1];
C=[1 0];
D=0;

%%% - intervallo di tempo
TT = 0:0.1:10; % da 0 a 10 secondi con passo 0.1

%%% - calcolo la G
 %modello = ss (A ,B ,C , D ) ;
 %G = tf ( modello ) 
 s  = tf('s');
 G = C*inv(s*eye(2) - A)*B + D

%%% - calcolo dei poli

p = pole ( G ) ;
z = zero ( G ) ;

%%% - mostrare poli sul piano complesso
%figure;
%pzmap ( G ) ;
%grid on;
%title('Poli della funzione di trasferiemento');

%%% - mostrare risposta a gradino
%figure;
%YY = step (G , TT ) ;
%plot ( TT , YY ) ;
%title('Risposta a gradino nell%intervallo [0, 10]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specifiche
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ampiezze gradini
WW=3;
DD=3;

% errore a regime
e_star = 0.05;

% attenuazione disturbo sull'uscita
A_d = 50;
omega_d_min = 0.00001;
omega_d_max = 0.1;


% attenuazione disturbo di misura
A_n = 60;
omega_n_min = 1.5e4;
omega_n_max = 1e6;

% Sovraelongazione massima e tempo d'assestamento all'1%
S_star = 8;
T_star = 0.1;

% Margine di fase
Mf_esp = 45;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Diagramma di Bode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
bode(G,{omega_plot_min,omega_plot_max});
title('Diagrammi di Bode G');
grid on; zoom on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Regolatore statico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% valore minimo prescritto per L(0)
mu_s_error = (DD+WW)/e_star;
mu_s_dist  = 10^(A_d/20);

% guadagno minimo del regolatore ottenuto come L(0)/G(0)
G_0 = abs(evalfr(G,j*0)); %guadagno statico della G
G_omega_d_MAX = abs(evalfr(G,j*omega_d_max));

mu_s_RR_err=mu_s_error/G_0;          %regolatore che assicura la specifica sull'errore
mu_s_RR_dist=mu_s_dist/G_omega_d_MAX;%regolatore che assicura la specifica sul disturbo
RR_s =1.25*max(mu_s_RR_err,mu_s_RR_dist); 
%%%%
%CONSIDERAZIONE : nonostante abbiamo scelto RR_S come il massimo fra mu_s_RR_err e mu_s_RR_dist la nostra FdT non rispetta 
%                 per poco il vincolo dato dall' attenuazione del disturbo sull'uscita,
%                 per questo abbiamo deciso di moltiplicare per 1.25 il
%                 regolatore statico in modo da garantire il non
%                 attraversamento della FdT nella zona proibita
%%%%


% Sistema esteso
GG_e = RR_s*G;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diagrammi di Bode di Ge con mappatura specifiche
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mappatura delle specifiche, stampa della FdT insieme alle specifiche

figure(2);
hold on;

% Calcolo specifiche S% => Margine di fase
xi_star = abs(log(S_star*0.9/100))/sqrt(pi^2 + log(S_star*0.9/100)^2);
M_f = max(xi_star*100,Mf_esp);

% Specifiche su d
Bnd_d_x = [omega_d_min; omega_d_max; omega_d_max; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_max; omega_n_max; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 150; 150];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione critica)
omega_Ta_min = 1e-4;
omega_Ta_max = 460/(M_f*T_star);
Bnd_Ta_x = [omega_Ta_min; omega_Ta_max; omega_Ta_max; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}" ;"G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(G,{omega_plot_min,omega_plot_max});
grid on; zoom on;



% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_max;
omega_c_max = omega_n_min;

phi_up = M_f - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
Bnd_Mf_y = [phi_up; phi_up; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);


%specifica riguardo al margine di fase
omega_c_max=omega_n_min; %min pulsaz. dove posso incontrare disturbo di misura

patch_M_f_x = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
patch_M_f_y = [M_f-180; M_f-180; -270; -270];
patch(patch_M_f_x, patch_M_f_y, 'g','FaceAlpha',0.2,'EdgeAlpha',0);

Legend_arg = ["G(j\omega)"; "M_f"]; % Legenda colori
legend(Legend_arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%zone proibite, stampa della G_e
figure(3);
hold on;

patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0); % Specifiche su d
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0); % Specifiche su n
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0); % Specifiche tempo d'assestamento (minima pulsazione critica)

% Legenda colori
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG_e,{omega_plot_min,omega_plot_max});
title('Diagrammi di Bode G estesa');
grid on; zoom on;

% Legenda colori
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0); % Specifiche sovraelongazione (margine di fase)
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

patch(patch_M_f_x, patch_M_f_y, 'g','FaceAlpha',0.2,'EdgeAlpha',0); %specifica riguardo al margine di fase
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Regolatore dinamico con progettazione di rete anticipatrice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mf_star = M_f+10;
omega_c_star = 450;

mag_omega_c_star_dB = abs(evalfr(GG_e,j*omega_c_star));
arg_omega_c_star    = rad2deg(angle(evalfr(GG_e,j*omega_c_star)));

M_star = 1/mag_omega_c_star_dB;
phi_star = Mf_star - 180 - arg_omega_c_star;

tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha_tau = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180));
alpha = alpha_tau / tau;

if min(tau,alpha) < 0
    fprintf('Errore: parametri rete anticipatrice negativi');
    return;
end

R_d_anticipatrice=(1 + tau*s)/(1 + alpha * tau*s);
figure(4);
bode(R_d_anticipatrice);
title('Diagrammi di Bode Rete anticipatrice');
grid on; zoom on;


RR_d = R_d_anticipatrice/(1+s/4000);

RR = RR_s*RR_d
LL = RR*G; % funzione di anello


%return;
figure(5);
hold on;

%%%ripetiamo le specifiche delle zone proibite
% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_arg);


%figure;
%grid on; zoom on;
%step(RR);

%figure;
%grid on; zoom on;
%step(LL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Check prestazioni in anello chiuso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funzione di sensitività complementare
FF = LL/(1+LL);

% Risposta al gradino
figure(6);

T_simulation = 2*T_star;
[y_step,t_step] = step(WW*FF, T_simulation);

plot(t_step,y_step,'b');
title('Grafico della risposta a gradino');
xlabel('t[s]');
ylabel('y(t)');
grid on, zoom on, hold on;

LV = evalfr(WW*FF,0);

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'1%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.01),LV*(1+0.01),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check disturbo in uscita
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funzione di sensitività
SS = 1/(1+LL);
figure(7);

% Simulazione disturbo in uscita sul sistema lineare 
omega_d = 0.01;
tt = 0:1e-2:2e2;
dd = DD*sin(omega_d*tt)+DD*sin(2*omega_d*tt)+DD*sin(3*omega_d*tt);
y_d = lsim(SS,dd,tt);

hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
title('Grafico attenuazione disturbo in uscita');
xlabel('t[s]');
grid on
legend('d(t)','y_d(t)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check disturbo di misura
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funzione di sensitività complementare
FF = LL/(1+LL);
figure(8);

% Simulazione disturbo di misura sul sistema lineare
omega_n = 1e5;
NN      = 3;
tt = 0:1e-5:2*1e-3;
nn = NN*sin(omega_n*tt)+NN*sin(2*omega_n*tt)+NN*sin(3*omega_n*tt);
y_n = lsim(FF,nn,tt);

hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
title('Grafico attenuazione disturbo di misura');
xlabel('t[s]');
grid on
legend('n(t)','y_n(t)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Animazione del sistema controllato

% parametri del disegno (circonferenza)
r=0.5; %raggio
%Wrad = 0:.02:2*pi;
Wrad = 0:.02:2*pi;
Wx = r*cos(Wrad);
Wy = r*sin(Wrad);

for yIndex= 1 : 100
    figure(9);
    clf;
    %bx e by sono le coordinate del punto B che si muove lungo la circonferenza
    theta=y_step(yIndex);
    bx=0.5*cos(theta)+0.5;
    by=-0.5*sin(theta)+1;   %rotazione in senso orario

    %Dx e Dy sono le coordinate del punto D, dipendenti da theta, calcolate
    %tramite un sistema in cui vengono fissate le dimensioni dei segmenti
    %BD=sqrt(1.25) e DC=1
    dx=-(-4*(bx^3)+6*(bx^2)+10*bx-4*bx*(by^2)-4*by*sqrt(-(bx^4)+6*(bx^3)-9*(bx^2)-2*(bx^2)*(by^2)+6*bx*(by^2)+5-(by^4))-6*(by^2)-15)/(2*(4*(bx^2)-12*bx+4*(by^2)+9));
    dy=sqrt(1-(dx-1.5)^2);
    
    axis([0 3.5  -0.5 2.3])
    hold on;
    patch(Wx+0.5, Wy+1,'w','FaceColor','#F3A949');                          % Circonferenza
    LineBD = line([bx,  dx],[by,dy], 'Color','#D7491B','LineWidth',6);      % Segmento 1 della Figura 1 nella traccia
    LineDG = line([dx,  dx+1],[dy,dy], 'Color','#D7491B','LineWidth',6);    % Segmento 2 della Figura 1 nella traccia
    LineCD = line([1.5,  dx],[0,dy], 'Color','#D7491B','LineWidth',6);      % Segmento 3 della Figura 1 nella traccia
    LineGE = line([1.5+1,  dx+1],[0,dy],'Color','#D7491B', 'LineWidth',6);  % Segmento 4 della Figura 1 nella traccia
    
    terra=rectangle('Position',[0,-0.5,3.5,0.5], 'FaceColor','#5DBB63');    % terreno
    plot(0.5,1,'x');                                                        % centro della circonferenza
    title('Interfaccia grafica - SISTEMA CONTROLLATO')
    
    % pausa prima del prossimo frame
    pause(0.0001);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Animazione del sistema non controllato

[y_step_NC,t_step_NC] = step(WW*G, T_simulation*100);

for yIndex= 1 : length(y_step_NC)
    figure(10);
    clf;
    %bx e by sono le coordinate del punto B che si muove lungo la circonferenza
    theta=y_step_NC(yIndex);
    bx=0.5*cos(theta)+0.5;
    by=-0.5*sin(theta)+1;   %rotazione in senso orario

    %Dx e Dy sono le coordinate del punto D, dipendenti da theta, calcolate
    %tramite un sistema in cui vengono fissate le dimensioni dei segmenti
    %BD=sqrt(1.25) e DC=1
    dx=-(-4*(bx^3)+6*(bx^2)+10*bx-4*bx*(by^2)-4*by*sqrt(-(bx^4)+6*(bx^3)-9*(bx^2)-2*(bx^2)*(by^2)+6*bx*(by^2)+5-(by^4))-6*(by^2)-15)/(2*(4*(bx^2)-12*bx+4*(by^2)+9));
    dy=sqrt(1-(dx-1.5)^2);
    
    axis([0 3.5  -0.5 2.3])
    hold on;
    patch(Wx+0.5, Wy+1,'w','FaceColor','#F3A949');                          % Circonferenza
    LineBD = line([bx,  dx],[by,dy], 'Color','#D7491B','LineWidth',6);      % Segmento 1 della Figura 1 nella traccia
    LineDG = line([dx,  dx+1],[dy,dy], 'Color','#D7491B','LineWidth',6);    % Segmento 2 della Figura 1 nella traccia
    LineCD = line([1.5,  dx],[0,dy], 'Color','#D7491B','LineWidth',6);      % Segmento 3 della Figura 1 nella traccia
    LineGE = line([1.5+1,  dx+1],[0,dy],'Color','#D7491B', 'LineWidth',6);  % Segmento 4 della Figura 1 nella traccia
    
    terra=rectangle('Position',[0,-0.5,3.5,0.5], 'FaceColor','#5DBB63');    % terreno
    plot(0.5,1,'x');                                                        % centro della circonferenza
    title('Interfaccia grafica - SISTEMA NON CONTROLLATO')
    
    % pausa prima del prossimo frame
    pause(0.0001);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%