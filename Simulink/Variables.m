% Variables
clear, clc;

%% Subsistema mecánico
b_l=0.1;
g=9.80665;
l_l=0.5;
l_cm=0.25;
m=1;
m_l=0;
J_l=0.0833;
k_l=g*(m*l_cm+m_l*l_l);
r=120;
J_m=1.4e-5;
b_m=1.5e-5;
Jeq=J_m+J_l/(r^2);
beq=b_m+b_l/(r^2);
T_l=0;
T_leq=T_l/r;

%% Subsistema electromagnético
Pp=3;
Lq=5.8e-3;
Ld=6.6e-3;
Lls=0.8e-3;
lambda_m=0.016;
Rs=1.02;
Rs_40=1.02;    % A temperatura 40ºC
Rs_115=1.32;   % A temperatura 115ºC

%% Subsistema térmico
Rs_REF=1.02;    % No encuentro ese valor en la guía
Cts=0.818;
Rts_amb=146.7;
alpha_cu=3.9e-3;
Ts_REF=40;      % Valor del punto 4.c de la guía

%% Estabilidad
%La fórmula es s^2+2*xita*omega_m*s+omega_m^2=0
b_eq=[(1.5e-5)+((0.1-0.03)/(120^2)), (1.5e-5)+((0.1)/(120^2)), (1.5e-5)+((0.1+0.03)/(120^2))];
R_s=[1.02, 1.32];
P_p=3;
J_eq=[1.4e-5+0.0833/(120^2), 1.4e-5+0.0833/(120^2), 1.4e-5+(0.0833+0.375)/(120^2)];
omega_minimo=[];
omega_maximo=[];
xita_minimo=[];
xita_maximo=[];
Gvqs_Rs_40=[];
Gvqs_Rs_115=[];
for i = 1:3
omega_minimo(i)=sqrt((b_eq(i)*R_s(1)+(3/2)*(P_p^2)*(lambda_m^2))/(J_eq(i)*Lq));
xita_minimo(i)=(J_eq(i)*R_s(1)+b_eq(i)*Lq)/(J_eq(i)*Lq*2*omega_minimo(i));
omega_maximo(i)=sqrt((b_eq(i)*R_s(2)+(3/2)*(P_p^2)*(lambda_m^2))/(J_eq(i)*Lq));
xita_maximo(i)=(J_eq(i)*R_s(2)+b_eq(i)*Lq)/(J_eq(i)*Lq*2*omega_maximo(i));
% Funcion de transferencia theta/v_qs
Gvqs_Rs_40=tf([3/2*Pp*lambda_m],[Lq*J_eq(1) (R_s(1)*J_eq(1)+Lq*b_eq(i)) (b_eq(i)*R_s(1)+3/2*Pp^2*lambda_m^2) 0]);
Gvqs_Rs_115=tf([3/2*Pp*lambda_m],[Lq*J_eq(2) (R_s(2)*J_eq(2)+Lq*b_eq(i)) (b_eq(i)*R_s(2)+3/2*Pp^2*lambda_m^2) 0]);
end
omega_minimo;
omega_maximo;
xita_minimo;
xita_maximo;
%Gráficas de polos y ceros a lazo abierto
Gvqs_min_Rs_40=tf([3/2*Pp*lambda_m],[Lq*J_eq(1) (R_s(1)*J_eq(1)+Lq*b_eq(1)) (b_eq(1)*R_s(1)+3/2*Pp^2*lambda_m^2) 0]);
Gvqs_nom_Rs_40=tf([3/2*Pp*lambda_m],[Lq*J_eq(1) (R_s(1)*J_eq(1)+Lq*b_eq(2)) (b_eq(2)*R_s(1)+3/2*Pp^2*lambda_m^2) 0]);
Gvqs_max_Rs_40=tf([3/2*Pp*lambda_m],[Lq*J_eq(1) (R_s(1)*J_eq(1)+Lq*b_eq(3)) (b_eq(3)*R_s(1)+3/2*Pp^2*lambda_m^2) 0]);
Gvqs_min_Rs_115=tf([3/2*Pp*lambda_m],[Lq*J_eq(2) (R_s(2)*J_eq(2)+Lq*b_eq(1)) (b_eq(1)*R_s(2)+3/2*Pp^2*lambda_m^2) 0]);
Gvqs_nom_Rs_115=tf([3/2*Pp*lambda_m],[Lq*J_eq(2) (R_s(2)*J_eq(2)+Lq*b_eq(2)) (b_eq(2)*R_s(2)+3/2*Pp^2*lambda_m^2) 0]);
Gvqs_max_Rs_115=tf([3/2*Pp*lambda_m],[Lq*J_eq(2) (R_s(2)*J_eq(2)+Lq*b_eq(3)) (b_eq(3)*R_s(2)+3/2*Pp^2*lambda_m^2) 0]);
% Obtener los polos
poles_min_Rs_40 = pole(Gvqs_min_Rs_40);
poles_nom_Rs_40 = pole(Gvqs_nom_Rs_40);
poles_max_Rs_40 = pole(Gvqs_max_Rs_40);
poles_min_Rs_115 = pole(Gvqs_min_Rs_115);
poles_nom_Rs_115 = pole(Gvqs_nom_Rs_115);
poles_max_Rs_115 = pole(Gvqs_max_Rs_115);
% Polos y ceros adicionales
zero1 = -175.8621;
zero2 = -227.5862;
%pzplot(Gvqs_min_Rs_40, 'g')
%hold on
%pzplot(Gvqs_nom_Rs_40, 'b')
%hold on
%pzplot(Gvqs_max_Rs_40, 'm')
%plot(real(zero1), imag(zero1), 'ro'); % Cero adicional en rojo
%legend('Min','Nom','Max')
%title('G(s)=Theta(s)/V_{qs}(s) y cero de G(s)=Theta(s)/T_{l}(s), con Rs_{min}')
%grid on
%figure 
%pzplot(Gvqs_min_Rs_115, 'g')
%hold on
%pzplot(Gvqs_nom_Rs_115, 'b')
%hold on
%pzplot(Gvqs_max_Rs_115, 'm')
%plot(real(zero2), imag(zero2), 'ro'); % Cero adicional en rojo
%legend('Min','Nom','Max')
%title('G(s)=Theta(s)/V_{qs}(s) y cero de G(s)=Theta(s)/T_{l}(s), con Rs_{max}')
%grid on

%% Variables para analizar la observabilidad y controlabilidad
syms lambda P_p J_eq R_s L_q b_eq L_d

%% Observabilidad
C1=[1,0,0,0];
A=[0,1,0,0;
   0,-(b_eq/J_eq),-((3*P_p*lambda)/(2*J_eq)),0;
   0,-((lambda*P_p)/J_eq),-(R_s/L_q),0;
   0,0,0,-R_s/L_d];
observ1=[C1;
        C1*A;
        C1*A^2;
        C1*A^3];
disp('El rango medido desde theta sin sensor en corriente es: ')
disp(rank(observ1));
C2=[0,1,0,0];
observ2=[C2;
        C2*A;
        C2*A^2;
        C2*A^3];
disp('El rango medido desde omega sin sensor en corriente es: ')
disp(rank(observ2));
%Midiendo i_ds para el caso de theta como para omega
C3=[1,0,0,0;
    0,0,0,1];
observ3=[C3;
        C3*A;
        C3*A^2;
        C3*A^3];
disp('El rango medido desde theta con sensor en corriente es: ')
disp(rank(observ3));
C4=[1,0,0,0;
    0,1,0,0;
    0,0,0,1];
observ4=[C4;
        C4*A;
        C4*A^2;
        C4*A^3];
disp('El rango medido desde omega con sensor en corriente y en theta es: ');
disp(rank(observ4));

%% Controlabilidad
B=[0;
   0;
   1/L_q;
   0];
control1=[B, A*B, (A^2)*B, (A^3)*B];
disp('El rango medido de la matriz de controlabilidad es: ');
disp(rank(control1));
B1=[0,0;
   0,0;
   1/L_q,0;
   0,1/L_d];
control2=[B1, A*B1, (A^2)*B1, (A^3)*B1];
disp('El rango medido de la matriz de controlabilidad con consigna de control v_ds agregada es: ');
disp(rank(control2));

%% Modulador de torque
pi=5000;%[rad/s]
Rq=Lq*pi;%[Ω]
Rd=Ld*pi;%[Ω]
R0=Lls*pi;%[Ω]
Jeq_nom=1.4e-5+(0.0833+0.375)/(120^2);
beq_nom=1.5e-5+(0.1)/(120^2);
G_i_Rs_115=tf([0 0 1],[0 Lq/Rq +1]);
%figure
%pzmap(G_i_Rs_115)
%pole(G_i_Rs_115)
%title('Polos')
%grid on
%legend('polos en 5000')

%% Controlador
% Metodo Sintonia Serie
n=2.5;
wpos=800;
wint=(1/n)*wpos;
J_eq=[1.4e-5+0.0833/(120^2), 1.4e-5+0.0833/(120^2), 1.4e-5+(0.0833+0.375)/(120^2)];
for i=1:3
ba(i)=J_eq(i)*n*wpos;
ksa(i)=ba(i)*wpos;
ksia(i)=ksa(i)*wint;
end
GthetaPID1=tf([ba(1) ksa(1) ksia(1)],[J_eq(1) ba(1) ksa(1) ksia(1)]);
GthetaPID2=tf([ba(2) ksa(2) ksia(2)],[J_eq(2) ba(2) ksa(2) ksia(2)]);
GthetaPID3=tf([ba(3) ksa(3) ksia(3)],[J_eq(3) ba(3) ksa(3) ksia(3)]);
figure
pzmap(GthetaPID1, 'g');
zero(GthetaPID1);
pole(GthetaPID1);
hold on;
pzmap(GthetaPID2, 'b');
zero(GthetaPID2);
pole(GthetaPID2);
hold on
pzmap(GthetaPID3, 'm')
zero(GthetaPID3);
pole(GthetaPID3);
legend('Min','Nom','Max')
title('Polos y zeros con Metodo de Sintonia Serie')
grid on
figure
pzmap(GthetaPID1, 'r');
zero(GthetaPID1);
pole(GthetaPID1);
hold on;
pzmap(G_i_Rs_115, 'g');
pole(G_i_Rs_115);
hold on;
pzplot(Gvqs_min_Rs_40, 'm');
hold on;
pzplot(Gvqs_nom_Rs_115, 'b');
legend('PID','Modulador','Planta Rs_40','Planta Rs_115')
title('Comparación de polos y zeros del controlador, del modulador y de la planta')
grid on

%% Observador 
Ke_theta = 6400;
Ke_w = 3200^2;

%% Verificación
% Ganancias en observador modificadas
Ke_theta = 3*3200;
Ke_w = 3*3200^2;
Ke_i = 3200^3;

% Degradación de sensores
wn_theta = 2000*3;
wn_i = 6000*3;
taub_Ts = 20;
% Degradación de actuadores
wsv=30000;

ts=0.001;