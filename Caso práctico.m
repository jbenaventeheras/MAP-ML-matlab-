clear;clc; close all;

Frogs = readtable('Frogs_MFCCs.csv');

% Mostramos las clases disponibles
display(categories(categorical(Frogs.Genus)));

% Cogemos los MFCCs y el Genero (En numerico)
MFCCs = table2array(Frogs(:, 1:22));
Genero = double(categorical(Frogs.Genus));

% Sacamos la X y la Y de los géneros Dendropsophus y Leptodactylus
XGenus = MFCCs(Genero==3 | Genero==5,:);
YGenus = Genero(Genero==3 | Genero==5);

% Individualizamos los resultados
Dendropsophus = (XGenus(YGenus==3,:));
Leptodactylus = (XGenus(YGenus==5,:));

% Mostramos en un histograma cada MFCC y los comparamos
figure('Name','Histogramas MFCCs','NumberTitle','off');
for (i = 1:22)
 subplot(6,4,i);
 histogram(Dendropsophus(:,i));
 hold on;
 histogram(Leptodactylus(:,i));
end
legend('Dendropsophus','Leptodactylus');

% Bimodalidad clara en MFCC 5 (Tambien MFCC3)
Den_5 = (XGenus(YGenus==3,5));
Lep_5 = (XGenus(YGenus==5,5));
Den_5_mean = mean(Den_5);
Lep_5_mean = mean(Lep_5);
Den_5_std = std(Den_5);
Lep_5_std = std(Lep_5);

% Caso 5
x0= -1:0.001:1;

p_x_H0 = normpdf(x0,Lep_5_mean,Lep_5_std) % Izquierda del umbral
p_x_H1 = normpdf(x0,Den_5_mean,Den_5_std) % Deracha del umbral


% Probabilidades
p_H0 = length(Lep_5)/length(YGenus)
p_H1 = length(Den_5)/length(YGenus)

% ML
syms x eta;
p_x_H0_ML = normpdf(x,Lep_5_mean,Lep_5_std);
p_x_H1_ML = normpdf(x,Den_5_mean,Den_5_std);
eta = solve(p_x_H0_ML==p_x_H1_ML,x);
ML = eval(eta)
% Falsa alarma y error ML
PFA_ML = eval(int(p_x_H0_ML,ML(1),1))  % ----- PFa
PM_ML = eval(int(p_x_H1_ML,-1,ML(1)))  % ----- PM
PD_ML = 1 - PM_ML                      % ----- PD

% MAP
syms x2 eta2;
p_x_H0_MAP = normpdf(x2,Lep_5_mean,Lep_5_std);
p_x_H1_MAP = normpdf(x2,Den_5_mean,Den_5_std);
eta2 = solve(p_x_H0_MAP*p_H0==p_x_H1_MAP*p_H1,x2);
MAP = eval(eta2)
% Falsa alarma y error MAP
PFA_MAP = eval(int(p_x_H0_MAP,ML(1),1))  % ----- PFa
PM_MAP = eval(int(p_x_H1_MAP,-1,ML(1)))    % ----- PM
PD_MAP = 1 - PM_MAP                        % ----- PD
 

% Dibujamos su densidad espectral de probabilidad normalizada
% entre -1 y 1 (Donde se situan los datos habilitados)
figure('Name','PDF H0 y H1 / Caso 5','NumberTitle','off'); hold on;

area(x0,p_x_H0,'FaceColor',[1 0 1],'FaceAlpha',0.3);
area(x0,p_x_H1,'FaceColor',[0 1 0],'FaceAlpha',0.3);
plot([ML(1) ML(1)],[0 4],'LineWidth',1.5);
plot([MAP(1) MAP(1)],[0 4],'LineWidth',1.5);

xlabel('x'); ylabel('Probability density function MFFC5');
legend('Leptodactylus','Dendropsophus','\eta_{ML}','\eta_{MAP}');


% ROC ML
umbrales = linspace(-1,1);
for iu=1:length(umbrales)
    PFAs(iu)=int(p_x_H0_ML,umbrales(iu),1);
    PMs(iu)=int(p_x_H1_ML,-1,umbrales(iu));
end
PDs=1-PMs;
figure('Name','ROC ML','NumberTitle','off'); hold on; 
plot(PFAs,PDs,'.-');
xlabel('P_{FA}'); ylabel('P_D');
xlim([0 1]); ylim([0 1]);
plot(PFA_ML,PD_ML,'ko');


% ROC MAP
umbrales2 = linspace(-1,1);
for iu2=1:length(umbrales2)
    PFAs2(iu2)=int(p_x_H0_MAP,umbrales2(iu2),1);
    PMs2(iu2)=int(p_x_H1_MAP,-1,umbrales2(iu2));
end
PDs2=1-PMs2;
figure('Name','ROC MAP','NumberTitle','off'); hold on; 
plot(PFAs2,PDs2,'.-');
xlabel('P_{FA}'); ylabel('P_D'); 
xlim([0 1]); ylim([0 1]);
plot(PFA_MAP,PD_MAP,'ro');

%Perror = PM*p_H0+PFA*p_H1;
