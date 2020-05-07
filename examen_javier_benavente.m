clear; clc;close all

data=readtable('CTG.csv'); %variableNames=data.Properties.VariableNames;
X=table2array(data(:,1:end-1));
y=table2array(data(:,end));
X=X(y==2 | y==3,:); y=y(y==2 | y==3); % y=2: estado fetal sospechoso; y=3: estado fetal patológico
x=X(:,11); % usar solo la variable MLTV (mean value of long term variability)


fetal_sos = (x(y==2));
fetal_pat = (x(y==3));
fetal_sos_mean = mean(fetal_sos);
fetal_pat_mean = mean(fetal_pat);
fetal_sos_std = std(fetal_sos);
fetal_pat_std = std(fetal_pat);

%sacamos histogram de fetal sos y fetal pat para ver su forma 
histogram(fetal_sos,'Normalization','pdf');
hold on 
histogram(fetal_pat,'Normalization','pdf');

%observamos que feta_sos se aproxima a gausiana y fetal_pat a exp
landa = 1/fetal_pat_mean %estimacion maxima verosimilitud
x0= 0:0.001:30;
%MUNDO X
p_x_H0_muestras = landa*exp(-landa*x0) % Izquierda del umbral-fetal_pat
p_x_H1_muestras = normpdf(x0,fetal_sos_mean,fetal_sos_std) % Derecha del umbral-fetal_sos


%MUNDO Y
% Probabilidades de una clase u otra Clase/Total
p_H0 = length(fetal_pat)/length(x)
p_H1 = length(fetal_sos)/length(x)

% ML
syms x1 eta;
p_x_H0 = landa*exp(-landa*x1) % Izquierda del umbral-fetal_pat
p_x_H1 = normpdf(x1,fetal_sos_mean,fetal_sos_std) % Derecha del umbral-fetal_sos
p_x_H0_ML = p_x_H0;
p_x_H1_ML = p_x_H1;
eta = solve(p_x_H0_ML==p_x_H1_ML,x1);
ML = eval(eta)

% Falsa alarma y error ML al haber dos cortes integramos a ML(2)
PFA_ML = eval(int(p_x_H0,ML(1),ML(2)))  % ----- PFa
PM_ML = eval(int(p_x_H1,0,ML(1)))  % ----- PM
PD_ML = 1 - PM_ML                      % ----- PD

% MAP
syms x2 eta2;
p_x_H0_MAP = landa*exp(-landa*x2);
p_x_H1_MAP = normpdf(x2,fetal_sos_mean,fetal_sos_std);
eta2 = solve(p_x_H0_MAP*p_H0==p_x_H1_MAP*p_H1,x2);
MAP = eval(eta2)

% Falsa alarma y error MAP al haber dos cortes integramos a ML(2)
PFA_MAP = eval(int(p_x_H0_MAP,MAP(1),ML(2)))  % ----- PFa
PM_MAP = eval(int(p_x_H1_MAP,0,MAP(1)))    % ----- PM
PD_MAP = 1 - PM_MAP                        % ----- PD


%dibujamos ML
%dibujamos p_x_H0 y p_x_H1 junto con MAP y ML
plot(x0,p_x_H0_muestras)%,'FaceColor',[1 0 1],'FaceAlpha',0.3);
plot(x0,p_x_H1_muestras)%,'FaceColor',[0 1 0],'FaceAlpha',0.3);
plot([ML(1) ML(1)],[0 0.3],'LineWidth',1.5);
plot([MAP(1) MAP(1)],[0 0.3],'LineWidth',1.5);
plot([ML(2) ML(2)],[0 0.3],'LineWidth',1.5);

%calculamos probabilidad de error
Perror_MAP = PM_MAP*p_H1+PFA_MAP*p_H0;
Perror_ML = PM_ML*p_H1+PFA_ML*p_H0;

%explicacion: como podemos comprobar nos sale Perror_MAP = 0.19, Perro_ML=
%0.20, se comprueba que ambas son mas bajas que cuando las etiquetas
%estaban sin separar pensado con un poco de logica es el resultado esperado
%ya que al etiquetar por separado el error es menor, es importante
%aprovechar las etiqeutas ya que estas llevan coste alto de tiempo y dinero y ya que estan es
%bueno aprovecharse de ello.






%sin comprobar
%% ROC ML
umbrales = linspace(0,30);
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


%% ROC MAP
umbrales2 = linspace(0,30);
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
