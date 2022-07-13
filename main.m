%Main for male-female
%Fit model (ODE_URT_air.m) to data
%Mean parameter estimate output in variable "parm"
%95% confidence intervals output in variable "conf"
%% Input or load Data
%MALE 2 DATA
%AEROSOL DATA
data1=[
1	1.60205999
2	1.25527251
3	0.3];
data1RNA=[1   4.14921911;
2   3.58046924;
3   3.41161971;
4   2.2622175;
5   2.48144263;
10  1.66];
%%NASAL WASH DATA
data=[1 5.115611;
2   4.612784;
3   2.161368;
4    0.7];
dataRNA=[1    6.59703666;
2   6.68727919;
3   6.45757915;
4   7.03008134;
5   6.36657984;
10  2.40959502];
%weight function coefficients
w=[-9.47180860959682e-05,0.00170936857675351,-0.00744632883134376,0.0138627418347039,1.00079225481647];
xi=1.092739037;

%% Organize and scale data to correct units

xdata = data(:,1); %time of infectious URT
xdata_v=[xdata ones(length(xdata),1)]; 
xdata1 = data1(:,1); %time of infectious aerosol
xdata1_v=[xdata1 2*ones(length(xdata1),1)];
xdata2 = dataRNA(:,1); %time of RNA in URT
xdata2_v=[xdata2 3*ones(length(xdata2),1)];
xdata21 = data1RNA(:,1); %time of RNA in aerosol
xdata21_v=[xdata21 4*ones(length(xdata21),1)];
%infectious virus data swab
ydata_inf = 10/0.7*10.^(data(:,2));%10/0.7 scaling to TCID50 per mL %%infectious virus URT
%infectious virus data aerosol
ydata1_inf = 2.5/0.7*10.^(data1(:,2));%2.5/0.7 scaling to TCID50 per mL %%infectious virus aerosol
%non-infectious virus swab
ydata_TOT = 10*10.^(dataRNA(:,2));%scaling to per mL %%RNA URT
%non-infectious virus aerosol
ydata1_TOT = 2.5*10.^(data1RNA(:,2));%scaling to per mL %%RNA aerosol

%MAKE X and Y
X=[xdata_v;xdata1_v;xdata2_v;xdata21_v]; %time matrix
Y=log10([ydata_inf;ydata1_inf;ydata_TOT;ydata1_TOT]); %log data matrix

%% Least-squares optimization
%minimize the least square distance between log10 model output (model_out.m) and
%Y (log10 data) at times X (data time points)
parm0 =[log(50);log(4);log(1e-2);log(1e-3);log(20);]; %initial parameter guess
[parm,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(X,Y,@model_out,parm0);
parm=exp(parm)';
ssq=norm(R)% norm of residuals

%% 95% Confidence interval
conf = nlparci(parm,R,'covar',CovB);
conf=exp(conf)';

%% Plot fits
y0=[1e+7*xi 0 0 1e+5/0.7 0  0  0]; %IC 
choice=odeset('AbsTol',10^-9,'RelTol',10^-6);
[t,yfit] = ode15s(@ODE_URT_air,[0 10],y0, choice);

figure(1);
subplot(2,2,1)
semilogy(xdata,ydata_inf,'ro',t,yfit(:,4),'r');
axis([0 10 1 1e+9])
hold on
semilogy(xdata2,ydata_TOT,'bd',t, yfit(:,5),'b');
axis([0 10 1 1e+9])
hold on
subplot(2,2,2)
semilogy(xdata1,ydata1_inf,'ro',t,yfit(:,6),'r');
axis([0 10 1 1e+5])
hold on
semilogy(xdata21,ydata1_TOT,'bd',t,yfit(:,7),'b');
axis([0 10 1 1e+5])
hold on

