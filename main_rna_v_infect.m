%Main for RNA v infect model
%Fit model (func_RNA_v_infect.m) to data
%Mean parameter estimate output in variable "parm"
%95% confidence intervals output in variable "conf"
%% Enter or load data

%CONTACT DATA
data=[8.54	3.539
11.14	5.289
10.42	3.706
9.36	1.789
7.77	2.456
11.38	6.039
10.83	3.769
10.28	1.789
10.07	4.206
11.31	4.956
10.45	3.289
9.64	1.789
7.93	3.71
10.69	6.071
10.01	4.623
9.24	1.789
8.24	3.21
10.57	5.289
9.84	2.956
9.02	1.789
8.17	4.07
10.85	6.956
10.37	4.841
9.76	1.789
];
%%
%MAKE X and Y
X=10.^data(:,1); %RNA matrix
Y=data(:,2); %infectious virus matrix

%% Least-squares optimization
% %minimize the least square distance between function output (func_RNA_v_infect) and
% %Y (log10 data) at times X (data time points)
parm0 =[0.3; 20;]; %Vm, h, Km 
[parm,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(X,Y,@func_RNA_v_infect,parm0);
ssq=norm(R)% norm of residuals

%%
%95% Confidence interval
conf = nlparci(parm,R,'covar',CovB);

parm=[parm(1) exp(parm(2))];
conf=[conf(1,:)' exp(conf(2,:)')];

%% Plot fits
global Vm h Km
subplot(2,2,1)
f2=@(v) (Vm.*v.^h./(Km.^h+v.^h));
lims=[1 10^12];
fplot(f2,lims) 
hold on
plot(X,Y,'ko')


