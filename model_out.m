%model output function for male-female study
function F=model_out(parm,X)

global beta k delta p c phi1 phi2 w xi c2 p2 dl

%fixed parameters
beta=5e-6; k = 4; c=10; c2=1; dl=1;

%parameters to fit
p=exp(parm(1)); delta=exp(parm(2)); phi1=exp(parm(3)); phi2=exp(parm(4)); p2=exp(parm(5)); 

%MALE 2
%weight function coefficients
w=[-9.47180860959682e-05,0.00170936857675351,-0.00744632883134376,0.0138627418347039,1.00079225481647];
xi=1.092739037;

%initial conditions
y01=[1e+7*xi 0 0 1e+5/0.7 0 0 0];%scaling to TCID50 per mL

%solving the ODE on interval 0 to 10
choice=odeset('AbsTol',10^-9,'RelTol',10^-6);
sol = ode15s(@ODE_URT_air, [0 10], y01, choice); 
n_X=length(X(:,1));

%format and log10-scale model solution
F=zeros(1,1);
for i=1:n_X
    if X(i,2)==1 %infectious virus
        val=deval(sol,X(i,1));
        F(i,1)=val(4,:);
    elseif X(i,2)==2 %infectious aerosol
        val=deval(sol,X(i,1));
        F(i,1)=val(6,:);
    elseif X(i,2)==3 %total virus
        val=deval(sol,X(i,1));
        F(i,1)=val(5,:);
    else  %total aerosol
        val=deval(sol,X(i,1));
        F(i,1)=val(7,:);
    end
end
F=log10(F);
end