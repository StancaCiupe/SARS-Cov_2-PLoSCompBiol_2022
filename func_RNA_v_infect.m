%infectious virus as a function of RNA
function F=func_RNA_v_infect(parm,X)
global Vm h Km
Vm=8.5; h=parm(1); 
Km=exp(parm(2));

F=Vm.*(X).^h./(Km.^h+X.^h);
end

