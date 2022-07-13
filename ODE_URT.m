function dy = Covid_hamster_ode_Nat(t,y)% w is a vector with the associated 5th deg polynomial coefficients
global beta k delta p c w c2 p2
dy = zeros(5,1);
dy(1) = -beta*y(1)*y(4);
dy(2)=beta*y(1)*y(4)-k*y(2);
dy(3)=k*y(2)-delta*y(3);
dy(4)=p*y(3)/(polyval(w,t))-c*y(4);%polyval evaluates a 5th deg polynonial with the given coefficients at time t
dy(5)=p2*y(3)/(polyval(w,t))-(c2)*y(5);%polyval evaluates a 5th deg polynonial with the given coefficients at time t


