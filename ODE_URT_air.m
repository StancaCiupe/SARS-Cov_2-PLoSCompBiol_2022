function dy = ODE_URT_air(t,y)% w is a vector with the associated 5th deg polynomial coefficients
global beta k delta p c phi1 phi2 w c2 p2 dl
dy = zeros(7,1);
dy(1) = -beta*y(1)*y(4);
dy(2)=beta*y(1)*y(4)-k*y(2);
dy(3)=k*y(2)-delta*y(3);
dy(4)=p*y(3)/(polyval(w,t))-(c+phi1)*y(4);%polyval evaluates a 5th deg polynonial with the given coefficients at time t
dy(5)=p2*y(3)/(polyval(w,t))-(c2+phi2)*y(5);%polyval evaluates a 5th deg polynonial with the given coefficients at time t
dy(6)=phi1*y(4)-(c2+dl)*y(6);
dy(7)=phi2*y(5)-c2*y(7);%dropped phi1*y(4) 