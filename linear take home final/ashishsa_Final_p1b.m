%Integrate the Differential Equation Between the values 0 and 1 with
%initial condition y0=0
[t,y]=ode45(@ashishsa_Final_p1a,[0 1],0);
%Defining the Function
y_original=@(t)(t*exp(-t^2));
%Find the Error Of the Problem at the point with precision of 5 digits
error=round(abs(y(length(y))-y_original(1)),5);