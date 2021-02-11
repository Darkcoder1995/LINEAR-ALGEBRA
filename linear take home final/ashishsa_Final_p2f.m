%load Theta1 Matrix
load('theta1.mat')
%Call the function 2e and find the fitted values and then plot these values
%by passing them to function 2c 
f=ashishsa_Final_p2e(t,theta1)
ashishsa_Final_p2c(t,f,0.6,0.8,0.1)
