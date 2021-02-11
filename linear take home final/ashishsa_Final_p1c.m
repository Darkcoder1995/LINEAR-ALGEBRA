%SOURCE::: GARY DARGUSH LECTURE 21
function [t,y] = ashishsa_Final_p1c(N)
    %Generate N linearly spaced vectors betweeen 0,1 
    t=linspace(0,1,N);
    %take the value of function y at point 1 as 0
    y(1)=0;
    for i=1:N-1
        %find the time difference
        time_diff=t(i+1)-t(i);
        if t(i) == 0
            %for the initial time step initialize the value of expression 1
            exp1=1;
        else
            %find the expression for all the values of time steps
            exp1=ashishsa_Final_p1a(t(i),y(i));
        end
        exp2 = ashishsa_Final_p1a(t(i)+time_diff,y(i)+time_diff*exp1);
        y(i+1) = y(i)+0.5*time_diff*(exp1+exp2);
    end
end