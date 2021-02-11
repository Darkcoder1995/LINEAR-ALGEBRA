%source:::GARY DARGUSH LECTURE 21
function [x,theta2] = ashishsa_Final_p2c(t,theta1, L1,L2,h)
    %First we calculate the corresponding values of x and theta2 from
    %function2 
    [x,theta2]=ashishsa_Final_p2b(theta1,L1,L2,h);
    
    for i=1:length(x)-1
        %Finding the velocity from the change in displacement and dividing it with the difference in the time interval
        %using Forward Finite Difference Method
        v(i)=(x(i+1)-x(i))/(t(i+1)-t(i));
    end
    %Here in the above method we will be able to find the change in
    %velocity for all points except the last point so we need to use the
    %Backward Finite Difference Method to find the velocity at the last point
    %Interval.
    k=length(x);
    v(k)=(x(k)-x(k-1))/(t(k)-t(k-1));
    
    for i=1:length(v)-1
        %Finding the acceleration from the change in velocity and dividing it with the difference in the time interval
        %using Finite Forward Finite Difference Method
        a(i)=(v(i+1)-v(i))/(t(i+1)-t(i));
    end
    %Here in the above method we will be able to find the change in
    %acceleration for all points except the last point so we need to use the
    %Backward Finite Difference to find the velocity at the last point
    %Interval.
    k1=length(v);
    a(k1)=(v(k1)-v(k1-1))/(t(k1)-t(k1-1));
   
   
    m=2;
    n=3;
    figure(1)
    
    %plotting t vs theta1
    
    subplot(m,n,1)
    plot(t,theta1)
    title('t vs theta1')
    xlabel('t')
    ylabel('theta1')
    
    %plotting t vs theta2
    
    subplot(m,n,2)
    plot(t,theta2)
    title('t vs theta2')
    xlabel('t')
    ylabel('theta2')
    
    %plotting t vs x
    
    subplot(m,n,3)
    plot(t,x)
    title('t vs x')
    xlabel('t')
    ylabel('x')
    
    %plotting t vs v
    
    subplot(m,n,4)
    plot(t,v)
    title('t vs v')
    xlabel('t')
    ylabel('v')
    
    %plotting t vs a
    
    subplot(m,n,5)
    plot(t,a)
    title('t vs a')
    xlabel('t')
    ylabel('a')
end