%source:::ALGORITHM FOR NEWTON SOLVER GARY DARGUSH LECTURE 23,24
function [x,theta2]=ashishsa_Final_p2b(theta1,L1,L2,h)
    for i=1:length(theta1)
        %loop throught the values of theta1 and find the corresponding
        %values of Theta2 and x Vector
        [x(i),theta2(i)]=ashishsa_Final_subsection(theta1(i),L1,L2,h);
    end
end

function [Final_X,Final_Theta2] = ashishsa_Final_subsection(theta1,L1,L2,h)
    %Initialize the value of Epsilon.Initialize Random Value of x
    %vector,Theta2(angle) and then find the Initial Value of the function. This gives us Random
    %Initial Value of the postion Of Piston
    eps=1e-6;
    x(1)=1;
    theta2(1)=6;
    f_initial=[x(1)-L1*cos(theta1)-L2*cos(theta2(1));h-L1*sin(theta1)-L2*sin(theta2(1))];
    %run a loop between 1 and 100
    for(i=1:100)
       if(abs(f_initial)>eps)
           %calculate the value of Jacobian Matrix.
            J=[1,L2*sin(theta2(i));0,-L2*cos(theta2(i))];
            %For the iteration one,the first value is f_initial
            if i==1
                F1=f_initial;
            else
                %After plugging in Random Value for the First Value and
                %then calculate the values for each value of Theta in the Theta Matrix
                F1=[x(i)-L1*cos(theta1)-L2*cos(theta2(i));h-L1*sin(theta1)-L2*sin(theta2(i))];
            end
            %According to the Newtons Method J*A=F and solving the system
            %of J,A,F we get A=inverse of J*F. Here A=delta(x) which is the
            %changes in the values of Displacement(x) and Angle(theta2)
            A=-J\F1;
            %update the values of x and theta with each
            %increment in time. This gives us the corresponding values of x
            %and Theta2 for a given value of x and theta1
            x(i+1) = x(i)+A(1);
            theta2(i+1) = theta2(i) + A(2);
       else
           %If the Initial Value of the Function does not satisfy the condition Break Out of the Loop. 
           break
       end 
    end
    %In the End We for each theta1 we are finding final theta2 and x and returning it.
    Final_Theta2 = theta2(length(theta2));
    Final_X=x(length(x));
end