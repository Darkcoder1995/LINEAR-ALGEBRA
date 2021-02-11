%SOURCE:::GARY DARGUSH LECTURE 25, DAVID SALAC LECTURE 25
function [fp] = ashishsa_Final_p2e(t,theta1)
    %Initialize the function which is a equation containing t,c
    f=@(t,c) c(1)*t/(t+c(2));
    %Initialize the values of c=(c1,c2)=(5,2)--Random Value. Also we
    %initialize the value of c_previous as (0,0). We want to find the convergance in the values of
    %c so we first initialize Epsilon and then continue until the convergence reaches Epsilon.
    s=0;
    c=[5,2];
    eps=1e-6;
    c_prev=[0,0];
    while(abs(c_prev-c)>eps)
        for i = 1:length(t)
            %We calculate the value of Jacobian for the Matrices -- partial
            %derivatives with respect to c1 and c2
            J(i,1)=t(i)/(t(i)+c(2));
            J(i,2)=-c(1)*t(i)/((t(i)+c(2))^2);
        end
            %We calculate the  Residuals for every data point(R) and Sum of the Squares of the Residuals(s)
        s=0;
        for i = 1:length(t)
            r(i)=theta1(i)-f(t(i),c);
            s=s+(r(i)^2);
        end
        %We Now Calculate The Value Of Gradient and The Hessian Matrix 
         g=2*(J.')*r.';
         H=2*(J.')*J;
         %Initialize the values of c_prev,c_prev as c(1) and c(2)
         c_prev(1)=c(1);
         c_prev(2)=c(2);
         %Now we calculate the changes in the value of c after each
         %iteration.
         del_c=H\g;
         %From the values of del_c we find the convergence in the values of c(1) and c(2)
         %and update the values of c1 and c2 for the next iteration
         c(1) = c(1) + 0.1*del_c(1);
         c(2) = c(2) + 0.1*del_c(2);
    end
    %After getting the converged c value we find the fitted value for each
    %time step.
    for i = 1:length(t)
        fp(i,1) = f(t(i),c);
    end
end