function yp = ashishsa_Final_p1a(t,y)
    %If T value is 0 return Derivative as 1 Else Return Derivative as the
    %Following Initial Value Condition
    if t==0
        yp = 1;
    else
       yp = y*(-2*t+1/t);
    end
end