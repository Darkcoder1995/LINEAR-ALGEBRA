function [mu,V] = ashishsa_final_p1(n)
%INSPIRATION:=LECTURE NOTES: MONTE CARLO METHODS
    
    %GENERATE NORMALLY DISTRIBUTED FORCE VECTOR FOR n VALUES USING normrnd
    F=normrnd(1000,20,[1,n]);
    
    %GENERATE NORMALLY DISTRIBUTED ELASTIC MODULUS VECTOR FOR n VALUES USING normrnd
    E=normrnd(210e9,2e9,[1,n]);
    
    %GENERATE NORMALLY DISTRIBUTED DIAMETER OF BEAM VECTOR FOR n VALUES USING normrnd
    d=normrnd(1e-2,2.5e-4,[1,n]);
    
    %CALCULATE THE DEFLECT AND ANGLE FROM THE RANDOMLY GENERATED FORCE
    %VECTOR, ELASTIC MODULUS VECTOR,DIAMETER OF BEAM VECTOR BY PASSING
    %THESE RANDOMLY GENERATED VECTORS INTO eas501_final_beamDeflection
    %WHICH IS AN CONTENT OBSCURED FILE
    [deflect,angle]=eas501_final_beamDeflection(F,E,d);
    
    %CALCULATE THE MEAN OF THE DEFLECT VECTOR FROM THE DEFLECT VECTOR
    deflect_mean=sum(deflect)/n;
    
    %CALCULATE THE MEAN OF THE ANGLE VECTOR FROM THE ANGLE VECTOR
    angle_mean=sum(angle)/n;
    
    %DELCATE mu VECTOR WHICH CONTAINS THE MEAN OF THE DEFLECT VECTOR
    %AND THE ANGLE VECTOR 
    mu=[deflect_mean;angle_mean];
    
    %DECLARE THE COVARIANCE MATRIX AS A ZERO VECTOR
    V=zeros(2,2);
    
    %LOOP FROM 1 TO NUMBER OF SAMPLES AND FILL IN THE VALUES FOR
    %CO-VARIANCE FOR DEFLECT MEAN, ANGLE_MEAN
    for i=1:n
        V(1,1)=V(1,1)+(deflect(i)-deflect_mean)^2;
        V(1,2)=V(1,2)+(deflect(i)-deflect_mean)*(angle(i)-angle_mean);
        V(2,2)=V(2,2)+(angle(i)-angle_mean)^2;
    end
    
    %DIVIDE THE VALUES BY (NUMBER_OF_SAMPLES-1) TO GET UNBIASED MEAN
    V(1,1)=V(1,1)/(n-1);
    
    %DIVIDE THE VALUES BY (NUMBER_OF_SAMPLES-1) TO GET UNBIASED MEAN
    V(1,2)=V(1,2)/(n-1);
    
    %WE SEE THAT THE NON-DIAGONAL ELEMENTS ARE EQUAL
    V(2,1)=V(1,2);
    
    %DIVIDE THE VALUES BY (NUMBER_OF_SAMPLES-1) TO GET UNBIASED MEAN
    V(2,2)=V(2,2)/(n-1);
end

