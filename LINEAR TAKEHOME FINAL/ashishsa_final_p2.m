function [L,scores] = ashishsa_final_p2(f)
%INSPIRATION:=LECTURE NOTES: PCA 
    %LOAD THE RATINGS FILE
    load(f)
    
    %FIND THE SUM OF THE ELEMENTS COLUMN WISE IE., WE HAVE 9 COLUMNS
    sum_columns=sum(ratings,1);
    
    %FIND SIZE OF RATINGS MATRIX
    [m,n]=size(ratings);
    
    %FIND THE MEAN OF THE ELEMENTS IN EACH COLUMN ROW WISE IE., WE HAVE 250
    %ROWS IN EACH COLUMN SO WE DIVIDE COLUMN SUM BY 249 TO GET UNBIASED
    %ESTIMATE
    mean_columns=(sum_columns/m-1);
    
    %CALCULATE THE COVARIANCE MATRIX. WE USE THE FOLLOWING FUNCTION TO
    %CALCULATE IT
    covariance_matrix=custom_covariance_calculator(mean_columns,ratings,m,n);

    %CALCULATE THE EIGEN VALUES AND THE EIGEN VECTORS USING THE EIGS
    %FUNCTION
    [V1,D]=eigs(covariance_matrix,n);
    
    %THE EIGEN VALUES ARE THE DIAGONAL ELEMENTS OF THE D VARIABLE AND THEY ARE RETURNED AS A VECTOR 
    for i=1:n
        L(i)=D(i,i);   
    end
    
    %CALCULATE THE SCALED MATRIX WHICH IS GIVEN BY THE FOLLOWING FUNCTION
    scaled_matrix=scaled_matrix_calculator(ratings,mean_columns,n);
    
    %CALCULATE THE SCORES FOR EACH ELEMENT ALONG EACH OF THE 9 PRINCIPAL
    %DIRECTIONS BY MULTIPLYING THE EIGEN VECTORS WITH THE SCALED MATRIX
    scores=scaled_matrix*V1;
end

%FUNCTION TO CALCULATE COVARIANCE
function [V]=custom_covariance_calculator(mean_columns,ratings,m,n)
    %INITIALIZE THE COVARIANCE MATRIX
    V=zeros(n,n);
    
    %CALCULATE THE VALUES FOR THE PRINCIPAL DIAGONAL IE., IT FILLS
    %THE VALUES V(1,1),V(2,2),V(3,3),V(4,4),V(5,5),V(6,6),V(7,7),V(8,8),V(9,9)
    for i=1:n
        V(i,i)=V(i,i)+sum((ratings(:,i)-mean_columns(i)).^2);
        V(i,i)=V(i,i)/(m-1);
    end
    
    %CALCULATE THE VALUES FOR THE SECOND PRINCIPAL DIAGONAL IE., IT FILLS
    %THE VALUES V(1,2),V(2,3),V(3,4),V(4,5),V(5,6),V(6,7),V(7,8),V(8,9)
    for i=2:n
        V(i-1,i)=sum((ratings(:,i-1)-mean_columns(i-1)).*(ratings(:,i)-mean_columns(i)));
        V(i-1,i)=V(i-1,i)/(m-1);
    end
    
    %CALCULATE THE VALUES FOR THE THIRD PRINCIPAL DIAGONAL IE., IT FILLS
    %THE VALUES V(1,3),V(2,4),V(3,5),V(4,6),V(5,7),V(6,8),V(7,9)    
    for i=3:n
        V(i-2,i)=sum((ratings(:,i-2)-mean_columns(i-2)).*(ratings(:,i)-mean_columns(i)));
        V(i-2,i)=V(i-2,i)/(m-1);
    end
    
    %CALCULATE THE VALUES FOR THE FOURTH PRINCIPAL DIAGONAL IE., IT FILLS
    %THE VALUES V(1,4),V(2,5),V(3,6),V(4,7),V(5,8),V(6,9)
    for i=4:n
        V(i-3,i)=sum((ratings(:,i-3)-mean_columns(i-3)).*(ratings(:,i)-mean_columns(i)));
        V(i-3,i)=V(i-3,i)/(m-1);
    end
    
    %CALCULATE THE VALUES FOR THE FIFTH PRINCIPAL DIAGONAL IE., IT FILLS
    %THE VALUES V(1,5),V(2,6),V(3,7),V(4,8),V(5,9)
    for i=5:n
        V(i-4,i)=sum((ratings(:,i-4)-mean_columns(i-4)).*(ratings(:,i)-mean_columns(i)));
        V(i-4,i)=V(i-4,i)/(m-1);
    end
    
    %CALCULATE THE VALUES FOR THE SIXTH PRINCIPAL DIAGONAL IE., IT FILLS
    %THE VALUES V(1,6),V(2,7),V(3,8),V(4,9)
    for i=6:n
        V(i-5,i)=sum((ratings(:,i-5)-mean_columns(i-5)).*(ratings(:,i)-mean_columns(i)));
        V(i-5,i)=V(i-5,i)/(m-1);
    end
    
    %CALCULATE THE VALUES FOR THE SEVENTH PRINCIPAL DIAGONAL IE., IT FILLS
    %THE VALUES V(1,7),V(2,8),V(3,9)
    for i=7:n
        V(i-6,i)=sum((ratings(:,i-6)-mean_columns(i-6)).*(ratings(:,i)-mean_columns(i)));
        V(i-6,i)=V(i-6,i)/(m-1);
    end
    
    %CALCULATE THE VALUES FOR THE EIGTH PRINCIPAL DIAGONAL IE., IT FILLS
    %THE VALUES V(1,8),V(2,9)
    for i=8:n
        V(i-7,i)=sum((ratings(:,i-7)-mean_columns(i-7)).*(ratings(:,i)-mean_columns(i)));
        V(i-7,i)=V(i-7,i)/(m-1);
    end
    
    %CALCULATE THE VALUES FOR THE NINTH PRINCIPAL ELEMENT IE., IT FILLS
    %THE VALUES V(1,9)
    V(1,n)=sum((ratings(:,1)-mean_columns(1)).*(ratings(:,n)-mean_columns(n)));
    V(1,n)=V(1,n)/(m-1);
    
    %AT THE END OF THIS OPERATION WE GET A UPPER DIAGONAL COVARIANCE MATRIX
    %FILL IN THE VALUES FOR THE LOWER DIAGONAL ELEMENTS FOR THE COVARIANCE
    %MATRIX BY CALLING THE FOLLOWING FUNCTION AS NON DIAGONAL ELEMENTS
    %FOLLOW A PATTERN IN COVARIANCE MATRIX SO ITS NOT LOGICAL TO
    %RECALCULATE IT. INSTEAD WE REPLACE THE ELEMENTS FROM UPPER DIAGONAL
    %MATRIX TO FILL UP THE LOWER DIAGONAL PART OF THE MATRIX.
    V1=lower_diagonal_fill(V,n);
    
    V=V1;
end

%SUB-FUNCTION TO CALCULATE COVARIANCE
function [V]=lower_diagonal_fill(V,n)
    %FILL IN THE VALUES FOR THE LOWER DIAGONAL ELEMENTS WE OBSERVE THAT IN
    %THE COVARIANCE MATRIX V(1,2)=V(2,1),V(1,3)=V(3,1).... IE., FIRST
    %ROW=FIRST COLUMN
    for i=1:n-1
        V(i+1,1)=V(1,i+1);
    end
    
    %FILL IN THE VALUES FOR THE LOWER DIAGONAL ELEMENTS WE OBSERVE THAT IN
    %THE COVARIANCE MATRIX V(2,3)=V(3,2),V(2,4)=V(4,2).... IE., SECOND
    %ROW=SECOND COLUMN    
    for i=2:n-1
        V(i+1,2)=V(2,i+1);
    end
    
    %FILL IN THE VALUES FOR THE LOWER DIAGONAL ELEMENTS WE OBSERVE THAT IN
    %THE COVARIANCE MATRIX V(3,4)=V(4,3),V(3,5)=V(5,3).... IE., THIRD
    %ROW=THIRD COLUMN    
    for i=3:n-1
        V(i+1,3)=V(3,i+1);
    end
    
    %FILL IN THE VALUES FOR THE LOWER DIAGONAL ELEMENTS WE OBSERVE THAT IN
    %THE COVARIANCE MATRIX V(4,5)=V(5,4),V(4,6)=V(6,4).... IE., FOURTH
    %ROW=FOURTH COLUMN    
    for i=4:n-1
        V(i+1,4)=V(4,i+1);
    end
    
    %FILL IN THE VALUES FOR THE LOWER DIAGONAL ELEMENTS WE OBSERVE THAT IN
    %THE COVARIANCE MATRIX V(5,6)=V(6,5),V(5,7)=V(7,5).... IE., FIFTH
    %ROW=FIFTH COLUMN    
    for i=5:n-1
        V(i+1,5)=V(5,i+1);
    end
    
    %FILL IN THE VALUES FOR THE LOWER DIAGONAL ELEMENTS WE OBSERVE THAT IN
    %THE COVARIANCE MATRIX V(6,7)=V(7,6),V(6,8)=V(8,6).... IE., SIXTH
    %ROW=SIXTH COLUMN   
    for i=6:n-1
        V(i+1,6)=V(6,i+1);
    end
    
    %FILL IN THE VALUES FOR THE LOWER DIAGONAL ELEMENTS WE OBSERVE THAT IN
    %THE COVARIANCE MATRIX V(7,8)=V(8,7),V(7,9)=V(9,7).... IE., SEVENTH
    %ROW=SEVENTH COLUMN   
    for i=7:n-1
        V(i+1,7)=V(7,i+1);
    end
    
    %FILL IN THE VALUES FOR THE LOWER DIAGONAL ELEMENTS WE OBSERVE THAT IN
    %THE COVARIANCE MATRIX V(9,8)=V(8,9)IE., NINTH ROW=NINTH COLUMN    
    V(n,n-1)=V(n-1,n);
end

%FUNCTION TO CALCULATE THE SCALED MATRIX
function [scaled_mat]=scaled_matrix_calculator(ratings,mean_columns,n)
    
    %DEFINE SCALED MATRIX OF THE SAME SIZE AS RATINGS
    scaled_mat=zeros(size(ratings));
    
    %USING FOR LOOP FILL IN THE VALUE OF EACH COLUMN-COLUMN MEAN INTO THE
    %SCALED_MATRIX
    for i=1:n
        scaled_mat(:,i)=ratings(:,i)-mean_columns(i);
    end
end
