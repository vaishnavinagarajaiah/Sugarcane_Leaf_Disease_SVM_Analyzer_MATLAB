function [ Ib ] = binimage(I, tl,tu)
    [m,n]=size(I);
    %remove % sybole to print the floating values of variable I 
    %I    
    for i=1:m
        for j=1:n
            if ((I(i,j)>tl) && (I(i,j)<=tu))
                Ib(i,j)=1;
            else
                Ib(i,j)=0;
            end
            
        end
    end


end