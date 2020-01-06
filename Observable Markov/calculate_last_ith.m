function [ output ] = calculate_last_ith(data,k)
% calculate transition probabilities for data using last k values of the
% data this is used when the index is bigger than k
% data
output = [];
    per = permn([0,1,2],k+1);
    for i = 1: 3^(k+1)
        count1 = 0;
        count2 = 0;
        for j = 1:size(data,2)
            if(sum(abs(data(:,j)-int32(fliplr(per(i,:))')))==0)
                count1  = count1+1;
            end
        end
       
        if k~=0
            for j = 1:size(data,2)
                if(sum(abs(data(1:k,j)-int32(fliplr(per(i,2:k+1))')))==0)
                    count2  = count2+1;
                end
            end 
        else
%             size(data,2)
           count2 = size(data,2);
        end
        if count2 == 0
            output = [output, 0]; 
        else
            output = [output, count1/count2]; 
        end
         
       
    end
    

end

