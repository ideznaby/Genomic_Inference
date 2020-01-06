function [ M ] = probe_calculate( k,data )
% calculates the transition probabilities in markov chain for given order
% as k and given data.
M = zeros(size(data,1),3^(k+1));

    for i = 1:size(data,1)
        if i<k+1
            
            temp = probe_calculate(i-1,data(1:i,:));
            
            row = [];
            for l = 0:2
                row = [row,repmat(temp(i,l*size(temp,2)/3+1:(l+1)*size(temp,2)/3),1,3^(k+1)/3^i) ];
            end
            
            M(i,:) = row;
        else
            M(i,:) = calculate_last_ith(data(i-k:i,1:end),k);
        end
        
    end
    
end

