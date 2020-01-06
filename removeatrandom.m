function [ knowndata ] = removeatrandom( data, percent )
%This function removes given percentage of data by percent at random and replace these random selected data with -1
knowndata = data;
unknowns = randperm(size(data,1) * size(data,2));
unknowns = unknowns(1:ceil(percent * (size(data,1) * size(data,2))));
for i = 1:size(data,1)
    for j = 1:size(data,2)
        if sum(find((i-1) * size(data,2) + j == unknowns)) > 0
            knowndata(i,j) = -1;
        end
    end
end
end

