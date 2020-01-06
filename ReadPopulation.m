function [ out, ib, tt ] = ReadPopulation( populationpath, IDs )
%ReadPopulation Reads the population data and intersects it with family
%data
file = fopen(populationpath);
format = repmat('%d ',1,100);
D = textscan(file,format);
[tt,ia,ib]=intersect(D{1},IDs(:,1));
D = cell2mat(D);
out = D(ia,2:end);
end

