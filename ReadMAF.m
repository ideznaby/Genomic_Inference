function [ out, ib ] = ReadMAF( MAFpath, IDs, extracted )
%Reads the MAF data from MAFpath and intersects it with the given IDs
file = fopen(MAFpath);
if ~extracted
    ignore = '%*s ';
    skip = repmat(ignore,1,9);
    s = ['%*s %d ',skip,'%f %f %f'];
    D = textscan(file,s,'HeaderLines',1);
else
    s = '%d %f %f %f';
    D = textscan(file,s,'HeaderLines',0);
end
[tt ia ib]=intersect(D{1},IDs(:,1));
A = D{2};
B = D{3};
C = D{4};
out = [A(ia,:) B(ia,:) C(ia,:)];
end

