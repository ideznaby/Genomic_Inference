function [ familytree ] = ReadFamilyTree( familytree_path )
%READFAMILYTREE Reads the family tree from the file given
%   INPUT familytree_path: The path to family tree file
%   OUTPUT familytree: A matrix with 3 columns first one shows the father,
%   the second one the mother and third one the child
fid1 = fopen(familytree_path,'r');
i = 1;
while ~feof(fid1)
    line = fgets(fid1); %# read line by line
    if line(1) == '#'
        continue;
    end
    splitted = strsplit(line);
    if strcmp(splitted(end),'') || strcmp(splitted(end), '\n')
        splitted = splitted(1:end-1);
    end
    if length(splitted) == 1
        familytree = zeros(str2num(line),3);
    else
        familytree(i, :) = [str2num(splitted{1}), str2num(splitted{2}), str2num(splitted{3})];
        i = i+1;
    end
end
fclose(fid1);

end

