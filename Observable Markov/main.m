fileID = fopen('outputwithdl.txt');
data = textscan(fileID,'%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d');
%%
DATA = cell2mat(data);
% DATA = DATA(1:100,:)
size(DATA)
%%
data0 = probe_calculate(0,DATA);
%%
data1 = probe_calculate(1,DATA);
%%
data2 = probe_calculate(2,DATA);
%%
data3 = probe_calculate(3,DATA);
%%
save(data0,'data0');
