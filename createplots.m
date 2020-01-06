A = zeros(30,11);
for i=1:30
A(i,5:end) = smooth(allorder1error(i,5:end));
end
for i=1:30
A(i,1:4) = smooth(allorder1error(i,1:4));
end
B = zeros(30,11);
for i=1:30
B(i,5:end) = smooth(allorder1errorpheno(i,5:end));
end
for i=1:30
B(i,1:4) = smooth(allorder1errorpheno(i,1:4));
end
C = zeros(30,11);
for i=1:30
C(i,5:end) = smooth(allorder2errorpheno(i,5:end));
end
for i=1:30
C(i,1:4) = smooth(allorder2errorpheno(i,1:4));
end
D = zeros(30,11);
for i=1:30
D(i,5:end) = smooth(allorder3errorpheno(i,5:end));
end
for i=1:30
D(i,1:4) = smooth(allorder3errorpheno(i,1:4));
end
E = zeros(30,11);
for i=1:30
E(i,5:end) = smooth(allorder4errorpheno(i,5:end));
end
for i=1:30
E(i,1:4) = smooth(allorder4errorpheno(i,1:4));
end
F = zeros(30,11);
for i=1:30
F(i,5:end) = smooth(allHMMerrorpheno(i,5:end));
end
for i=1:30
F(i,1:4) = smooth(allHMMerrorpheno(i,1:4));
end

figure;
plot(mean(A),'b-*')
hold on;
plot(mean(B),'r-+')
hold on;
plot(mean(C),'g-s')
hold on;
plot(mean(D),'c-o')
hold on;
plot(mean(E),'m-^')
hold on;
plot(mean(F),'k-v')
hold on;

xlabels = {'0','GP3', 'GP4', 'P6', 'C7', 'C8', 'C9', 'C10', 'C11', 'GP2', 'GP1'};
set(gca,'XtickLabel',xlabels(1:11))
title('error');

A = zeros(30,11);
for i=1:30
A(i,5:end) = smooth(allorder1ent(i,5:end));
end
for i=1:30
A(i,1:4) = smooth(allorder1ent(i,1:4));
end
B = zeros(30,11);
for i=1:30
B(i,5:end) = smooth(allorder1entpheno(i,5:end));
end
for i=1:30
B(i,1:4) = smooth(allorder1entpheno(i,1:4));
end
C = zeros(30,11);
for i=1:30
C(i,5:end) = smooth(allorder2entpheno(i,5:end));
end
for i=1:30
C(i,1:4) = smooth(allorder2entpheno(i,1:4));
end
D = zeros(30,11);
for i=1:30
D(i,5:end) = smooth(allorder3entpheno(i,5:end));
end
for i=1:30
D(i,1:4) = smooth(allorder3entpheno(i,1:4));
end
E = zeros(30,11);
for i=1:30
E(i,5:end) = smooth(allorder4entpheno(i,5:end));
end
for i=1:30
E(i,1:4) = smooth(allorder4entpheno(i,1:4));
end
F = zeros(30,11);
for i=1:30
F(i,5:end) = smooth(allHMMentpheno(i,5:end));
end
for i=1:30
F(i,1:4) = smooth(allHMMentpheno(i,1:4));
end
figure;
plot(mean(A),'b-*')
hold on;
plot(mean(B),'r-+')
hold on;
plot(mean(C),'g-s')
hold on;
plot(mean(D),'c-o')
hold on;
plot(mean(E),'m-^')
hold on;
plot(mean(F),'k-v')
hold on;

xlabels = {'0','GP3', 'GP4', 'P6', 'C7', 'C8', 'C9', 'C10', 'C11', 'GP2', 'GP1'};
set(gca,'XtickLabel',xlabels(1:11))
title('entropy');
%%

A = zeros(20,11);
for i=1:20
A(i,5:end) = smooth(allorder3errorpheno(i,5:end));
end
for i=1:20
A(i,1:4) = smooth(allorder3errorpheno(i,1:4));
end
B = zeros(20,11);
for i=1:20
B(i,5:end) = smooth(all40MC3errorpheno(i,5:end));
end
for i=1:20
B(i,1:4) = smooth(all40MC3errorpheno(i,1:4));
end
C = zeros(30,11);
for i=1:30
C(i,5:end) = smooth(all30MC3errorpheno(i,5:end));
end
for i=1:30
C(i,1:4) = smooth(all30MC3errorpheno(i,1:4));
end
D = zeros(30,11);
for i=1:30
D(i,5:end) = smooth(all20MC3errorpheno(i,5:end));
end
for i=1:30
D(i,1:4) = smooth(all20MC3errorpheno(i,1:4));
end
E = zeros(20,11);
for i=1:20
E(i,5:end) = smooth(all10MC3errorpheno(i,5:end));
end
for i=1:20
E(i,1:4) = smooth(all10MC3errorpheno(i,1:4));
end

figure;
plot(mean(A),'b-*','LineWidth',2)
hold on;
plot(mean(B),'r-+','LineWidth',2)
hold on;
plot(mean(C),'g-s','LineWidth',2)
hold on;
plot(mean(D),'c-o','LineWidth',2)
hold on;
plot(mean(E),'m-v','LineWidth',2)
hold on;
grid on;
xlabel('Revealed Family Members');
ylabel('Estimation Error');
legend('50% of SNPs revealed','60% of SNPs revealed','70% of SNPs revealed','80% of SNPs revealed','90% of SNPs revealed')

xlabels = {'0','GP3', 'GP4', 'P6', 'C7', 'C8', 'C9', 'C10', 'C11', 'GP2', 'GP1'};
set(gca,'XtickLabel',xlabels(1:11))
%%
A = zeros(20,11);
for i=1:20
A(i,5:end) = smooth(allorder3entpheno(i,5:end));
end
for i=1:20
A(i,1:4) = smooth(allorder3entpheno(i,1:4));
end
B = zeros(20,11);
for i=1:20
B(i,5:end) = smooth(all40MC3entpheno(i,5:end));
end
for i=1:20
B(i,1:4) = smooth(all40MC3entpheno(i,1:4));
end
C = zeros(10,11);
for i=1:10
C(i,5:end) = smooth(all30MC3entpheno(i,5:end));
end
for i=1:10
C(i,1:4) = smooth(all30MC3entpheno(i,1:4));
end
D = zeros(10,11);
for i=1:10
D(i,5:end) = smooth(all20MC3entpheno(i,5:end));
end
for i=1:10
D(i,1:4) = smooth(all20MC3entpheno(i,1:4));
end
E = zeros(20,11);
for i=1:20
E(i,5:end) = smooth(all10MC3entpheno(i,5:end));
end
for i=1:20
E(i,1:4) = smooth(all10MC3entpheno(i,1:4));
end

figure;
plot(mean(A),'b-*','LineWidth',2)
hold on;
plot(mean(B),'r-+','LineWidth',2)
hold on;
plot(mean(C),'g-s','LineWidth',2)
hold on;
plot(mean(D),'c-o','LineWidth',2)
hold on;
plot(mean(E),'m-v','LineWidth',2)
hold on;
grid on;
xlabel('Revealed Family Members');
ylabel('Normalized Entropy');
legend('50% of SNPs revealed','60% of SNPs revealed','70% of SNPs revealed','80% of SNPs revealed','90% of SNPs revealed')

xlabels = {'0','GP3', 'GP4', 'P6', 'C7', 'C8', 'C9', 'C10', 'C11', 'GP2', 'GP1'};
set(gca,'XtickLabel',xlabels(1:11))

%%

A = zeros(20,11);
for i=1:20
A(i,5:end) = smooth(allHMMerrorpheno(i,5:end));
end
for i=1:20
A(i,1:4) = allHMMerrorpheno(i,1:4);
end
B = zeros(20,11);
for i=1:20
B(i,5:end) = smooth(all40hHMMerrorpheno(i,5:end));
end
for i=1:20
B(i,1:4) = all40hHMMerrorpheno(i,1:4);
end
C = zeros(20,11);
for i=1:20
C(i,5:end) = smooth(all30hHMMerrorpheno(i,5:end));
end
for i=1:20
C(i,1:4) = all30hHMMerrorpheno(i,1:4);
end
D = zeros(11,11);
for i=1:20
D(i,5:end) = smooth(all20hHMMerrorpheno(i,5:end));
end
for i=1:20
D(i,1:4) = all20hHMMerrorpheno(i,1:4);
end
E = zeros(20,11);
for i=1:20
E(i,5:end) = smooth(all10hHMMerrorpheno(i,5:end));
end
for i=1:20
E(i,1:4) = all10hHMMerrorpheno(i,1:4);
end

figure;
plot(mean(A),'b-*','LineWidth',2)
hold on;
plot(mean(B),'r-+','LineWidth',2)
hold on;
plot(mean(C),'g-s','LineWidth',2)
hold on;
plot(mean(D),'c-o','LineWidth',2)
hold on;
plot(mean(E),'m-v','LineWidth',2)
hold on;
grid on;
xlabel('Revealed Family Members');
ylabel('Estimation Error');
legend('50% of SNPs revealed','60% of SNPs revealed','70% of SNPs revealed','80% of SNPs revealed','90% of SNPs revealed')

xlabels = {'0','GP3', 'GP4', 'P6', 'C7', 'C8', 'C9', 'C10', 'C11', 'GP2', 'GP1'};
set(gca,'XtickLabel',xlabels(1:11))
%%
A = zeros(20,11);
for i=1:20
A(i,5:end) = smooth(allHMMentpheno(i,5:end));
end
for i=1:20
A(i,1:4) = smooth(allHMMentpheno(i,1:4));
end
B = zeros(20,11);
for i=1:20
B(i,5:end) = smooth(all40hHMMentpheno(i,5:end));
end
for i=1:20
B(i,1:4) = smooth(all40hHMMentpheno(i,1:4));
end
C = zeros(20,11);
for i=1:20
C(i,5:end) = smooth(all30hHMMentpheno(i,5:end));
end
for i=1:20
C(i,1:4) = smooth(all30hHMMentpheno(i,1:4));
end
D = zeros(20,11);
for i=1:20
D(i,5:end) = smooth(all20hHMMerrorpheno(i,5:end));
end
for i=1:20
D(i,1:4) = all20hHMMentpheno(i,1:4);
end
E = zeros(20,11);
for i=1:20
E(i,5:end) = smooth(all10hHMMentpheno(i,5:end));
end
for i=1:20
E(i,1:4) = smooth(all10hHMMentpheno(i,1:4));
end

figure;
plot(mean(A),'b-*','LineWidth',2)
hold on;
plot(mean(B),'r-+','LineWidth',2)
hold on;
plot(mean(C),'g-s','LineWidth',2)
hold on;
plot(mean(D),'c-o','LineWidth',2)
hold on;
plot(mean(E),'m-v','LineWidth',2)
hold on;
grid on;
xlabel('Revealed Family Members');
ylabel('Normalized Entropy');
legend('50% of SNPs revealed','60% of SNPs revealed','70% of SNPs revealed','80% of SNPs revealed','90% of SNPs revealed')

xlabels = {'0','GP3', 'GP4', 'P6', 'C7', 'C8', 'C9', 'C10', 'C11', 'GP2', 'GP1'};
set(gca,'XtickLabel',xlabels(1:11))