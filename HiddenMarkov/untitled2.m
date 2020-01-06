genDistance = smallgeneticmapchr22combinedb36;
haplotype = smallCEU;
%%
N = 146
Ne = 11418;
mutate = computeMutate(N);
%%
alpha = zeros(N,N,3);
predict_ = zeros(1,1000);
%%
for i = 1:1000
    
    [probs,alpha] = hmm(N,Ne,smallCEU(:,2:N+1),genDistance(:,2),predict_(1:i),mutate,genDistance(:,1),alpha);
    probs
  %  alpha
    if (hidden(i)~=-1)
        predict_(i) = hidden(i);
%         predict_(i)
    else
        predict_(i) = find(probs == max(probs))-1;
    end

end