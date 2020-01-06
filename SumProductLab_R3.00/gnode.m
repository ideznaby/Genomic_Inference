classdef gnode < factor_node
    %gnodes are the nodes corresponding to the g nodes in the paper they are
    %showing the correlation between nodes we implemented observable markov model
    %for the output of these nodes
    properties
        % add your properties here
        markovorder;
        trans;
        numofSNPs;
        isHMM;
        haplotype;
        genDistance;
        knownfamilymembers;
    end
    
    methods
        function s = gnode(id)
            s = s@factor_node(id);
        end
        
        function msg = factor_fun(s, in_msg, from_id, to_id)
            % add your code here
            % ..................
            global gnodeprobs;
            
            N = 146;% Some Parameters necessary for HMM (Recombination model) to find out more please take a look at this paper: Quantifying Genomic Privacy via Inference Attack with High-Order SNV Correlations
            Ne = 11418;
            
            if ~s.isHMM
                %Calculate the parobabilities using Markov Chain algorithm
                msg = [0 0 0];
                index = mod(to_id - 10000000,s.numofSNPs);
                per = permn([1 2 3],s.markovorder+1);
                for i=1:size(per,1)
                    product = 1;
                    for j= 2:size(per,2)
                        previndex = s.numofSNPs - (index + 1)  + (j-1);
                        msgincome = 1/3;
                        if(previndex < 0 || previndex > size(in_msg,1))
                        else
                            msgincome = in_msg{previndex}(per(i,j));
                        end
                        product = product * msgincome;
                    end
                    msg(per(i,1)) = msg(per(i,1)) + product * s.trans(index+1,i);
                end
            else
                %Calculate the probabilities using HMM
                alpha = zeros(N,N,3);
                predict_ = zeros(1,1000);
                 mutate = computeMutate(N);
                 msg = [0 0 0];
                 index = mod(to_id - 10000000,s.numofSNPs) + 1;
                 if index == s.numofSNPs
                     for j= 1:index
                         [probs{j},alpha] = hmm(N,Ne,s.haplotype(:,2:N+1),s.genDistance(:,2),predict_(1:j),mutate,s.genDistance(:,1),alpha);
                         msgindex = s.numofSNPs - j;
                         probs{j} = probs{j}'; 
                         if msgindex <= 0
                         else
                            preprobs{j} = probs{j} .* in_msg{msgindex};
                         end
                         if j~=index
                             if preprobs{j}(1) == preprobs{j}(2) && preprobs{j}(1) == preprobs{j}(3)
                                 [~,maxindx] = max(in_msg{msgindex});
                                 predict_(j) = maxindx - 1;
                             else
                                 predict_(j) = find(preprobs{j} == max(preprobs{j}))-1;
                             end
                         end
                         probs{j} = probs{j} + 0.0001; % add a little bit noise so it does not contradict with other probabilities calculated using mendel's law and MAFs
                         probs{j} = probs{j} / sum(probs{j});
                     end
                     msg = probs{j};
                     gnodeprobs = probs;
                 else
                     msg = gnodeprobs{index};
                 end
            end
            sum_msg = sum(msg);
            if sum_msg == 0
                msg = ones(size(msg));
            end
            msg = (msg/sum(msg));
        end
    end
end
