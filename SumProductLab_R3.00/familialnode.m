classdef familialnode < factor_node
    %these nodes correspond to familial relationship between individuals it is
    %f nodes in the paper.
    properties
        Fr %Mendelian inheritance probabilities the first index
        FrChild %Mendelian inheritance probabilities for child
        MAF %minor allele frequency
        childnodeid %id of the child
        num_of_individuals %num of individuals we considering
        numofSNPs %num of SNPs we have
        knownfamilymembers %list of known families
    end
    
    methods
        function s = familialnode(id)
            s = s@factor_node(id);
        end
        
        function msg = factor_fun(s, in_msg, from_id, to_id)
            % add your code here
            % ..................
            msg = [0 0 0];
            if length(from_id) == 1
                msg = s.MAF;
            else
                for k = 1:3
                    for i = 1:3
                        for j = 1:3
                            if to_id == s.childnodeid
                                msg(k) = msg(k) + (in_msg{1}(i) * in_msg{2}(j) * s.Fr(i,j,k));
                            else
                                if from_id(2) == s.childnodeid
                                        msg(k) = msg(k) + (in_msg{1}(i) * in_msg{2}(j) * (s.FrChild(j,i,k)));
                                else
                                        msg(k) = msg(k) + (in_msg{1}(j) * in_msg{2}(i) * (s.FrChild(j,i,k)));
                                end
                            end
                        end
                    end
                end
            end
            % re-scale msg
            if sum(in_msg{1} == ([1 1 1] / 3)) == 3 && sum(in_msg{2} == ([1 1 1] / 3)) == 3
                msg =[1 1 1];
            end
            sum_msg = sum(msg);
            if sum_msg == 0
                msg = ones(size(msg));
            end
            msg = msg/sum(msg);
        end
    end
end
