classdef equ_node < factor_node
    %EQU equ node with no pre-defined number of port
    %   EQU is a subclass of factor node.  It functions as EQU3 but for more
    %   than three ports.
    %
    % Description of factor_fun
    %   input: in_msg, from_id, to_id
    % then msg = in_msg{1}.*in_msg{2}.* .....
    
    properties
        MAF
        knownfamilymembers
        withHMM;
    end
    
    methods
        function s = equ_node(id)
            s = s@factor_node(id);
        end
        
        function msg = factor_fun(s, in_msg, from_id, to_id)
            msg = ones(1,size(in_msg,2));
            for i = 1:size(in_msg,1)
                try
                    msg = msg.*in_msg{i};
                catch
                    disp('error');
                end
            end
            sum_msg = sum(msg);
            if sum_msg == 0
                if ~isempty(msg)
                    disp('Warning: Contradiction in messages received! (Most probabily because there is a contradiction in family data)');
                end
                msg = ones(size(msg));
            end
            msg = msg/sum(msg);
            if ~s.withHMM
                % add noise to MAF so in case it contradicts with the
                % model, it does not make the output msg all zeros
                noisyMAF = s.MAF + 0.0001;
                noisyMAF = noisyMAF / sum(noisyMAF);
            else
                noisyMAF = s.MAF;
            end
            if ~isempty(msg)
                sument = 0;
                sumentMAF = 0;
                for j=1:3
                    ent = log2(msg(j)) * msg(j);
                    entMAF = log2(noisyMAF(j)) * noisyMAF(j);
                    if msg(j) == 0
                        ent = 0;
                        continue;
                    end
                    if noisyMAF(j) == 0
                        entMAF = 0;
                        continue;
                    end
                    sument = sument + ent;
                    sumentMAF = sumentMAF + entMAF;
                end
                entropy = sument / log2(3);
                entropyMAF = sumentMAF / log2(3);
                if entropyMAF > entropy
                    msg = noisyMAF;
                end
            end
        end
    end
end
