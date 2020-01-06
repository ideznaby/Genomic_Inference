function [probs,alpha_] = hmm( N,Ne,haplotype,geneticDist, prefixSeq,mutateMatrix,snpList,alpha )

prefixLen = size(prefixSeq,2);
prefixProb = 1;
 %class maybe it's better to get it from input;

if prefixLen == 1
    m1 = repmat(haplotype(1,1:N),N,1);
    m2 = m1';
    rowIdx = m1+m2;
    rowIdx = rowIdx + 1;
    for columnIdx =1:3
        alpha(:,:,columnIdx) = reshape(mutateMatrix(rowIdx,columnIdx),N,N)./(N*N);  % this is bugos
    end
else
    if geneticDist(prefixLen) == 0.0
        geneticDist(prefixLen) = 1e-10;
    end

       p = exp(-biasCorrection(4*Ne*geneticDist(prefixLen),1)*1.0/N);
       q = (1-p)/N;
       alpha_j = alpha(:,:,prefixSeq(prefixLen-1)+1);
       term_1 = p*p*alpha_j;
       
       sum_vec1 = sum(alpha_j,1);
       sum_vec2 = sum(alpha_j,2)';
       sum_mat1 = repmat(sum_vec1,N,1);
       sum_mat2 = repmat(sum_vec2,N,1);
       sum_mat2 = sum_mat2';
       term_2_3 = p*q*(sum_mat1+sum_mat2);
       prefixProb = sum(sum_vec1);
       %term 4
       term_4 = q*q*prefixProb;
       
       %alpha_{j+1}
       term = term_1+term_2_3+term_4;
       m1 = repmat(haplotype(prefixLen,1:N),N,1);
       m2 = m1';
       rowIdx = m1+m2;
       rowIdx = rowIdx + 1;
       
       for columnIdx = 1:3
           alpha(:,:,columnIdx) = reshape(mutateMatrix(rowIdx,columnIdx),N,N).*term; 
       end
       alpha = alpha./prefixProb;
       
      
       

    
    
    
    
    
    
end
 probs = [sum(sum(alpha(:,:,1)));sum(sum(alpha(:,:,2)));sum(sum(alpha(:,:,3)))];
alpha_ = alpha;

end

