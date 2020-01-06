function [ mutate ] = computeMutate( N )

theta = sum(1./[1:N-1])^-1;
r = 0.0179940463182*0.5/(0.0179940463182+N);
mutate =[(1-r)^2 , 2*r*(1-r) , r^2; r*(1-r), r^2+(1-r)^2, r*(1-r);r^2, 2*r*(1-r), (1-r)^2];
end

