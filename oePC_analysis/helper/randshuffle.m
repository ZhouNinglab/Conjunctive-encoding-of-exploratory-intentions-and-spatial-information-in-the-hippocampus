function randomVec = randshuffle(origVec)
%Randomly permute the n*m vector along the row. 
%Keeping the mean of each row unchanged
%   
[n,m] = size(origVec);
newVec = zeros(n,m);
for i=1:n
   temp = origVec(i,:);
   idx = randperm(m);
   newVec(i,:) = temp(idx);
    
    
end
randomVec = newVec;
end

