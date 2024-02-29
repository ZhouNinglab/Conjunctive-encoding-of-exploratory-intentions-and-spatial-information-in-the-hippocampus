function shiftVec = randshift(origVec)
% Randomly permute the n*m vector along the row. 
% Keeping the mean of each row unchanged
%   
[n,m] = size(origVec);
newVec = zeros(n,m);
cellindx = (1:m);
for i=1:n
   temp = origVec(i,:);
   sf = randi([1,m],1);
   cellindx = cellindx + sf;
   for j=1:m
      if cellindx(j) > m
          cellindx(j) = cellindx(j)-m;
      end
   end
   newVec(i,:) = temp(cellindx); 
    
end
shiftVec = newVec;

% outputVec = origVec;
% outputArg2 = inputArg2;
end

