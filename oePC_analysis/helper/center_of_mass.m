function [C_theta, C_bin ] = center_of_mass(data,binNum)
%Compute Center of Mass in circular tracks
%   data - binned neuronal activity of one cell, single column or row array
%   data - must be N *1 or 1 * N array
%   binNum - N
%   C_theta - theta position of COM reported as radius values, ranging from 0 to 2*pi 
%   C_theta should be used as the main output result
%   C_bin - for a reference as to which bin number the COM is located at

    if length(data) ~= binNum
        disp('Error! Lengths of neural data and bin number do not match!')
        return
    end

    if size(data,1) == 1 && size(data,2) > 1
         data = data';
    elseif size(data,1) > 1 && size(data,2) > 1
        disp('Error! Neural data must be one-dimensional!')
        return
    end
    %a1 = transpose(linspace(-pi/binNum,pi*(binNum-1)/binNum,binNum))
    %a1 = transpose(linspace(2*pi/binNum,2*pi*(binNum-1)/binNum,binNum));  % generate linear angular bin centers
    a2 = linspace(-pi,pi,binNum+1);
    a1 = linspace(-pi,pi,binNum)';

    [x2,y2] = pol2cart(a1,data);
    pgon = polyshape(x2,y2);
    [xc,yc] = centroid(pgon);
    [C_theta,~] = cart2pol(xc,yc);
    C_bin = discretize(C_theta, a2);

end