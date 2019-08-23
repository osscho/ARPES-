function y1 = saturate1(matrix1,saturation1)
%saturate matrix1 values at the value saturation*max_value
size1 = size(matrix1)
imax = size1(1);
jmax = size1(2);
matrix1 = matrix1-min(min(matrix1));
max1 = max(max(matrix1))
for i = 1:imax
    for j = 1:jmax
        saturation1*max1
       if matrix1(i,j) > saturation1*max1
            matrix1(i,j) = saturation1*max1;
       end
    end
end
max1 = max(max(matrix1))
y1 = matrix1;
end