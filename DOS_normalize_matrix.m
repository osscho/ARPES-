function y1 = DOSnorm1(matrix1)
%saturate matrix1 values at the value saturation*max_value
size1 = size(matrix1)
imax = size1(1);
jmax = size1(2);
matrix1 = matrix1-min(min(matrix1));
max1 = max(max(matrix1))
    for i = 1:imax
        MCD1(i) = sum(matrix1(i,:));
        matrix1(i,:) = matrix1(i,:)/MCD1(i);
    end
    matrix1 = matrix1-min(min(matrix1));
y1 = matrix1;
end