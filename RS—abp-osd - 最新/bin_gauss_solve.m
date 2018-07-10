function [mat] = bin_gauss_solve(mat)
% solver for system of linear equations in binary arithmetic using Gaussian
% elimination algorithm.
%mat =[1,0,1,0,1,0,0;1,0,1,1,0,1,0;0,1,1,1,0,0,0;1,0,0,1,0,0,1];
    [m n] = size(mat); % read the size of the original matrix A
    for i=1:m
        if mat(i,i)==0
            for j=(i+1):m
                if mat(j,i)==1%swap
                   tmp = mat(i,:);
                   mat(i,:) = mat(j,:);
                   mat(j,:) = tmp;
                   break;
                end
            end
        end
        if mat(i,i)==1
            rowOne = find(mat((i+1):m,i))+i;
            matOne = ones(length(rowOne),1)*mat(i,:);
            mat(rowOne,:) = bitxor(mat(rowOne,:),matOne);
        end
    end
    for i = 1 : n
         j = find(mat(i:m, i), 1); % finds the FIRST 1 in i-th column starting at i
 
         if isempty(j)
             %error('More than one solution.');
             continue;
         else
             j = j + i - 1;    % we need to add i-1 since j starts at i
             temp = mat(j, :);      % swap rows
            mat(j, :) = mat(i, :);
            mat(i, :) = temp;
            % add i-th row to all rows that contain 1 in i-th column
            % starting at j+1 - remember up to j are zeros
             for k = find(mat( (j+1):m, i ))' 
                 mat(j + k, :) = bitxor(mat(j + k, :), mat(i, :));
             end
         end
     end

    %remove all-zero rows
%     mat = mat( sum(mat,2)>0 ,:);
    
%     if any(sum( mat(:,1:n) ,2)==0)  % no solution because matrix A contains
%     	error('No solution.');      % all-zero row, but with nonzero right hand side
%     end
    % calculate final solution
%     x = zeros(n, 1);     % just an initialization
%      for i = n : -1 : 1  % go back from n to 1
%        x(i) = bitand(dot(mat(i, i:n), x(i:n)) + mat(i, n + 1), 1);
%      end
    for i=m:-1:2
        if mat(i,i)==1
            rowOne = find(mat(1:i-1,i));
            matOne = ones(length(rowOne),1)*mat(i,:);
            mat(rowOne,:) = bitxor(mat(rowOne,:),matOne);
        end
    end
end