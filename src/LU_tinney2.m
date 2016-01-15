function [L,U,bus_order] = LU_tinney2(A)
n = size(A,1);
L = eye(n);
U = zeros(n);
bus_numbers = 1:n; % For storing the order in which buses are reduced
bus_order(1:n) = 0;

% Perform initial valence ordering:
valence(1:n) = 0;
for l = 1:n
    valence(l) = nnz(A(l,:))-1;
end
[valence,ind] = sort(valence);    % Sort by valence to obtain an index
A = A(ind,ind); % Re-order A matrix
bus_numbers = bus_numbers(ind); % Shuffle bus vector to reflect sorting
bus_order(1) = bus_numbers(1);  % Store first bus to be reduced

% Perform LU factorization with Tinney-2 ordering:
for i = 1:(n-1) % For each node reduced    
    U(i,i:n) =  A(i,i:n);   % Store elements in i-th row to right of i as
                            % i-th row of upper matrix
    for j = i+1:n   % Each element in the i-th column
        A(j,i) = A(j,i)/A(i,i); % Divide each element in column i below i
                                % itself by the element (i,i)
        L(j,i) = A(j,i);    % Store factor in element of lower matrix
        A(j,j) = A(j,j) - (A(j,i)*A(i,j));  % Subtract the product of (j,i)
                                            % found above from diagonal
                                            % element (j,j)
        for k = i+1:n   % Loop over other non-zero off-diagonal
                        % elements in row i.
            if k ~= j
                A(j,k) = A(j,k)-(A(j,i)*A(i,k));    % Adjust other off-
                                                    % diagonal elements
            end
        end
    end
    
    % Now that node i has been reduced, sort (i+1:n,i+1:n) submatrix:
    A_sub = A((i+1):n,(i+1):n);
    % Filter rounding errors before determining valence:
%     nz = A_sub.*A_sub>1e-5;
%     [m,n] = find(nz);
%     A_sub = sparse(m,n,A_sub(nz));
%     A_sub = full(A_sub);
    clear valence
    clear ind
    % Sort submatrix by valence:
    valence(1:n-i) = 0;
    for l = 1:n-i
        valence(l) = nnz(A_sub(l,:));   % If i is 1, we want valence of 2-
                                        % n, and so on.
    end
    [valence,ind] = sort(valence);    % Sort by valence to obtain index
    A_sub = A_sub(ind,ind); % Use index to sort submatrix
    
%    bus_numbers = bus_numbers(i+1:n);   % Remove number of last bus 
                                        % reduced  
    bus_numbers = bus_numbers(ind);   % Use the new index to re-shuffle
                                        % bus numbers
    bus_order(i+1) = bus_numbers(1);   % Store number of next bus
    A((i+1):n,(i+1):n) = A_sub; % Store submatrix back in A for next reduction
end

U(n,n) = A(n,n);    % Store A's last element in U
% A_f = A;  % Diagnostic tool
% Filter potential rounding errors:
nz = U.*U>1e-5;
[i,j] = find(nz);
U = sparse(i,j,U(nz));

end