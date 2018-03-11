% Called function %
function [x] = back_substitute(M)
	
	[rows, cols] = size(M);
		
	M_ut = upper_triangular(M)
	x = back_sub(M_ut);

end

% Convert matrix to upper triangular form				%
%	- This can lead to a significant loss of accuracy	%
function [M] = upper_triangular (M)

	[rows, cols] = size(M);
	reps = min(rows, cols);
	
	for it = 1:reps
		for it_r = 1 + it:rows
			coef = M(it_r, it) / M(it, it);
			M(it_r, it:end) = M(it_r, it:end) - (coef * M(1, it:end)); 
		end
	end
	
end

% Back substitution algorithm %
function [x] = back_sub(M)

	A = M(:, 1:end - 1);
	[rows, cols] = size(A);
	reps = min(rows, cols);
	x(1:cols) = 0;
	
	b = M(:, end)';
		
	cnt = 1;
	
	for it_r = reps:-1:1
		tmp = 0;
		if cnt > 1
			tmp = sum(A(it_r, it_r + 1:end) .* flip(x(1:cnt - 1)));
		end
			
		x(cnt) = (b(it_r) - tmp) / A(it_r, it_r);
		cnt = cnt + 1;		
	end
	
	x = flip(x)';

end
