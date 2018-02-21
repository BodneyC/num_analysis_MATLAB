x = forward_sub([1 0 0 0 1;-4 1 0 0 -2;5 4 1 0 14;3 -2 4 1 5])

function [x] = forward_sub(M)

	A = M(:, 1:end - 1);
	b = M(:, end);
	
	n = size(A, 1);
	x = 1:n;

	for k = 1:n
		tmp = 0;
		if k > 1
			tmp = sum(A(k, 1:k - 1) .* x(1:k - 1));
		end
		x(k) = (b(k) - tmp) / A(k, k);
	end

end
