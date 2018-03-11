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
