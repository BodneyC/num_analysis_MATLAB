function [L, U, P] = ob_LUP(W, p)

	n = size(W, 1);
	
	LU = W(p, 1:n);
	L = LU;
	U = LU;
	
	P = zeros(n);
	
	for it = 1:n
		L(it, it) = 1;
		L(it, it + 1:end) = 0;
		
		U(it, 1:it - 1) = 0;
		
		P(it, p(it)) = 1;
	end
	
end
