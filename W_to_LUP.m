[L_1, U_1, P_1] = ob_LUP([1/2 1/7 -6/7 6/7 -6/7;4 7 7 6 57;-1/2 3/7 -1/2 -1 -1;1/2 -7/2 -23/2 1 -65/2], [2 4 1 3]')

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
