[W_1, x_1, p_1] = gauss_elim_piv([-.4 2.8 2.3 .8 -.2; .5 .2 -1.3 2.9 .6;2.4 2.1 1.3 -.1 4.1;-2.8 -.4 .2 -1 -3], 'naive')
[W_2, x_2, p_2] = gauss_elim_piv([-.4 2.8 2.3 .8 -.2; .5 .2 -1.3 2.9 .6;2.4 2.1 1.3 -.1 4.1;-2.8 -.4 .2 -1 -3], 'scaled')
[W_3, x_3, p_3] = gauss_elim_piv([-1.6 -3 1.3 .2 1.5;1.6 .5 2.3 .1 -5.8;-.4 2 1.4 2.1 1.9;-1.1 1.1 -.1 2.6 10], 'partial')

function [W, x, p] = gauss_elim_piv(W, strat)

	if nargin < 2
		strat = 'no_piv';
	end
		
	n = size(W, 1);
	p = 1:n;	
	
	% Scaling %
	d = max(abs(W(:,1:end - 1)), [], 2);
	
	for k = 1:n	
		j = piv_strat(W, k, p, d, strat);
		
		p([j k]) = p([k j]);
			
		W(p(k + 1:n), k) = W(p(k + 1:n), k) / W(p(k), k);
		
		W(p(k + 1:n), k + 1:end) = W(p(k + 1:n), k + 1:end) - (W(p(k + 1:n), k) * W(p(k), k + 1:end));
	end
	
	A = W(:,1:end - 1);
	b = W(:, end);
	
	x = back_sub(A, b, p);
	
end

% Choose a strat %
function [j] = piv_strat(M, k, p, d, strat)
	
	n = size(M, 1);
	
	switch strat
		case 'naive'
			j = find(M(p(k:n), k));
			j = j(1) + (k - 1);
		case 'partial'
			[~, j] = max(abs(M(p(k:n), k)));
			j = j + (k - 1);
		case 'scaled'
			[~, j] = max(abs(M(p(k:n), k)) ./ d(k:n));
			j = j + (k - 1);
		otherwise
			j = k;
	end
		
end

% Back substitution algorithm %
function [x] = back_sub (A, b, p)

	n = size(A, 1);
	x = 1:n;

	for k = n:-1:1
		tmp = 0;
		if k < n
			tmp = sum(A(p(k), k + 1:end) .* x(k+1:n));
		end
		x(k) = (b(p(k)) - tmp) / A(p(k), k);
	end

end