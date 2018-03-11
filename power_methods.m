clc; clear;

[four_mu_max, four_z_max, four_mu_n, four_k] = power_method    ([4 -1 2;-1 3 -2;2 -2 1], [0.7 -0.5 0.45]', 50)
[four_mu_min, four_z_min, four_mu_n, four_k] = inv_power_method([4 -1 2;-1 3 -2;2 -2 1], [0.7 -0.5 0.45]', 50, 'naive')

[five_mu_max, five_z_min, five_mu_n, five_k] = power_method    ([4 -1 1;-1 3 -2;1 -2 3], [0.7 -0.5 0.45]', 50)
[five_mu_min, five_z_min, five_mu_n, five_k] = inv_power_method([4 -1 1;-1 3 -2;1 -2 3], [0.7 -0.5 0.45]', 50, 'naive')

function [mu, zk, mu_n, it] = power_method(A, z0, its, tol)

	zk = z0;
	mu_n = zeros(its, 1);
	it = 1;
	mu_p = inf;
	mu = 1;
	
	if nargin == 3
		tol = 0.00001;
	end
	
	while it <= its && abs(mu_p - mu) > tol
		mu_p = mu;
		
		y = A * zk;
		[~, mu_idx] = max(abs(y));
		mu = y(mu_idx);
		zk = y / mu;
		
		mu_n(it) = abs(mu - mu_p) / abs(mu);
		it = it + 1;
	end
	
	mu_n(it:end) = [];
	it = it - 1;
	
end

function [lambda, z, mu_n, it] = inv_power_method(A, z, its, strat, tol)

	if nargin == 3
		strat = 'no_sc';
	end
	
	if nargin == 4
		tol = 0.00001;
	end
	
	[W, p] = gauss_elim_piv([A z], strat);
	[L, U, p] = ob_LUP(W, p);
	mu_n = zeros(its);
	mu_p = Inf;
	mu = 1;
	it = 1;

	while it <= its && abs(mu_p - mu) > tol
		mu_p = mu;
		y_bar = forward_sub(L, z);
		y = back_sub(U, y_bar);
		[~, mu_idx] = max(abs(y));
		mu = y(mu_idx);
		z = (1 / mu) * y';
		mu_n(it) = abs(mu - mu_p) / abs(mu);
		it = it + 1;
	end

	lambda = 1 / mu;
	mu_n(it:end) = [];
	mu_n = mu_n';
	it = it - 1;

end

function [x] = back_sub (A, b)

	n = size(A, 1);
	x(1:n) = 0;

	for k = n:-1:1
		tmp = 0;
		if k < n
			tmp = sum(A(k, k + 1:end) .* x(k+1:n));
		end
		x(k) = (b(k) - tmp) / A(k, k);
	end

end

function [x] = forward_sub (A, b)

	n = size(A, 1);
	x(1:n) = 0;

	for k = 1:n
		tmp = 0;
		if k > 1
			tmp = sum(A(k, 1:k - 1) .* x(1:k - 1));
		end
		x(k) = (b(k) - tmp) / A(k, k);
	end

end

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

function [W, p] = gauss_elim_piv(W, strat)

	if nargin < 2
		strat = 'no_sc';
	end
	n = size(W, 1);
	p = 1:n;

	d = max(abs(W), [], 2);
	for k = 1:n
		j = piv_strat(W, k, p, d, strat);
		p([j k]) = p([k j]);
		W(p(k + 1:n), k) = W(p(k + 1:n), k) / W(p(k), k);
		W(p(k + 1:n), k + 1:end) = W(p(k + 1:n), k + 1:end) - (W(p(k + 1:n), k) * W(p(k), k + 1:end));
	end
	
end

function [j] = piv_strat(W, k, p, d, strat)

	n = size(W, 1);
	switch strat
		case 'naive'
			j = find(W(p(k:n), k));
			j = j(1) + (k - 1);
		case 'partial'
			[~, j] = max(abs(W(p(k:n), k)));
			j = j + (k - 1);
		case 'scaled'
			j = max(abs(W(p(k:n), k)) ./ d(k:n));
			j = j + (k - 1);
		otherwise
			j = k;
	end
	
end
