function result = compMultiFoxHIntegrand(y, params)
z = params{1};  mn = params{2};pq = params{3};
c = params{4};  d = params{5};  a = params{6};  b = params{7};
m = mn(1,:);    n = mn(2,:);    p = pq(1,:);    q = pq(2,:);

npoints =  size(y,1); dims = size(y,2); s = 1j * y;
lower = zeros(1,dims);  upper = zeros(1,dims);
for dim_l= 1:dims
    if ~isempty(b{dim_l})
        bj = b{dim_l}(1,:);        Bj = b{dim_l}(2,:);
        bj = bj(1:m(dim_l+1));        Bj = Bj(1:m(dim_l+1));
        lower(dim_l) = -min(bj / Bj);
    else
        lower(dim_l) = -100;
    end
     if ~isempty(a{dim_l})
        aj = a{dim_l}(1,:);        Aj = a{dim_l}(2,:);
        aj = aj(1:n(dim_l+1));        Aj = Aj(1:n(dim_l+1));
        upper(dim_l) = min((1-aj) / Aj);
    else
        upper(dim_l) = 1;
     end
end
mindist = norm(upper - lower);  sigs = 0.5 * (upper + lower);
for j = 1:n(1)
    num = 1 - c(j,1) - sum(c(j,2:end) .* lower);
    cnorm = norm(c(j,2:end),'fro');
    newdist = abs(num) / (cnorm + eps);
    if newdist < mindist
        mindist = newdist;
        sigs = lower + 0.5 * num .* c(j,2:end) ./ (cnorm * cnorm);
    end
end
s = s+sigs;
s1 = [ones(npoints, 1), s];
coder.varsize('prod_gam_num');  coder.varsize('prod_gam_denom');

prod_gam_num = (1 + 0i);    prod_gam_denom = (1 + 0i);
for j = 1:n(1)
    prod_gam_num = prod_gam_num .* ggamma(1-sum(bsxfun(@times, s1, c(j,:)), 2));
end
for j = 1:q(1)
    prod_gam_denom = prod_gam_denom .* ggamma(1-sum(bsxfun(@times, s1, d(j,:)), 2));
end
for j = 1+n(1):p(1)
    prod_gam_denom = prod_gam_denom .* ggamma(sum(bsxfun(@times, s1, c(j,:)), 2));
end
for dim_l = 1:dims
    for j = 1:n(dim_l + 1)
        prod_gam_num = prod_gam_num .* ggamma(1 - a{dim_l}(1,j) - a{dim_l}(2,j) * s(:, dim_l));
    end
    for j = 1:m(dim_l + 1)
        prod_gam_num = prod_gam_num .* ggamma(b{dim_l}(1,j) + b{dim_l}(2,j) * s(:, dim_l));
    end
    for j = 1+n(dim_l + 1) : p(dim_l + 1)
        prod_gam_denom = prod_gam_denom .* ggamma(a{dim_l}(1,j) + a{dim_l}(2,j) * s(:, dim_l));
    end
    for j = 1+m(dim_l + 1) : q(dim_l + 1)
        prod_gam_denom = prod_gam_denom .* ggamma(1 - b{dim_l}(1,j) - b{dim_l}(2,j) * s(:, dim_l));
    end
end 
zs = z.^(-s);
result = (prod_gam_num ./ prod_gam_denom) .* prod(zs, 2) ./ (2 * pi)^dims;
end