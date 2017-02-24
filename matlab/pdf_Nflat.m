function p = pdf_Nflat(mu, sig, m1, m2, data_N)
    p = bsxfun(@times, (sig').^(-data_N), exp( bsxfun(@rdivide, -(m2 -2*m1*mu +data_N*(mu.^2))/2, (sig').^2) ));
end
