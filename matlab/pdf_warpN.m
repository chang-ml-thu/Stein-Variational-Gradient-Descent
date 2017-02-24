function p = pdf_warpN(w1, w2, m1, m2, data_N, sig0, sigx)
    w2sq = w2.^2;
    w1sq = w1.^2;
    mu = bsxfun(@plus, w1, w2sq');
    p = exp( - bsxfun(@plus, w1sq, w2sq') / (2 * sig0^2) - ( m2 - 2*m1*mu + data_N*mu.^2 ) / (2 * sigx^2) );
end
