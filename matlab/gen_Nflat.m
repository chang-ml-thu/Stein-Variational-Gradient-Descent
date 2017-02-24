function gen_Nflat( N, mu, sig )
    fid = fopen(sprintf('data_Nflat_%d_%.0f_%.0f.txt', N, mu, sig), 'w');
    fprintf(fid, '%d %.6f %.6f\n', N, mu, sig);
    fprintf(fid, ' %.6f', mu + sig * randn(1, N));
    fprintf(fid, '\n');
    fclose(fid);
end