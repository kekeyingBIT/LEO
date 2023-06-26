function [ob_nIBI] = ParsePilot(rec_sig,N_block,l_nIBI,leo,fr)

ob_nIBI = zeros(l_nIBI,prod(leo.N_r),N_block,leo.N_p);
for np = 1:leo.N_p
    for blk = 1:N_block
        % parse nIBI part for JADCE
        ts_idx = (blk-1)*fr.n_symbol_pblock +...
            ( (blk-1)*fr.n_ts_pblock+1 : blk*fr.n_ts_pblock );
        ob_TS = rec_sig(ts_idx,:,np);
        ob_nIBI(:,:,blk,np) = ob_TS(end-l_nIBI+1:end,:);
    end
end

end
