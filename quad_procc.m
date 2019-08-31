
function [mean_nz, mean_var,rept ] = quad_procc(crp_fret, crp_procc)
% Mean of the cell part from the cropped region
    idx = find(crp_procc==1); 
    val = crp_fret(idx);
    val = val(:);
    mean_var = sum(val)/length(idx);
    mean_nz = nz_mean(val,1);
    rept = val;
end