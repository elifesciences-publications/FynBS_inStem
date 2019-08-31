function [ op ] = nz_mean( ip_arr, dim )
%NZ_MEAN Summary of this function goes here
%   This function gives the value of nz mean of a ip array
ip_arr = double(ip_arr);
ip_arr(find(ip_arr==0))=NaN;
op=nanmean(ip_arr,dim);

end

