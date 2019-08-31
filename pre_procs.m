function [ op_im ] = pre_procs( im,op_var,thresh )
%PRE_PROCS Pre Processing the Cell Masking

if size(im,3)>1
    im=rgb2gray(im);
else
    im=im;
end
 switch nargin
        case 3
           im=imbinarize(im,thresh);
        case 2
            im=imbinarize(im);
           
    end




im=imfill(im,'holes');
im=bwareaopen(im,op_var);
op_im=im;

end

