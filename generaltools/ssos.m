function out = ssos(in,dim)

if nargin > 1
    out = sqrt(sum(in.*conj(in),dim));
    
else
    
    if size(in,4)==1
        if size(in,3)==1
            out = sqrt(sum(in.*conj(in),2));
        else
            out = sqrt(sum(in.*conj(in),3));
        end
    else
        out = sqrt(sum(in.*conj(in),4));
    end
    
end