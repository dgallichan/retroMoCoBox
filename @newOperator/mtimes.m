function res = mtimes(a,b)
% res = mtimes(FT, x)
%

if a.adjoint
    res = a.At(b);
else
    res = a.A(b);
end



    
