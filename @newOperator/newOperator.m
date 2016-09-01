function  res = newOperator(A,At)

% attempt at making operator that uses given A and At as the mtimes and
% mtimes (adjoint) functions
%
% Daniel Gallichan, October 2012

res.adjoint = 0;
res.A = A;
res.At = At;
res = class(res,'newOperator');

