function out = gs(data,smthfact, interpfact)
%function out = gs(data,smthfact, interpfact)
% smooth 2D data with a gaussian
if nargin < 3
    interpfact = 2;
end

if nargin < 2
    smthfact = 0.5;
end

dimX = size(data,1);
dimY = size(data,2);

if smthfact
%     data = filter2(fspecial('gaussian',9,smthfact),data);
    data = filter2(gauss(smthfact,max(3,round(smthfact*3))),data);
end

if interpfact
    data = imresize(data,round([dimX*interpfact dimY*interpfact]),'bicubic');
end

out = data;