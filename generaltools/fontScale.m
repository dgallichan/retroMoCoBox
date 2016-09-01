function fontScale(scale)
% function fontScale(scale)

H = gcf;
allText   = findall(H, 'type', 'text');
allAxes   = findall(H, 'type', 'axes');
allFont   = [allText; allAxes];
fontSize = get(allFont,'FontSize');

if ~iscell(fontSize)
    fontSize = num2cell(fontSize);
end
   
newFontSize = LocalScale(fontSize, scale, 2);
set(allFont,{'FontSize'},newFontSize);
set(allFont,'FontName','Open Sans')

function newArray = LocalScale(inArray, scale, minValue)
n = length(inArray);
newArray = cell(n,1);
for k=1:n
  newArray{k} = max(minValue,scale*inArray{k}(1));
end