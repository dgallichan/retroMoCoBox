function out = perlinNoise1D(Npts,weights)
% function out = perlinNoise1D(Npts,weights)
%
% I first heard about 'Perlin noise' from the MRF Nature paper
% (www.nature.com/articles/nature11971). I then made this Matlab code by
% translating the pseudo-code from http://freespace.virgin.net/hugo.elias/models/m_perlin.htm
% which used to be an article about Perlin noise - but now seems to be inaccessible. 
% According to the Wikipedia discussion around the Wiki article on Perlin noise, this isn't
% strictly-speaking Perlin noise because it is 'value' noise and not
% 'gradient' noise - but for most MRI purposes this probably doesn't
% matter. 

n = length(weights);
if n==1 
    n = weights;
    weights = 2.^[0:n-1];
end

xvals = linspace(0,1,Npts);

total = zeros(Npts,1);
for i = 1:n
    frequency = 2^(i-1);    
    this_Npts = round(Npts/frequency);
    if this_Npts > 1
        total = total + weights(i) * interp1(linspace(0,1,this_Npts).',rand(this_Npts,1),xvals.','pchip');
    else
        disp(['Maxed out at octave ' num2str(i)]);
    end
end
total = total-min(total);
total = total/max(total);
out = total;

end



