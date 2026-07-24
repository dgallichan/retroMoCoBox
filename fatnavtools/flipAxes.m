function out = flipAxes(dataToFlip,flip_xyz) 

if isempty(flip_xyz)
    flip_xyz = [0 0 1]; % current best guess :)
end

out = dataToFlip;
if flip_xyz(1)
    out = flip(out,1);
end
if flip_xyz(2)
    out = flip(out,2);
end
if flip_xyz(3)
    out = flip(out,3);
end

