function [zr] = rectify(z)
    zr = z .* (z>0);
end