function [X,Y,Z] = imgrid(fov,dim)
% fov: 3-element vector containing x,y,z field of view
% dim: 3-element vector containing x,y,z dimensions
x = linspace(-fov(1)/2,fov(1)/2,dim(1));
y = linspace(-fov(2)/2,fov(2)/2,dim(2));
z = linspace(-fov(3)/2,fov(3)/2,dim(3));
[X,Y,Z] = ndgrid(x,y,z);
end

