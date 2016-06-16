function [ nplanes ] = GetNPlanes(min_z, max_z)

plane_info = load('Minerva_Planes.txt');

fid_modules = plane_info(((plane_info(:,6) >= min_z) & (plane_info(:,6) <= max_z)),4);

nplanes = length(fid_modules);

end

