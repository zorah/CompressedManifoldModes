function [ A ] = vertex_mass( M )
%VERTEX_MASS Calculates the vertex mass for every vertex of M in M.VERT.

num_vertices = size(M.VERT, 1);
faces = M.TRIV;

tri_areas = calc_tri_areas(M);

% compute area of Voronoi cells
A = zeros(num_vertices,1);
for j = 1:num_vertices
    for i = 1:3
        ind_j = find(faces(:,i) == j);
        for l = 1:size(ind_j,1)
            face_index = ind_j(l);
            A(j) = A(j) + tri_areas(face_index)/3;
        end        
    end
end

end

