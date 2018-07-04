function [ A, L ] = calc_cotan_mass( N )
%CALC_COTAN_MASS Calculates mass (A) and stiffness (L) matrix of N.
% based on code from Emanuele Rodolà

vertices = N.VERT;
faces = N.TRIV;

num_vertices = size(vertices,1);
num_faces = size(faces,1);

% check consistency
adjacency_matrix = sparse([faces(:,1); faces(:,2); faces(:,3)], ...
    [faces(:,2); faces(:,3); faces(:,1)], ...
    ones(3 * num_faces, 1), ...
    num_vertices, num_vertices, 3 * num_faces);
if any(any(adjacency_matrix > 1))
    error('Triangles must be oriented consistently.')
end
clear adjacency_matrix;

% first compute inner face angles and squared edge lengthes
pp = zeros(num_faces,3);
qq = zeros(num_faces,3);
angles = 0*faces;
squared_edge_length = 0*faces;

for i=1:3
    i1 = mod(i-1,3)+1;
    i2 = mod(i  ,3)+1;
    i3 = mod(i+1,3)+1;
    pp = vertices(faces(:,i2),:) - vertices(faces(:,i1),:);
    qq = vertices(faces(:,i3),:) - vertices(faces(:,i1),:);
    % normalize the vectors
    pp = pp ./ repmat( max(sqrt(sum(pp.^2,2)),eps), [1 3] );
    qq = qq ./ repmat( max(sqrt(sum(qq.^2,2)),eps), [1 3] );
    % compute angles
    angles(:,i1) = acos(sum(pp.*qq,2));
    squared_edge_length(:,i1) = sum((vertices(faces(:,i2),:) - vertices(faces(:,i3),:)).^2,2);
end
clear pp qq;

%then compute L
L = sparse(num_vertices,num_vertices);
for i=1:3
    i1 = mod(i-1,3)+1;
    i2 = mod(i  ,3)+1;
    i3 = mod(i+1,3)+1;
    L = L + sparse(faces(:,i1),faces(:,i2),-cot(angles(:,i3)),...
        num_vertices,num_vertices,num_faces);
end

%
L = 1/2 * (L + L');
L = sparse(1:num_vertices,1:num_vertices,-sum(L,2),num_vertices,num_vertices,...
    num_vertices) + L;

A = vertex_mass(N);
A = sparse(1:num_vertices, 1:num_vertices,A);

end

