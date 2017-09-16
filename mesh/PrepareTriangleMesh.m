function mesh = PrepareTriangleMesh(points,elements)
% Prepare triangle mesh

mesh.nop = size(points,1);
mesh.noe = size(elements,1);
mesh.p = points;
mesh.e = elements;

%% For each node, get the triangles containing this node
maxt = 20; % max number of triangles per node
t = zeros(mesh.nop,maxt);
n = zeros(mesh.nop,1);
e1 = elements(:,1); e2 = elements(:,2); e3 = elements(:,3);
for i = 1:mesh.noe
    n(e1(i)) = n(e1(i))+1;
    t(e1(i),n(e1(i))) = i;
    
    n(e2(i)) = n(e2(i))+1;
    t(e2(i),n(e2(i))) = i;
    
    n(e3(i)) = n(e3(i))+1;
    t(e3(i),n(e3(i))) = i;
end
nmax = max(n);
mesh.ntri = t(:,1:nmax);
mesh.ntri_n = n;

end