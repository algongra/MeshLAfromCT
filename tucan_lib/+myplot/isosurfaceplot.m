function h1 = isosurfaceplot(X,Y,Z,var,C,threshold)

[f v] = isosurface(Y,X,Z,var,threshold);

h1 = patch('Faces',f,'Vertices',v);
isonormals(Y,X,Z,var,h1);
set(h1,'facecolor','interp','edgecolor','none');
set(gca, 'XDir', 'reverse');

isocolors(Y,X,Z,C,h1);

return
end
