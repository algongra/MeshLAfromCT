function mesh = Mesh_CheckNormalsOrientation(mesh)
% Checks if the normal orientation is inwards or outwards by comparing the total
% surface with vertices displaced towards the normal direction by the sqrt of
% the mean area
    v_ori = mesh.v;
    f = mesh.f;

    % Compute original areas and total area
    [areas,normals,vprods] = compute_vecprod(v_ori,f);
    displ = sqrt(mean(areas));
    TotArea_ori = sum(areas);

    % Initialize counter and accumulation variables
    counter = zeros([size(v_ori,1) 1]);
    v_displ = zeros(size(v_ori));

    % Compute displacement vector by accumulating normal contribution from each
    % face, then normalize using counter variable
    for i_f = 1:size(f,1)
        v_displ(f(i_f,1),:) = v_displ(f(i_f,1),:) + normals(i_f,:)*displ;
        counter(f(i_f,1)) = counter(f(i_f,1)) + 1;
        v_displ(f(i_f,2),:) = v_displ(f(i_f,2),:) + normals(i_f,:)*displ;
        counter(f(i_f,2)) = counter(f(i_f,2)) + 1;
        v_displ(f(i_f,3),:) = v_displ(f(i_f,3),:) + normals(i_f,:)*displ;
        counter(f(i_f,3)) = counter(f(i_f,3)) + 1;
    end
    %v_displ = v_displ./counter;
    v_displ = bsxfun(@rdivide,v_displ,counter);

    [areas,normals,vprods] = compute_vecprod(v_ori+v_displ,f);
    TotArea_dis = sum(areas);


    if TotArea_dis > TotArea_ori
        mesh.norm_out = true;
    else
        mesh.norm_out = false;
    end
return
