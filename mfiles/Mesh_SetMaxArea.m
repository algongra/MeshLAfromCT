function mesh = Mesh_SetMaxArea(mesh,MaxArea)
% Split in 3 smaller triangles, faces bigger than a MaxArea threshold

v = mesh.v;
f = mesh.f;

if isfield(mesh,'v_edge_ordered')
    v_edge_ordered = mesh.v_edge_ordered;
    ShiftVEdge = true;
else
    ShiftVEdge = false;
end

if isfield(mesh,'face_type')
    face_type = mesh.face_type;
    UseFaceType = true;
else
    UseFaceType = false;
end

InitialTrianglesNumber = size(f,1);
disp(['Initial n. of vertices: ' num2str(size(v,1))]);
disp(['Initial n. of triangles: ' num2str(InitialTrianglesNumber)]);

% Compute all areas of the triangles
areas = elemvolume(v,f);
i_f = 0;
split = 0;
while i_f<length(f(:,1))
    i_f = i_f + 1;
    if areas(i_f) >= MaxArea
        split = split + 1;
        % Create new vertex
        new_v = mean(v(f(i_f,:),:),1);
        %%%% POSSIBLE BUG ???? %%%%
        last_old_v = min(f(i_f));

        % last_old_v = min(f(i_f,:));
        
        % Insert it right after the smallest vertex of the triangle
        v = [v(1:last_old_v,:);...
            new_v;...
            v(last_old_v+1:end,:)];
        
        %             % Remove it from the list (needed?)
        %             f(i_f,:) = [0 0 0];
        
        % Update all the triangle vertices since we added one vertex
        f(f>last_old_v) = f(f>last_old_v)+1;        

        if ShiftVEdge
            %%%%%%%% TESTING: Trying to shift the v_edge_ordered
            for i_edge = 1:length(v_edge_ordered)
                v_edge_ordered{i_edge}(v_edge_ordered{i_edge}>last_old_v) =... 
                    v_edge_ordered{i_edge}(v_edge_ordered{i_edge}>last_old_v)+1;
            end
        end

        
        % Set aside the triangle that has to be split
        new_f = f(i_f,:);
        if UseFaceType
            new_face_type = face_type(i_f);
        end
        
        % Add 3 triangle in place of the original one (decrease a lot
        % mesh quality...)
        f = [f(1:i_f-1,:);...
            [new_f' circshift(new_f,[0,-1])' repmat(last_old_v+1,[3 1])];
            f(i_f+1:end,:)];
        
        if UseFaceType
            face_type = [face_type(1:i_f-1);...
                repmat(new_face_type,[3 1]);...
                face_type(i_f+1:end)];
        end
        
        % f(i_f-2:i_f+4,:)
        % Add fake zero areas too
        areas = [areas(1:i_f-1); zeros(3,1); areas(i_f+1:end)];
    end
end
disp(['Triangles split in three: ' num2str(split)]);
if split > InitialTrianglesNumber/10
    warning('WARNING: I''m splitting more than 10%% of the initial triangles! Mesh quality will be poor. Increase MaxArea or keep more triangles!');
end

mesh.v = v;
mesh.f = f;
if ShiftVEdge
    mesh.v_edge_ordered = v_edge_ordered;
end
if UseFaceType
    mesh.face_type = face_type;
end
return
