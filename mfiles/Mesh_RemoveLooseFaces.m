function keep_face = Mesh_RemoveLooseFaces(f,keep_face)

    temp = keep_face;

    disp('Checking for loose faces...');
    %{
    % AGG_BEG: Saving this lines because they have useful comments according to
    %          Lorenzo.
    %          The new vectorized version is doing the same
    % List of indices that should be on the wall
    discard_face_ind = find(~keep_face);
    
    discard_face = reshape(f(repmat(~keep_face,[1 3])),...
        [length(discard_face_ind) 3]);

    counter = 0;
    tic
    for i_f = 1:length(discard_face_ind)              
        % Compare each f in the inlet (discard_face) with all the faces in
        % the inlet itself
        matching_faces1 = bsxfun(@eq,f(discard_face_ind(i_f),1),discard_face);
        matching_faces2 = bsxfun(@eq,f(discard_face_ind(i_f),2),discard_face);
        matching_faces3 = bsxfun(@eq,f(discard_face_ind(i_f),3),discard_face);
        
        % Count how many vertices for this face i_f are only in that one
        % face
        lone_vertices = double(sum(matching_faces1(:))<2) +...
            double(sum(matching_faces2(:))<2) +...
            double(sum(matching_faces3(:))<2);
        
        % Discard from wall if they have 2 lone_vertices (i.e. only 1
        % vertex connects that face to the rest of the faces)
        if lone_vertices >= 2
            keep_face(discard_face_ind(i_f)) = ~keep_face(discard_face_ind(i_f));
            counter = counter+1;
            disp(['Face ' num2str(discard_face_ind(i_f)) ' swapped']);
        end
    end
    if counter
        disp(['Removed ' num2str(counter) ' loosely connected faces']);
    end
    % AGG_END: Saving this lines because they have useful comments according to
    %          Lorenzo.
    %          The new vectorized version is doing the same
    %}
    

    keep_face = temp;
    %%%%%%%%%%%%%%%%%%%%%%% NEWER VECTORIZED VERSION %%%%%%%%%%%%%%%%%%%%%%    
    % List of indices that should be on th wall
    discard_face_ind = find(~keep_face);
    
    discard_face = reshape(f(repmat(~keep_face,[1 3])),...
        [length(discard_face_ind) 3]);
    
    matching_faces1 = bsxfun(@eq,f(discard_face_ind,1),reshape(discard_face,[1 size(discard_face)]));
    matching_faces2 = bsxfun(@eq,f(discard_face_ind,2),reshape(discard_face,[1 size(discard_face)]));
    matching_faces3 = bsxfun(@eq,f(discard_face_ind,3),reshape(discard_face,[1 size(discard_face)]));
       
    lone_vertices = double(sum(sum(matching_faces1,3),2)<2) +...
            double(sum(sum(matching_faces2,3),2)<2) +...
            double(sum(sum(matching_faces3,3),2)<2);
        
%     lone_vertices = double(sum(sum(matching_faces1,3),1)<2) +...
%             double(sum(sum(matching_faces2,3),1)<2) +...
%             double(sum(sum(matching_faces3,3),1)<2);
        
    keep_face(discard_face_ind( lone_vertices>=2 )) = ...
        ~ keep_face(discard_face_ind( lone_vertices>=2 ));
          
    if sum(lone_vertices>=2)
        disp(['Removed ' num2str(sum(lone_vertices>=2)) ' loosely connected faces']);
        discard_face_ind( lone_vertices>=2 )
    end
end
