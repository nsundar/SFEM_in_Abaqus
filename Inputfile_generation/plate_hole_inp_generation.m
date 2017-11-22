clear all
close all
clc

global L D

domain = @PlatewHoleDomain ;
L=1;D=1;
a = 0.25 ; % radius of the hole

% load polygonal mesh....
load PlateHolemesh.mat

numelem=size(element,1);
numnode=size(node,1);

% flags to print node numbers and element numbers
pflag.pMesh = 'yes' ;

% Mesh to plot the polygonal mesh...
if strcmp(pflag.pMesh,'yes')
    clf; axis equal; axis on; hold on;
    Element = Vertices(1:numelem)';                 %Only plot the first block
    MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
    PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
    ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
    ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
    patch('Faces',ElemMat,'Vertices',node,'FaceColor','w'); pause(1e-6)
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname=strcat('PlateHole.inp')

% abaqus file input name
fileID = fopen(fname,'w') ;

% material property...
E = 1000; nu = 0.3;
stressState='planestress';
ndof = 2 ;

props = [1 E nu];

% get boundary nodes
botNodes = find(node(:,2) == min(node(:,2))) ;
rightNodes = find(node(:,1) == max(node(:,1))) ;
topNodes = find(node(:,2) == max(node(:,2))) ;
leftNodes = find(node(:,1) == min(node(:,1))) ;

% the nodes are arranged in increasing node numbers...hence sort them
% out...
[t1,t2] = sort(node(rightNodes,2)) ;
rightNodes = rightNodes(t2) ;

% sort the top nodes
[t1,t2] = sort(node(topNodes,1)) ;
topNodes = topNodes(t2);

% write the connectivity for the rightNodes
nume = ( length(rightNodes) - 1)/(nnode - 1) ;

cnt = 0 ;
for in = 1:nume
    
    % loop over the number of nodes
    for inode = 1:nnode
        cnt = cnt + 1 ;
        rightEdge(in,inode) = rightNodes(cnt) ;
    end
    cnt = cnt - 1 ;
end

nume = ( length(topNodes) - 1)/(nnode - 1) ;
cnt = 0 ;
for in = 1:nume
    
    % loop over the number of nodes
    for inode = 1:nnode
        cnt = cnt + 1 ;
        topEdge(in,inode) = topNodes(cnt) ;
    end
    cnt = cnt - 1 ;
end


%---------- force vector computation
QuadType = 'Legendre' ;

% get integration points..
if strcmp(QuadType, 'Legendre')
W=[0.347854845137454 0.652145154862546 0.652145154862546 0.347854845137454]';
Q=[-0.861136311594053 -0.339981043584856 0.339981043584856 0.861136311594053]';
end

% initialization of force vector 
fvec = zeros(numnode*2,1);
for iel = 1:size(rightEdge,1)
    
    % get current element connectivity...
    econ = rightEdge(iel,:) ;
    nn = length(econ) ;
    
    gindx = econ.*2-1 ;
    gindy = econ.*2 ;
    
    % the corresponding nodal coordinates
    xynodes = node(econ,:) ;
    
    % loop over integration points...
    for igp = 1:size(W,1)
        pt = Q(igp,:) ;
        
        % One dimensional shape function
        shp(1,1)=0.5*(1-pt);
        shp(1,2)=0.5*(1+pt);
        shp(2,1)=-0.5;
        shp(2,2)=0.5;
        
        N = shp(1,:); dNdxi = shp(2,:);
        
        Ja = dNdxi*xynodes(:,2) ;
        gpt = N*xynodes ;
        
        [exact_stress,displace] = exact_plate_hole(gpt,a,stressState,E,nu) ;
        tx = exact_stress(1) ;
        ty = exact_stress(3) ;
        
        fvec(gindx,1) = fvec(gindx,1) + N'*tx*det(Ja)*W(igp) ;
        fvec(gindy,1) = fvec(gindy,1) + N'*ty*det(Ja)*W(igp) ;
    end
end


for iel = 1:size(topEdge,1)
    
    % get current element connectivity...
    econ = topEdge(iel,:); 
    nn = length(econ) ;
    
    gindx = econ.*2-1 ;
    gindy = econ.*2 ;
    
    % the corresponding nodal coordinates
    xynodes = node(econ,:) ;
    
    % loop over integration points...
    for igp = 1:size(W,1)
        pt = Q(igp,:) ;
        
      % One dimensional shape function
        shp(1,1)=0.5*(1-pt);
        shp(1,2)=0.5*(1+pt);
        shp(2,1)=-0.5;
        shp(2,2)=0.5;
        
        N = shp(1,:); dNdxi = shp(2,:);
        
        Ja = dNdxi*xynodes(:,1) ;
        gpt = N*xynodes ;
        
        if Ja < 0
            error('Jacobian less than zero')
        end
        
        [exact_stress,displace] = exact_plate_hole(gpt,a,stressState,E,nu) ;
        tx = exact_stress(3) ;
        ty = exact_stress(2) ;
        
        fvec(gindx,1) = fvec(gindx,1) + N'*tx*det(Ja)*W(igp) ;
        fvec(gindy,1) = fvec(gindy,1) + N'*ty*det(Ja)*W(igp) ;
    end
end

%==========================================================================
% find the number of sides of the element...and group them based on the
% number of sides...
m3 = 1; m4 = 1; m5 = 1; m6 = 1; m7 = 1; m8 = 1;
cnt = 1;
for iel = 1:numelem
    econ = element{iel};
    nn = length(econ);
    
    if nn == 3
        % three sided element
        P3(m3,:) = econ ;
        m3 = m3 + 1;
    elseif nn == 4
        % four sided element
        P4(m4,:) = econ ;
        m4 = m4 + 1 ;
    elseif nn == 5
        % five sided element
        P5(m5,:) = econ ;
        m5 = m5 + 1 ;
    elseif nn == 6
        % six sided element
        P6(m6,:) = econ;
        m6 = m6 + 1 ;
    elseif nn == 7
        % seven sided element
        P7(m7,:) = econ;
        m7 = m7 + 1 ;
    elseif nn == 8
        % eight sided element
        P8(m8,:) = econ;
        m8 = m8 + 1;
    end

    numSide(cnt,1) = nn;
    cnt = cnt + 1 ;
end
%    save (matname2, 'P3','P4', 'P5','P6','P7','P8');

maxEleSide = max(numSide);
minEleSide = min(numSide);

uniNumSide = unique(numSide);

%--------------- write to the file
fprintf(fileID,'*Heading\n');
fprintf(fileID,'4Nodes\n');
fprintf(fileID,'** Job name: Job-1 Model name: Plate\n');
fprintf(fileID,'** Generated by: Abaqus/CAE 6.12-3\n');
fprintf(fileID,'*Preprint, echo=NO, model=NO, history=NO, contact=NO\n');
fprintf(fileID,'**\n');
fprintf(fileID,'** PARTS\n');
fprintf(fileID,'**\n');
fprintf(fileID,'**Part, name=Part-1\n');
fprintf(fileID,'*Node,NSET=NALL\n');

% write the nodal coordinates...
for in = 1:numnode
    fprintf(fileID,'%d,  %16.14f,   %16.14f\n',in,node(in,1),node(in,2));
end

% element information
for iel = 1:length(uniNumSide)
    nel = uniNumSide(iel);
    fprintf(fileID,'*User element, nodes=%d, type=U%d, properties=3,coordinates=2, variables=7\n',nel,nel);
    fprintf(fileID,'1,2\n');
end

numElCnt = 1;
for iel = 1:length(uniNumSide)
    
    % get the number of nodes
    nel = uniNumSide(iel);
    
    if nel == 3
        Elmset = 'three';
        % write down for 3 sided element
        fprintf(fileID,'*Element, type=U%d,ELSET=%s\n',nel,Elmset);
        
        % loop over the elements in this category
        for iel = 1:size(P3,1)
            econ = P3(iel,:);
            fprintf(fileID,'%d,%d,%d,%d\n',numElCnt,econ);
            numElCnt= numElCnt + 1;
        end
        
        % material property for this element type
        fprintf(fileID,'*UEL Property, ELSET=%s\n',Elmset);
        fprintf(fileID,'%d, %6.4e,  %3.2f\n',props);
        
    elseif nel == 4
        Elmset = 'four';
        % write down for 3 sided element
        fprintf(fileID,'*Element, type=U%d,ELSET=%s\n',nel,Elmset);
        
        % loop over the elements in this category
        for iel = 1:size(P4,1)
            econ = P4(iel,:) ;
            fprintf(fileID,'%d,%d,%d,%d,%d\n',numElCnt,econ);
            numElCnt= numElCnt + 1;
        end
        
        % material property for this element type
        fprintf(fileID,'*UEL Property, ELSET=%s\n',Elmset);
        fprintf(fileID,'%d, %6.4e,  %3.2f\n',props);
        
    elseif nel == 5
        Elmset = 'five';
        % write down for 3 sided element
        fprintf(fileID,'*Element, type=U%d,ELSET=%s\n',nel,Elmset);
        
        % loop over the elements in this category
        for iel = 1:size(P5,1)
            econ = P5(iel,:) ;
            fprintf(fileID,'%d,%d,%d,%d,%d,%d\n',numElCnt,econ);
            numElCnt= numElCnt + 1;
        end
        
        % material property for this element type
        fprintf(fileID,'*UEL Property, ELSET=%s\n',Elmset);
        fprintf(fileID,'%d, %6.4e,  %3.2f\n',props);
        
    elseif nel == 6
        Elmset = 'six';
        % write down for 3 sided element
        fprintf(fileID,'*Element, type=U%d,ELSET=%s\n',nel,Elmset);
        
        % loop over the elements in this category
        for iel = 1:size(P6,1)
            econ = P6(iel,:) ;
            fprintf(fileID,'%d,%d,%d,%d,%d,%d,%d\n',numElCnt,econ);
            numElCnt= numElCnt + 1;
        end
        
        % material property for this element type
        fprintf(fileID,'*UEL Property, ELSET=%s\n',Elmset);
        fprintf(fileID,'%d, %6.4e,  %3.2f\n',props);
        
    elseif nel == 7
        Elmset = 'seven';
        % write down for 3 sided element
        fprintf(fileID,'*Element, type=U%d,ELSET=%s\n',nel,Elmset);
        
        % loop over the elements in this category
        for iel = 1:size(P7,1)
            econ = P7(iel,:) ;
            fprintf(fileID,'%d,%d,%d,%d,%d,%d,%d,%d\n',numElCnt,econ);
            numElCnt= numElCnt + 1;
        end
        
        % material property for this element type
        fprintf(fileID,'*UEL Property, ELSET=%s\n',Elmset);
        fprintf(fileID,'%d, %6.4e,  %3.2f\n',props);
        
    elseif nel == 8
        Elmset = 'eight';
        % write down for 3 sided element
        fprintf(fileID,'*Element, type=U%d,ELSET=%s\n',nel,Elmset);
        
        % loop over the elements in this category
        for iel = 1:size(P8,1)
            econ = P8(iel,:) ;
            fprintf(fileID,'%d,%d,%d,%d,%d,%d,%d,%d,%d,\n',numElCnt,econ);
            numElCnt= numElCnt + 1;
        end 
        
        % material property for this element type
        fprintf(fileID,'*UEL Property, ELSET=%s\n',Elmset);
        fprintf(fileID,'%d, %6.4e,  %3.2f\n',props);
    end
end


% create node sets...
fprintf(fileID,'*Nset, nset=Set-1\n');
fprintf(fileID,'%d,  %d,%d\n',1,  numnode, 1);
fprintf(fileID,'*Elset, elset=Set-1, generate\n');
fprintf(fileID,'%d,%d,%d\n',1,numelem,1);

fprintf(fileID,'*Step, name=Step-1,NLGEOM=YES\n');
fprintf(fileID,'*Static\n');
fprintf(fileID,'1., 1., 1., 1.\n');
fprintf(fileID,'** \n');
fprintf(fileID,'** BOUNDARY CONDITIONS\n');
fprintf(fileID,'** \n');
fprintf(fileID,'** Name: BC-1 Type: Displacement/Rotation\n');
fprintf(fileID,'*Boundary\n');

for i=1:length(botNodes)
    
    cn = botNodes(i);
    
    x = node(botNodes(i),1); y = node(botNodes(i),2);
    gpt = [x y];
    [exact_stress,displace] = exact_plate_hole(gpt,a,stressState,E,nu) ;
    
    fixed_u=displace(1) ;
    fixed_v=displace(2) ;

    fprintf(fileID,'%d,  %d,, %16.16f\n',cn,1,fixed_u);
    fprintf(fileID,'%d,  %d,, %16.16f\n',cn,2,fixed_v);

end

for i=1:length(leftNodes)
        
    cn = leftNodes(i);
    
    x = node(leftNodes(i),1); y = node(leftNodes(i),2);
    gpt = [x y];
    [exact_stress,displace] = exact_plate_hole(gpt,a,stressState,E,nu) ;
    
    fixed_u=displace(1) ;
    fixed_v=displace(2) ;

    fprintf(fileID,'%d,  %d,, %16.16f\n',cn,1,fixed_u);
    fprintf(fileID,'%d,  %d,, %16.16f\n',cn,2,fixed_v);

end

% force boundary condition
fprintf(fileID,'*CLOAD\n');
for in = 1:length(rightNodes)
    
    tin = rightNodes(in);
    
    tfx = fvec(2*tin-1,1);
    tfy = fvec(2*tin,1);
    
    fprintf(fileID,'%5d, %1d, %16.16f\n',tin,1,tfx) ;
    fprintf(fileID,'%5d, %1d, %16.16f\n',tin,2,tfy) ;
end

for in = 1:length(topNodes)
    
    tin = topNodes(in);
    
    tfx = fvec(2*tin-1,1);
    tfy = fvec(2*tin,1);
    
    fprintf(fileID,'%5d, %1d, %16.16f\n',tin,1,tfx) ;
    fprintf(fileID,'%5d, %1d, %16.16f\n',tin,2,tfy) ;
end

%--- output request
fprintf(fileID,'** OUTPUT REQUESTS\n');
fprintf(fileID,'** \n');
fprintf(fileID,'*Restart, write, frequency=0\n');
fprintf(fileID,'** \n');
fprintf(fileID,'** FIELD OUTPUT: F-Output-1\n');
fprintf(fileID,'** \n');
fprintf(fileID,'*Output, field, variable=PRESELECT\n');
fprintf(fileID,'** \n');
fprintf(fileID,'** HISTORY OUTPUT: H-Output-1\n');
fprintf(fileID,'** \n');
fprintf(fileID,'*Output, history, variable=PRESELECT\n');
fprintf(fileID,'*NODE PRINT, NSET=NALL\n');
fprintf(fileID,'U1,U2,CF1\n');
fprintf(fileID,'*End Step\n');

% close the file...all is done :)
fclose(fileID);
