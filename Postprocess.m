
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%show total displacement
figure('Renderer', 'opengl');%hardware acceleration
for tg = 0:1:0
    node = load(['resuNode_',num2str(tg),'.txt']);
    elem = load(['resuElem_',num2str(tg),'.txt']) + 1;
    disp = load(['resuDisp_',num2str(tg),'.txt']);
    face = zeros(6*size(elem,1),4);
    for ti=1:1:size(elem,1)
        face(6*(ti-1)+1,1:4) = [elem(ti,1),elem(ti,2),elem(ti,3),elem(ti,4)];
        face(6*(ti-1)+2,1:4) = [elem(ti,5),elem(ti,6),elem(ti,7),elem(ti,8)];
        face(6*(ti-1)+3,1:4) = [elem(ti,1),elem(ti,4),elem(ti,8),elem(ti,5)];
        face(6*(ti-1)+4,1:4) = [elem(ti,2),elem(ti,3),elem(ti,7),elem(ti,6)];
        face(6*(ti-1)+5,1:4) = [elem(ti,1),elem(ti,5),elem(ti,6),elem(ti,2)];
        face(6*(ti-1)+6,1:4) = [elem(ti,4),elem(ti,8),elem(ti,7),elem(ti,3)];
    end
    patch('Vertices', node(:,:), 'Faces', face(:,:), ...
        'FaceVertexCData', sqrt(disp(:,1).^2+disp(:,2).^2+disp(:,3).^2), 'FaceColor', 'interp');
    hold on;
end
colorbar;
colormap(jet);
axis equal;
view(30,40);
xlabel('x');
ylabel('y');
zlabel('z');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%show contact pressure
figure('Renderer', 'opengl');%hardware acceleration
for tg = 0:1:0
    intePoin=load(['resuIpga_',num2str(tg),'.txt']);
    scatter3(intePoin(:,1), intePoin(:,2), intePoin(:,3),25,intePoin(:,4),'filled');
    hold on;
end
axis equal;
colorbar;
colormap(jet);
