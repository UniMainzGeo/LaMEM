function VisualizeInitialMesh(X,Y,Z)
% Visualize the initial Mesh
%




plot3(squeeze(X(:,:,1  )),squeeze(Y(:,:,1  )),squeeze(Z(:,:,1  )),'k',squeeze(X(:,:,1  )).',squeeze(Y(:,:,1  )).',squeeze(Z(:,:,1  )).','k')

hold on
plot3(squeeze(X(:,:,end)),squeeze(Y(:,:,end)),squeeze(Z(:,:,end)),'k',squeeze(X(:,:,end)).',squeeze(Y(:,:,end)).',squeeze(Z(:,:,end)).','k')

hold on
plot3(squeeze(X(:,1,:  )),squeeze(Y(:,1,:  )),squeeze(Z(:,1,:)),'k',squeeze(X(:,1,:)).',squeeze(Y(:,1,:)).',squeeze(Z(:,1,:)).','k')

hold on
plot3(squeeze(X(:,end,:  )),squeeze(Y(:,end,:  )),squeeze(Z(:,end,:)),'k',squeeze(X(:,end,:)).',squeeze(Y(:,end,:)).',squeeze(Z(:,end,:)).','k')

hold on
plot3(squeeze(X(1,:,:  )),squeeze(Y(1,:,:  )),squeeze(Z(1,:,:)),'k',squeeze(X(1,:,:)).',squeeze(Y(1,:,:)).',squeeze(Z(1,:,:)).','k')

hold on
plot3(squeeze(X(end,:,:  )),squeeze(Y(end,:,:  )),squeeze(Z(end,:,:)),'k',squeeze(X(end,:,:)).',squeeze(Y(end,:,:)).',squeeze(Z(end,:,:)).','k')

axis equal
