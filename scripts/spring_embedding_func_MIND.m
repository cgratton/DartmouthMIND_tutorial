function [h nodepos] = spring_embedding_func_COW(graph,colors,Kc,L0)
% Spring-embedding algorithm
%[h nodepos] = spring_embedding_func_new(graph,colors,Kc,L0)
%Kc: try 1
%L0:try 25
% Kamada-Kawai algorithm, see e.g., https://arxiv.org/pdf/1201.3011.pdf


nodes = size(graph,1);

for n = 1:nodes
    row = zeros(1,nodes);
    row(n) = 1;
    degree(n) = sum(graph(n,~row));
end

outnodeind = find(degree==0);

%Define circle
r = L0;
theta = 0:2*pi/nodes:2*pi-2*pi/nodes;
nodepos = [r*sin(theta); r*cos(theta)]';

% Find shortest path using dijkstra algorithm
D = dijkstra(logical(graph),((1-graph) .* logical(graph)));
D(isinf(D))=0;

% Find desirable length between nodes
maxd = max(D(:));
L = L0/maxd;
I = L.*D;

% Strength of springs between nodes
K = Kc./(D.^2);
K(isinf(K))=0;

xm = repmat(nodepos(:,1),1,nodes);
xm(logical(eye(nodes))) = 0;
ym = repmat(nodepos(:,2),1,nodes);
ym(logical(eye(nodes))) = 0;
xother = repmat(nodepos(:,1),1,nodes)';
xother(logical(eye(nodes))) = 0;
yother = repmat(nodepos(:,2),1,nodes)';
yother(logical(eye(nodes))) = 0;

Im = I;
Km = K;

dEdxmstep = Km.*((xm-xother)-((Im.*(xm-xother))./((xm-xother).^2 + (ym - yother).^2).^.5));
dEdxmstep(isnan(dEdxmstep))=0;
dEdxm = sum(dEdxmstep,2);
dEdymstep = Km.*((ym-yother)-((Im.*(ym-yother))./((xm-xother).^2 + (ym - yother).^2).^.5));
dEdymstep(isnan(dEdymstep))=0;
dEdym = sum(dEdymstep,2);
deltam = sqrt((dEdxm).^2+(dEdym).^2)';

count = 1;
while max(deltam)>.05
    string{count} = ['Iteration ' num2str(count)];
    if count==1; fprintf('%s',string{count}); else fprintf([repmat('\b',1,length(string{count-1})) '%s'],string{count}); end
    if count>10000
        break
    end
    [maxval maxi] = max(deltam);
    xm = repmat(nodepos(maxi,1),nodes-1,1);
    ym = repmat(nodepos(maxi,2),nodes-1,1);
    xother = nodepos(:,1);
    xother(maxi) = [];
    yother = nodepos(:,2);
    yother(maxi) = [];
    Im = I(maxi,:)';
    Im(maxi) = [];
    Km = K(maxi,:)';
    Km(maxi) = [];
    dEdxm = sum(Km.*((xm-xother)-((Im.*(xm-xother))./((xm-xother).^2 + (ym - yother).^2).^.5)));
    dEdym = sum(Km.*((ym-yother)-((Im.*(ym-yother))./((xm-xother).^2 + (ym - yother).^2).^.5)));
    
    d2Edxm2 = sum(Km.*(1-((Im.*(ym-yother).^2)./((xm-xother).^2+(ym-yother).^2).^(3/2))));
    d2Edym2 = sum(Km.*(1-((Im.*(xm-xother).^2)./((xm-xother).^2+(ym-yother).^2).^(3/2))));
    
    d2Edxmdym = sum(Km.*((Im.*(xm-xother).*(ym-yother))./((xm-xother).^2+(ym-yother).^2).^(3/2)));
    d2Edymdxm = sum(Km.*((Im.*(xm-xother).*(ym-yother))./((xm-xother).^2+(ym-yother).^2).^(3/2)));
    
    A = [d2Edxm2 d2Edxmdym; d2Edymdxm d2Edym2];
    B = [-dEdxm;-dEdym];
    
    X = linsolve(A,B);
    deltax = X(1);
    deltay = X(2);
    nodepos(maxi,1) = nodepos(maxi,1)+deltax;
    nodepos(maxi,2) = nodepos(maxi,2)+deltay;
    
    xm = repmat(nodepos(:,1),1,nodes);
    xm(logical(eye(nodes))) = 0;
    ym = repmat(nodepos(:,2),1,nodes);
    ym(logical(eye(nodes))) = 0;
    xother = repmat(nodepos(:,1),1,nodes)';
    xother(logical(eye(nodes))) = 0;
    yother = repmat(nodepos(:,2),1,nodes)';
    yother(logical(eye(nodes))) = 0;
    
    Im = I;
    Km = K;
    
    dEdxmstep = Km.*((xm-xother)-((Im.*(xm-xother))./((xm-xother).^2 + (ym - yother).^2).^.5));
    dEdxmstep(isnan(dEdxmstep))=0;
    dEdxm = sum(dEdxmstep,2);
    dEdymstep = Km.*((ym-yother)-((Im.*(ym-yother))./((xm-xother).^2 + (ym - yother).^2).^.5));
    dEdymstep(isnan(dEdymstep))=0;
    dEdym = sum(dEdymstep,2);
    deltam = sqrt((dEdxm).^2+(dEdym).^2)';
    
    count = count+1;
end

h=figure('position',[1500 500 1000 750],'Color','white');
hold
axis([-r r -r r])

nodeSize = 5;

%Paste nodes not in network at bottom
outnodepos = -r:(2*r)/(length(outnodeind)-1):r;
for j = 1:length(outnodeind)
    nodepos(outnodeind(j),:) = [outnodepos(j) -r];
end

%Draw lines between nodes
for i = 1:size(graph,1)
    initcoordmat = repmat(nodepos(i,:),nnz(graph(i,:)),1);
    line([initcoordmat(:,1) nodepos(logical(graph(i,:)),1)]',[initcoordmat(:,2) nodepos(logical(graph(i,:)),2)]','Color',[105/255 105/255 105/255]);
end
    
%Mark and color node positions
int=1;
for i = 1:nodes
    plot(nodepos(int,1),nodepos(int,2),'MarkerFaceColor',[colors(i,1) colors(i,2) colors(i,3)],'Marker','o','MarkerSize',nodeSize,'MarkerEdgeColor','k','LineWidth',1);
    int = int+1;
end