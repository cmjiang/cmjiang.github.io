close all;clear all;
%% Generate random data set
ns = 50;
offset = 10;

S1 = rand(ns,2);
S1(:,1) = S1(:,1)*5;
S1(:,2) = S1(:,2)*10;

S2 = rand(ns,2);
S2(:,1) = S2(:,1)*5+offset;
S2(:,2) = S2(:,2)*10;

S3 = rand(ns,2);
S3(:,1) = S3(:,1)*5+offset*2;
S3(:,2) = S3(:,2)*10;

S4 = rand(ns,2);
S4(:,1) = S4(:,1)*5+offset*3;
S4(:,2) = S4(:,2)*10;

S = [S1;S2;S3;S4];
n = size(S,1);

%% Compute the weight matrix W
sigma = 2;
W = zeros(n,n);
for i = 1:n
    for j = 1:n
        W(i,j) = exp(-sqrt(sum((S(i,:)-S(j,:)).^2))/(2*sigma^2));
    end
end

%% Compute the Laplacian
d = sum(W,2);
D = diag(d);
L = diag(d)-W;

%% Solve the eigenvalue problem
[X,lambda] = eig(L,D);
lambdad = diag(lambda);

%% Cluster the eigvectors by K-means
numEigenvalues = 4;
Y = zeros(n,numEigenvalues);
for i = 1:n
    Y(i,:) = X(i,1:numEigenvalues)/...
        norm(X(i,1:numEigenvalues));
end

%% Try different numbers of clusters and plot
numClusters1 = 2;
c1 = kmeans(Y(:,1:numClusters1),numClusters1);
numClusters2 = 3;
c2 = kmeans(Y(:,1:numClusters2),numClusters2);
numClusters3 = 4;
c3 = kmeans(Y(:,1:numClusters3),numClusters3);
area = 25;
area2 = 10;

figure(1);
scatter(S(:,1),S(:,2),area,'filled');
%title('Original Data','FontSize', 15);

figure(2);
for i = 1:4
    subplot(2,2,i);
    scatter((1:n)',X(:,i),area2,'filled');
    title(strcat('x_',int2str(i)),'FontSize', 15);
end

% 2-way clustering
figure(3);
for i = 1:n
    switch c1(i)
        case 1
            scatter(S(i,1),S(i,2),area,'r','filled');
        case 2 
            scatter(S(i,1),S(i,2),area,'b','filled');
        case 3 
            scatter(S(i,1),S(i,2),area,'g','filled');
        case 4 
            scatter(S(i,1),S(i,2),area,'k','filled'); 
    end
    hold on;
end
%title('2-way Partitioning','FontSize', 15);

% 3-way clustering
figure(4);
for i = 1:n
    switch c2(i)
        case 1
            scatter(S(i,1),S(i,2),area,'r','filled');
        case 2 
            scatter(S(i,1),S(i,2),area,'b','filled');
        case 3 
            scatter(S(i,1),S(i,2),area,'g','filled');
        case 4 
            scatter(S(i,1),S(i,2),area,'k','filled'); 
    end
    hold on;
end
%title('3-way Partitioning','FontSize', 15);     

% 4-way clustering
figure(5);
for i = 1:n
    switch c3(i)
        case 1
            scatter(S(i,1),S(i,2),area,'r','filled');
        case 2 
            scatter(S(i,1),S(i,2),area,'b','filled');
        case 3 
            scatter(S(i,1),S(i,2),area,'g','filled');
        case 4 
            scatter(S(i,1),S(i,2),area,'k','filled'); 
    end
    hold on;
end
%title('4-way Partitioning','FontSize', 15);     
