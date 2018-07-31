dataset = [];
labels = [];

mu1 = [-2,42];
sigma1 = [12, 4 ; 4 , 8];
X = mvnrnd(mu1,sigma1,300);
dataset = [dataset; X];
labels = [labels; ones(size(X, 1), 1)];


mu1 = [23,42];
sigma1 = [12, 4 ; 4 , 8];
X = mvnrnd(mu1,sigma1,300);
dataset = [dataset; X];
labels = [labels; 2 .* ones(size(X, 1), 1)];

pert1 = size(dataset, 1);

mu1 = [-25,0];
sigma1 = [12, 4 ; 4 , 8];
X = mvnrnd(mu1,sigma1,200);
dataset = [dataset; X];
labels = [labels; 3 .* ones(size(X, 1), 1)];

mu1 = [-7,-6];
sigma1 = [12, 4 ; 4 , 8];
X = mvnrnd(mu1,sigma1,200);
dataset = [dataset; X];
labels = [labels; 4 .* ones(size(X, 1), 1)];

mu1 = [10,-6];
sigma1 = [12, 4 ; 4 , 8];
X = mvnrnd(mu1,sigma1,200);
dataset = [dataset; X];
labels = [labels; 5 .* ones(size(X, 1), 1)];

mu1 = [27,-6];
sigma1 = [12, 4 ; 4 , 8];
X = mvnrnd(mu1,sigma1,200);
dataset = [dataset; X];
labels = [labels; 6 .* ones(size(X, 1), 1)];

mu1 = [50,0];
sigma1 = [12, 4 ; 4 , 8];
X = mvnrnd(mu1,sigma1,200);
dataset = [dataset; X];
labels = [labels; 7 .* ones(size(X, 1), 1)];
dataset = [dataset, labels];

perturbation_point = pert1 + 1;


a= dataset(perturbation_point:end, :) ;
[m,n] = size(a);
idx = randperm(m) ;
b(:,:) = a(idx,:);
dataset(perturbation_point:end, :) = b;

labels = dataset(:, 3);
dataset = dataset(:, 1:2);


% plot(dataset(:,1), dataset(:,2), '.')
% axis equal
% hold on
% for i = 1:1000
%     plot(dataset(i,1), dataset(i,2), '.')
% end
clear a b i idx m mu1 n sigma1 X pert1 perturbation_point