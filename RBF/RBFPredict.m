function preTest = RBFPredict(para,Testin)

N     = size(Testin,1);
[~,M] = size(para);

preTest = zeros(N,M);

for i = 1 : N
    for j = 1 : M
        preTest(i,j) = RBFInterp(para{j},Testin(i,:));
    end
end

