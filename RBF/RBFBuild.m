function RBF_para = RBFBuild(Problem,TrainIn,TrainOut)

for i = 1 : Problem.M
    RBF_para{i} = RBFCreate(TrainIn, TrainOut(:,i), 'gaussian');
end
