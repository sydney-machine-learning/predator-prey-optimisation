
Findex=14;
dim=100;
[low,up]=test_functions_range(Findex);
lb=low*ones(1,dim);
ub=up*ones(1,dim);
T=1000;
rng(1,'philox');
prob = @test_functions;
[bestsol, bestfitness, BestFitIter, P, f] = QCMBO_2( prob, lb, ub, T, Findex);