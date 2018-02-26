
% Idea: Have two data sets, each with n points (n small for now). Assume
% that there is a bijective correspondence between data sets
%     ex: points are sampled from the same scene, same time, two
%         different cameras. Match points based on time of sampling
% Also assume that we've already determined some correspondences (m << n).
% For now try to determine the full set of correspondences by looking at
% all possible permutations and minimizing an energy function.

GenerateData;  % X_G contains 9 points. 3 classes, 3 points per class.
               % X_1, X_2 are projections into 2D

E = inf;
final_perm = [4 5 6 7 8 9];
for p = perms([4 5 6 7 8 9])'
    tempE = permEnergy(p, X_1, X_2);
    if( tempE < E)
        final_perm = p;
        E = tempE;
    end
end