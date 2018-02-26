
% If I want to save current results as well. Actually be careful about this
% though
% % % % save('PerturbationResults/bad.mat','E','Epert','V','Vpert','Vnorms','epsilon','K','Kpert')
%%%%

inLoc = '/home/gsiyer/Schoolwork/Chanussot/Matlab/MBO_Matlab/PerturbationResults/';
outLoc = '/home/gsiyer/Schoolwork/Chanussot/MeetingNotes/2018-01-18/PerturbationAnalysis/';
inFile = 'bad.mat';
outFile = 'Bad/';
numEvecs = 20;

A = load(strcat(inLoc,inFile));

imwrite(rescaleIm(A.K)*64,jet,strcat(outLoc,outFile,'normalClassification.png'))
imwrite(rescaleIm(A.Kpert)*64,jet,strcat(outLoc,outFile,'pertClassification.png'))

for i = 1:numEvecs
    if(i < 10)
        tempstr = strcat('evec0',num2str(i));
    else
        tempstr = strcat('evec',num2str(i));
    end
    imwrite(rescaleIm(reshape(A.V(:,i),[size(A.K,1) size(A.K,2)]))*64,jet,strcat(outLoc,outFile,tempstr,'.png'))
    imwrite(rescaleIm(reshape(A.Vpert(:,i),[size(A.K,1) size(A.K,2)]))*64,jet,strcat(outLoc,outFile,tempstr,'pert.png'))
end

