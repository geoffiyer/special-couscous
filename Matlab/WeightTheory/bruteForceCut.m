function [bestAssign] = bruteForceCut(W, method, verbose)

if(nargin < 3)
    verbose = false;
end

n = size(W,1);

bestAssign = ones(n,1);
bestCut = Inf;

if(strcmp(method,'ncut'))
    for i=0:2:(2^n-1)
        assign = dec2binvec(i,n);
        weightsCut = sum(sum(W(assign,~assign)));
        volA = sum(sum(W(assign,:)));
        volB = sum(sum(W(~assign,:)));
        testCut = weightsCut * (1/volA + 1/volB);
        if(testCut < bestCut)
            bestAssign = assign;
            bestCut = testCut;
        end
        if(verbose)
            disp(assign)
            disp(testCut)
        end
    end    
elseif(strcmp(method,'ratiocut'))
    for i=0:2:(2^n-1)
        assign = dec2binvec(i,n);
        weightsCut = sum(sum(W(assign,~assign)));
        volA = sum(assign);
        volB = sum(~assign);
        testCut = weightsCut * (1/volA + 1/volB);
        if(testCut < bestCut)
            bestAssign = assign;
            bestCut = testCut;
        end
        if(verbose)
            disp(assign)
            disp(testCut)
        end
    end 
else
    for i=1:2:(2^n-2)
        assign = dec2binvec(i,n);
        weightsCut = sum(sum(W(assign,~assign)));
        testCut = weightsCut;
        if(testCut < bestCut)
            bestAssign = assign;
            bestCut = testCut;
        end
        if(verbose)
            disp(assign)
            disp(testCut)
        end
    end
end

end
