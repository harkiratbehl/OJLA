function [ x, xv ] = BinIntMBest( f,A,b,Aeq,beq,options, M )
%INTMBEST Summary of this function goes here
%   Detailed explanation goes here
willplot = 0;
xdim = length(f);%num of edges
x = zeros(xdim, M);%size numofedges*100 for conditions
y = zeros(xdim, M);
xv = zeros(M,1);%size of 100*1 for values
yv = zeros(M,1);
Constraints = cell(M,1);%cell of size 100*1.. so mth m has constraints for mth solution
if(willplot)
    tstart = tic;
    tends=zeros(1,M);
end
for m=1:M
    if(m==1)
          [cx, value] = gurobi_ilp(f, A,b,Aeq,beq);
      cx = abs(round(cx));
      x(:,m) = cx;
      xv(m) = value;
    else
        [c, k] = min(yv(1:(m-1)));%min of nextbestsolutions of all solutions till now
        if(c == 1e20)
            m = m - 1;
            break;
        end
        x(:,m) = y(:,k);
        xv(m) = yv(k);%then this m is given the min of all nextbestsolutions
        diff1 = (x(:,m) ~= x(:,k));%edges which were not in orig soln but in its nextbest soln
        diff = find(diff1==1);%index of those edges
        Constraints{m} = [Constraints{k}; diff(1), x(diff(1), m)];%we want this constraint in mth solution
        Constraints{k} = [Constraints{k}; diff(1), ~x(diff(1), m)];%this edge shouldn't be there
        [value,cy] = CalcNextBestSolution(Constraints{k}, x(:,k),f, A, b, Aeq, beq,options);
        cy = abs(round(cy));
        y(:,k) = cy;
        yv(k) = value;
    end
    [value, cy] =  CalcNextBestSolution(Constraints{m}, x(:,m), f, A, b, Aeq, beq,options);
    cy = abs(round(cy));
    y(:,m) = cy;
    yv(m) = value;
    if(willplot)
        tends(m) = toc;
    end
end
M = m;
xv = xv(1:M);
x=x(:,1:M);
if(willplot)
    
    figure(99);
    plot(1:M,tends,'r+');
end
end

function [v,y] = CalcNextBestSolution(ce, xstar,f, A,b,Aeq,beq, options)
xdim = length(xstar);%numofedges
A=[A;xstar'];%add a dummy bb for that exist soln
b=[b;sum(xstar) - 1];%for that bb it's numoftrgr-1
[ch, ~] = size(ce);%number of constraints... for this trgt it's fixed
y = zeros(xdim, 1);
assigned = zeros(xdim, 1); 
v = 0;
if(~isempty(ce))
    y(ce(:,1)) = ce(:,2);
    assigned(ce(:,1)) = 1;
    
    b = b - A * y;
    beq = beq - Aeq * y;
    v = f' * y;
    
    f(ce(:,1)) = [];
    A(:,ce(:,1)) = [];
    Aeq(:,ce(:,1)) = [];
end
[y1, v1] = gurobi_ilp(f, A, b, Aeq, beq);
y1=abs(round(y1));
if(isempty(y1))
    v = 1e20;
else
    y(assigned == 0) = y1;
    v = v + v1;
end
end


