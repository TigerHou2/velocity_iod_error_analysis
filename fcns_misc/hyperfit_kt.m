function [x,y,R] = hyperfit_kt(XY)
%--------------------------------------------------------------------------
%  
%     Algebraic circle fit with "hyperaccuracy" (with zero 2nd order bias)
%     _kt stands for Kanatani, the author of the paper that improved hyperfit
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: [x y R] is the fitting circle: center (a,b) and radius R
%
%--------------------------------------------------------------------------

centroid = mean(XY);   % the centroid of the data set

X = XY(:,1) - centroid(1);  %  centering data
Y = XY(:,2) - centroid(2);  %  centering data
Z = X.*X + Y.*Y;
ZXY1 = [Z X Y ones(length(Z),1)];
[~,S,V]=svd(ZXY1,0);
if (size(S,1) < 4 || S(4,4)/S(1,1) < 1e-12)   %  singular case
    A = V(:,4);
else                         %  regular case
    R = mean(ZXY1);
    N = [8*R(1) 4*R(2) 4*R(3) 2; 4*R(2) 1 0 0; 4*R(3) 0 1 0; 2 0 0 0];
    
    % Kanatani improvement ---------------------------
    term = zeros(size(N));
    n = size(XY,1);
    M = zeros(4,4);
    
    for i = 1:n
        M = M + ZXY1(i,:)'*ZXY1(i,:)/n;
    end
    Minv = pinv(M,max(eig(M))*0.99);
    
    for i = 1:n
        ea = ZXY1(i,:)';
        V0ea = [4*Z(i), 2*X(i), 2*Y(i), 0; ...
                2*X(i), 1, 0, 0; ...
                2*Y(i), 0, 1, 0; ...
                0, 0, 0, 0];
        A = V0ea * Minv * (ea * ea');
        term = term + ( trace(Minv*V0ea) * (ea * ea') ...
                      + ea'*(Minv*ea) * V0ea ...
                      + A + A' ) / n^2;
    end
    
    N = N - term;
    % ------------------------------------------------
    
    W = V*S*V';
    [E,D] = eig(W*inv(N)*W);
    [~,ID] = sort(diag(D));
    Astar = E(:,ID(2));
    A = W\Astar;
end

Par = [ -(A(2:3))'/A(1)/2+centroid , ...
        sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];
x = Par(1);
y = Par(2);
R = Par(3);

end

