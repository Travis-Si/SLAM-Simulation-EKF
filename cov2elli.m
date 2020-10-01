function [X,Y] = cov2elli(x,P,ns,NP)
    % Ellipsoidal representation of multivariate Gaussian variables (2D). Different
    % sigma-value ellipses can be defined for the same covariances matrix. The most useful
    % ones are 2-sigma and 3-sigma
    %Ellipse points from mean and covariances matrix.
    %   [X,Y] = COV2ELLI(X0,P,NS,NP) returns X and Y coordinates of the NP
    %   points of the the NS-sigma bound ellipse of the Gaussian defined by
    %   mean X0 and covariances matrix P.
    %
    %   The ellipse can be plotted in a 2D graphic by just creating a line
    %   with line(X,Y).
    persistent circle
    if isempty(circle)
        alpha = 2*pi/NP*(0:NP);
        circle = [cos(alpha);sin(alpha)];
    end
     [R,D]=svd(P);
     d = sqrt(D);
    %circle -> aligned ellipse -> rotated ellipse -> ns-ellipse
     ellip = ns*R*d*circle;
    X = x(1)+ellip(1,:);
    Y = x(2)+ellip(2,:);
end