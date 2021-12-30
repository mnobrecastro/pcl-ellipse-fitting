% Miguel Nobre Castro
% mnobrecastro@gmail.com
% Created on Feb 17th, 2021
%
% Implementation of "Direct Least Square Fitting of Ellipses" (Fitzgibbon et al., 1999)

function [] = main()
clc; clear all;

% Generate dataset
N = 50;
model = [2.0, 1.0, 1.0, 1.0, pi/10]';
th = linspace(0, pi/3, N)';
% maxang = pi/6;
% th = linspace(maxang/N,maxang,N)';
pts = make_ellipse(model, th, 0.001); %0.1 leads to problems

% % % Dummy dataset
% % pts = [ 1.6823, -0.2154, -0.8957, 0.3129, 2.2218, 2.9022;
% %         2.1258, 1.5107, 0.3741, -0.1252, 0.4890, 1.6246;
% %       ]';

% Direct Least Square Fitting of Ellipses (Fitzgibbon et al., 1999)
idx = 1:N;
idx = idx(randperm(length(idx)));
n = 6; % Minimum number of seed points
con = fit_ellipse( pts(idx(1:n),1), pts(idx(1:n),2));

% Formulation transformation Conic -> Parametric (plot purposes)
par = real(conic2parametric(con));
th = linspace(0,2*pi)';
elp = make_ellipse(par, th, 0);

% Plot the Ellipse model
figure(1);
plot(pts(:,1), pts(:,2), '*');  % Draw dataset
hold on;
plot(pts(idx(1:n),1), pts(idx(1:n),2), 'ro');  %Draw seed pts
hold on;
plot(elp(:,1), elp(:,2), 'k'); % Draw Ellipse model
hold off;
% axis equal;
% axis([0, 3, 0, 3]);

disp(model);
disp(par);

end

function con = fit_ellipse(x,y)
    % "Direct Least Square Fitting of Ellipses" (Fitzgibbon et al., 1999)
    %
    % This function approximates an ellipse to a set of (x,y) points.

    % Build design matrix D
    D = [ x.*x, x.*y, y.*y, x, y, ones(size(x)) ]
    % Build scatter matrix S
    S = D'*D
    % Build 6x6 contraint matrix C
    C = zeros(6,6); C(1,3)=-2; C(2,2)=1; C(3,1)=-2
    % Solve generalized eigensystem: S*a = lambda*C*a
    [gevec, geval] = eig(S,C)
    % Find the negative eigenvalue
    [NegR, NegC] = find(geval<0 & ~isinf(geval))
    % Get fitted Conic eq. params a
    con = gevec(:,NegC);
end

function par = conic2parametric(con)
    % Converts Conic parameters to Parametric parameters

    % Conic params
    c_A=con(1); c_B=con(2); c_C=con(3); c_D=con(4); c_E=con(5); c_F=con(6);
    % Build aux 3x3 matrix M0    
    M0 = [  c_F, c_D/2, c_E/2;
            c_D/2, c_A, c_B/2;
            c_E/2, c_B/2, c_C ];
    % Build aux 2x2 matrix M
    M = [   c_A, c_B/2;
            c_B/2, c_C ];
    % Calculate the eigenvalues and eigenvectors of matrix M
    [eigvec, eigval] = eig(M);
    eigval = [eigval(1,1), eigval(2,2)];
    % Order the eigenvalues so that |lambda_1 - c_A| <= |lambda_1 - c_C|
    if (abs(eigval(1)-c_A) > abs(eigval(1)-c_C))
        aux = eigval(1);
        eigval(1) = eigval(2);
        eigval(2) = aux;
    end
    
    % Parametric eq. params
    p_a = sqrt( -det(M0) / (det(M)*eigval(1)) );
    p_b = sqrt( -det(M0) / (det(M)*eigval(2)) );
    p_h = (c_B*c_E-2*c_C*c_D) / (4*c_A*c_C-c_B^2);
    p_k = (c_B*c_D-2*c_A*c_E) / (4*c_A*c_C-c_B^2);
% %     p_t = acot( (c_A-c_C)/c_B )/2;
    p_t = (pi/2 - atan((c_A-c_C)/c_B))/2;
    par = [p_a, p_b, p_h, p_k, p_t]';
end

function con = parametric2conic(par)
    % Converts Parametric parameters to Conic parameters

    % Parametric eq. params
    p_a = par(1); p_b = par(2); p_h = par(3); p_k = par(4); p_t = par(5);
    % Conic eq. params
    c_A = (p_b*cos(p_t))^2 + (p_a*sin(p_t))^2;
    c_B = -2*cos(p_t)*sin(p_t) * (p_a^2-p_b^2);
    c_C = (p_b*sin(p_t))^2 + (p_a*cos(p_t))^2;
    c_D = -2*c_A*p_h - p_k*c_B;
    c_E = -2*c_C*p_k - p_h*c_B;
    c_F = -(p_a*p_b)^2 + c_A*p_h^2 + c_B*p_h*p_k + c_C*p_k^2;
    con = [c_A, c_B, c_C, c_D, c_E, c_F]';
end

function pts = make_ellipse(par, th, noise)
    % Generates the points of an ellipse using it Parametric 'par'
    % equation parameters, with a given rotation 'th' and noise.

    % Parametric eq. params
    p_a=par(1); p_b=par(2); p_h=par(3); p_k=par(4); p_t=par(5);
    % Calculate Ellipse data points
    pts = zeros(length(th),2);
    pts(:,1) = p_h + cos(p_t)* p_a*cos(th) - sin(p_t)* p_b*sin(th);
    pts(:,2) = p_k + sin(p_t)* p_a*cos(th) + cos(p_t)* p_b*sin(th);
    % Introduce noise in the data
    pts(:,1) = pts(:,1) + rand(size(pts(:,1)))*2*noise - noise;
    pts(:,2) = pts(:,2) + rand(size(pts(:,1)))*2*noise - noise;
end

