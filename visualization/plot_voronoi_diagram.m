% PLOT_VORONOI_DIAGRAM Plot the Voronoi diagram for a given set of sites.
%
%   Author: Valentin Hartmann
%
%   This script plots the Voronoi diagram of the sites specified in the
%   file sitesPath with the intersetion points of the hyperbolas specified
%   in intersectionsPath. These two files are created by the C++-program
%   'create_diagram'.
%   Adjust the following variables and (un)comment code if you e.g. only
%   want to plot the diagram and the points without the circles around them
%   representing their weight.
%   We use the same notation for the variables as in the bachelor thesis.

% the file containing the sites
sitesPath = 'sites.txt';
% the file containing the intersection points of the hyperbolas forming the
% boundaries of the Voronoi cells
intersectionsPath = 'intersections.txt';

% the number of points used to linearly approximate one hyperbola curve
% WARNING: A high resolution in combination with a large number of points
% can let your computer run out of memory.
resolution = 100000;

% the excerpt of the diagram in R^2 which should be plotted
% realArea = [0 1.6138403415682 0 0.861755371095];
realArea = [0 1 0 1];

% the size of the plot on the screen in pixels in the larger dimension
sizeOnScreen = 790;

% lineWidth = 1;
lineWidth = 0.1;

% END OF VARIABLES THAT NEED TO BE ADJUSTED




realAreaX = realArea(2) - realArea(1);
realAreaY = realArea(4) - realArea(3);
% the area on the screen in which the diagram is plotted (ratio must match
% that of realArea)
canvasArea = 10*ones(4, 1);

if (realAreaX > realAreaY)
    canvasArea(3) = canvasArea(3) + sizeOnScreen;
    canvasArea(4) = canvasArea(4) + sizeOnScreen*realAreaY/realAreaX;
else
    canvasArea(3) = canvasArea(3) + sizeOnScreen*realAreaX/realAreaY;
    canvasArea(4) = canvasArea(4) + sizeOnScreen;
end

sites = importdata(sitesPath, ' ');
intersections = importdata(intersectionsPath, ' ');

% EPS and INFINITY are used when a hyperbola is open on one or two sides.
% x -> 1/x will then e.g. be plotted in the range
% [transformed enpoint 1, INFINITY].
% Note that in the current version of the program computing the
% intersections this case does not occur, see the comments there.
EPS = 1e-5;
INFINITY = EPS^-1;

figure('position', canvasArea);
hold on;
axis(realArea);

% plot the points
plot(sites(:,1), sites(:,2), '.', 'Color', 'k', 'MarkerSize', 15);

% plot the circles representing the weights
% for i = 1 : size(sites, 1)
%     diameter = 2*sites(i,3);
%     if diameter ~= 0
%         rectangle('Position', [sites(i,1) - sites(i,3), sites(i,2) - sites(i,3), diameter, diameter], 'Curvature', [1, 1]);
%     end
% end

% plot the hyperbola segments
for i = 1 : size(intersections, 1)
    % p and q are the points between which the current hyperbola segment
    % lies.
    p = intersections(i, 1:2);
    q = intersections(i, 4:5);
    w_p = intersections(i, 3);
    w_q = intersections(i, 6);
    end_1 = intersections(i, 7:8);
    end_2 = intersections(i, 9:10);
    
    % We have a line instead of a hyperbola.
    if w_p == w_q
        x = linspace(end_1(1), end_2(1), resolution);
        y = linspace(end_1(2), end_2(2), resolution);
        plot(x, y, 'LineWidth', lineWidth);
        continue;
    end

    norm_qp = norm(q-p);
    diff_qp_1 = q(1) - p(1);
    if diff_qp_1 == 0
        sign_qp_1 = -1;
    else
        sign_qp_1 = sign(diff_qp_1);
    end

    cos_qp = abs(diff_qp_1) / norm_qp;
    sin_qp = sign_qp_1 * (q(2) - p(2)) / norm_qp;

    a = abs(w_p-w_q) / 2;
    b = sqrt(norm_qp^2/4 - a^2);

    M = (1/2) * (q + p);
    A = zeros(2, 2);
    A(1,1) = a*cos_qp + b*sin_qp;
    A(2,1) = a*sin_qp - b*cos_qp;
    A(1,2) = a*cos_qp - b*sin_qp;
    A(2,2) = a*sin_qp + b*cos_qp;
    A = (1/2) * A;

    A_inv = zeros(2, 2);
    A_inv(1,1) = cos_qp/a + sin_qp/b;
    A_inv(2,1) = cos_qp/a - sin_qp/b;
    A_inv(1,2) = sin_qp/a - cos_qp/b;
    A_inv(2,2) = sin_qp/a + cos_qp/b;


    % Does the respective endpoint exist or does the hyperbola go to
    % infinity? (i.e. our endpoint lies on the boundary of the realArea)
    infinit_1 = false;
    infinit_2 = false;
    if isnan(end_1(1))
        infinit_1 = true;
    end
    if isnan(end_2(1))
        infinit_2 = true;
    end

    % transform the endpoints to the graph of x -> 1/x or x -> -1/x
    if p(1) >= q(1)
        if w_p < w_q
            if infinit_1
                end_1(1) = EPS;
            else
                end_1 = A_inv*(end_1 - M)';
            end
            if infinit_2
                end_2(1) = INFINITY;
            else
                end_2 = A_inv*(end_2 - M)';
            end
        else
            if infinit_1
                end_1(1) = -INFINITY;
            else
                end_1 = A_inv*(end_1 - M)';
            end
            if infinit_2
                end_2(1) = -EPS;
            else
                end_2 = A_inv*(end_2 - M)';
            end
        end
    else
        if w_p > w_q
            if infinit_1
                end_1(1) = INFINITY;
            else
                end_1 = A_inv*(end_1 - M)';
            end
            if infinit_2
                end_2(1) = EPS;
            else
                end_2 = A_inv*(end_2 - M)';
            end
        else
            if infinit_1
                end_1(1) = -EPS;
            else
                end_1 = A_inv*(end_1 - M)';
            end
            if infinit_2
                end_2(1) = -INFINITY;
            else
                end_2 = A_inv*(end_2 - M)';
            end
        end
    end
    
    x = linspace(min(end_1(1), end_2(1)), max(end_1(1), end_2(1)), resolution);
    y = x.^-1;

    hyperbola = zeros(2, length(x));
    hyperbola(1, :) = x;
    hyperbola(2, :) = y;

    hyperbola = M'*ones(1, length(x)) + A*hyperbola;

    plot(hyperbola(1, :), hyperbola(2, :), 'LineWidth', lineWidth, 'Color', 'b');
end

% display a grid
% useful for diagrams stemming from images
% xResolution = 32;
% yResolution = 32;
% set(gca, 'XTick', 0:1/xResolution:1);
% set(gca, 'YTick', 0:1/yResolution:1);
% set(gca, 'XGrid', 'on');
% set(gca, 'YGrid', 'on');

% save the plot to a png file without the axes
axis off;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf, 'Color', [1 1 1]);
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperSize', [21.1506 21.1506]);
print('diagram.pdf', '-dpdf', '-r400');