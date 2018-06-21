displaygt.m is a Matlab routine for displaying
the ground-truth symmetry line on top of the image.
It provides a practical demo on interpreting
the actuall database (in folders S and M).

The folder S (respectively, M) contains single-symmetry
(respectively, multiple-symmetry) images and ground-truth data.

Each .mat file, when loaded into Matlab,
turns into a cell array called 'segments':
segments{1} contains information about the first symmetry axis,
segments{2} the second, and so on.
Use length(segments) to find out how many ground-truth segments there are.
Naturally, length(segments) == 1 for the images in folder S.

segments{i}, for any viable index i, is a 2x2 matrix (let's call it C)
containing the coordinates of the end-points of a symmetry line.
C(1,:) is the first point (let's call it p),
C(2,:) is the second point (let's call it q).

p(1) is the 'column' coordinate (or the x coordinate, in Matlab's coordinate system),
p(2) is the 'row' coordinate (or the y coordinate, in Matlab's coordinate system).