# Visualization
Drawing a Voronoi diagram with this software amounts to doing two things:
1. Generating the necessary data.
2. Plotting the diagram with MATLAB.


## Compilation
Open a shell and switch to this folder. Enter

    cmake .

followed by

    make

to create the program `create_diagram`.


## Usage

### Generating the Data
`create_diagram` will create two files:
- `<sites file>` contains the centers of the Voronoi cells and the corresponding weights.
- `<intersections file>` contains the points at which the hyperbolas that form the boundary of the Voronoi cells intersect.

There are three ways of creating Voronoi diagrams: with random points and weights or with points on a grid and prescribed weights or with completely arbitrary points and weights.

#### Random
    create_diagram <number of sites> <width> <height> <sites file> <intersections file>

- `<number of sites>` - the number of sites the diagram should contain
- `<width> <height>` - the sites will be contained in the rectangle [0,width] x [0,height]

#### Grid
    create_diagram <weight file> <sites file> <intersections file>

- `<weight file>` - a file containing the weights of the sites. It is assumed that the sites lie on a pixel grid contained in the square [0,1] x [0,1] which is given by the arrangement of the values in `<weight file>`. The file [weights.txt](samples/weights.txt) gives an example for such a file for 18 points in the rectangle [0,0.5] x [0,1].

#### Arbitrary
    create_diagram <sites file> <intersections file>

- `<sites file>` - this time not the output but the input file. It must contain three columns where each column describes one site: `<x-coordinate> <y-coordinate> <weight>`. We again provide a sample file: [sites.txt](samples/sites.txt).


### Plotting the Diagram
Once the sites and the intersections file have been created, the diagram is ready for being plotted. Open [plot_voronoi_diagram.m](plot_voronoi_diagram.m) in MATLAB and adjust it according to the instructions in the script. When being run, it should plot the desired Voronoi diagram.
