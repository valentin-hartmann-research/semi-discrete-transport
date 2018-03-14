# create_sites_file
This is a utility that takes a file A with double-valued columns in the format

        x y m
and a file B with one double-valued column

        w
and creates a copy C of A with `m` replaced by `w`. This can be used if you have a measure given by a list of points with associated masses, have computed the weight vector for the transport from another measure to this one and want to create the corresponding Voronoi diagram which requires the coordinates and weights in one file. C then contains the sites of the Voronoi diagram. Feel free to include this functionality directly in [create_diagram.cpp](../visualization/create_diagram.cpp).


## Compilation
Open a shell and switch to this folder. Enter

    g++ create_sites_file.cpp -o create_sites_file

to create the program `create_sites_file`.


## Usage
    create_sites_file <target measure> <normalizing factor> <weights> <sites>

- `<target measure>` - the file A from above containing the target measure
- `<weights>` - the file B from above containing the weights
- `<normalizing factor>` - a number >= the largest coordinate of the target's support. All coordinates of this support will be divided by the normalizing factor to make sure it is contained in [0,1]x[0,1]. If that's already the case, simply set it to 1.
- `<sites>` - the file C from above to which the Voronoi sites will be printed
