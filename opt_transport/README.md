# opt_transport
This program will compute the optimal weight vector for the transport between two measures contained in the square [0,1] x [0,1]. The source measure has to lie on a pixel grid while the support points of the target measure can also be in arbitrary position. The masses of the pixels are defined by files containing double numbers arranged in the same grid as the pixels. See [mu.txt](samples/mu.txt) as an example of such a file for a 32 x 32 grid. An example for a target measure with points in [0,3]x[0,3] not arranged in a grid is given by [nu.txt](samples/nu.txt). Each row represents a point and consists of the three entries

        x y m
where `x` and `y` are the _non-negative_ coordinates and `m` is the mass at `(x,y)`.


## Compilation
Open a shell and switch to this folder. Enter

    cmake .

followed by

    make

to create the program `opt_transport`.


## Usage
    opt_transport <-l/-g> <source> <target> [<normalizing factor for list target>] <weights> [<transport plan>] <regularization strength>

- `<-l/-g>` - `-l` if the target measure is given as a `l`ist of points and `-g` if it is given as a pixel `g`rid. In the former situation you have to also set the normalizing factor, otherwise this argument must be omitted.
- `<source>`, `<target>` - the files containing the masses of the source and target measure as described above
- `<normalizing factor for list target>` - a number >= the largest coordinate of the target's support. All coordinates of this support will be divided by the normalizing factor to make sure it is contained in [0,1]x[0,1]. If that's already the case, simply set it to 1.
- `<weights>` - the file to which the computed optimal weight vector should be printed
- `<transport plan>` (optional) - the file to which a rough transport plan should be printed. More precisely, the file will consist of three columns where each row has the format
        s t m
and describes that according to the transport plan derived from the optimal weight vector an amount of `m` mass will be transported from somewhere on the pixel (of the source measure) with index `s` to the pixel (of the target measure) with index `t`.
- `<regularization strength>` - the factor with which the L2-summand for stabilizing the optimization is multiplied. Higher value -> more stable, but less precise. Can usually be set to 0.

To adjust the output of the program, define or undefine the preprocessor macros at the beginning of [opt_transport.cpp](opt_transport,cpp) and [Source.cpp](measures/Source.cpp).
