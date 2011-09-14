Data Structures for Gravitational Microlensing
==============================================

For this project, we are dealing with large amounts of data.
Storing this data in an efficient way is incredibly important.

The most important data structures are for the lenses.
The lenses are composed of three floats: the x position, y position and the object's mass.
As tempting as it is to store these as a C struct, this would result in sub-optimal cache use.
For efficiency, each of the float types should have its own array.

_For a NULL lense, the mass is zero_

As an example, if we had 1000 lenses we would have three float arrays of length 1000.
The relevant attributes for the kth object would be found at x[k], y[k] and mass[k].
As they are all positioned next to each other, entire chunks will be loaded into memory simultaneously.
