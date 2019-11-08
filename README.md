This repository contains python scripts that models the behavior of a FODO Cell.
The folder title gives the type of variation used.

-"thin_lens" Folder
The scripts here model the FODO cell using the thin lens approximation. The
lattice elements are represented by the respective matrices of the drift, focusing lens,
and defocusing lens components.

-"continuous_lens" Folder
In the thin lens approximation the particle positions are updated before and after the
focusing and defocusing lens. But, there is no information inside the lens. This leads
to choppy regions on the x-x' plots. The continuous program remedies this by updating the
different matrices with a single transfer matrix.

-"thick_lens" Folder
The python script here models the FODO cell with a thickness given to the focusing elements. However, when the particle is in the lens the program does not evaluate the
behavior which is the reason for the choppy sections in the x-x' plots. 
