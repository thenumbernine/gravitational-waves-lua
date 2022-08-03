[![paypal](https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=KYWUWS86GSFGL)

Numerical simulation of hyperbolic formalisms

Made for getting some numerical relativity demos working

Requires my lua-ext, symmath-lua, and optionally lua-glapp, and the Malkia UFO ffi lua files for the OpenGL display.

There's a text-based version you can uncomment that will run it to the command-line if you don't want to bother set up the GL stuff.

`adm1d_v1.lua` is based on formalism described in Alcubierre's "Introduction to 3+1 Numerical Relativity" in the Toy 1+1 example chapter.

Specifically, a 3-variable hyperbolic system with lapse (alpha) and metric (g) separated from the state variable computations.


`adm1d5var.lua` is a hyperbolic sim based on Alcubierre's "The appearance of coordinate shocks in hyperbolic formalisms of General Relativity"

found [here](http://arxiv.org/pdf/gr-qc/9609015v2.pdf)


`adm2dspherical.lua` is the spherical solution in the same paper mentioned above.

Figured out that the paper's, "eigenfields" is an inner product with the variable basis and the rows that make up the inverse eigenvector matrix.


`adm3d.lua` is my start on a 3D implementation according to the gauge shock paper.

Maybe I'll also use Alcubierre's "Introduction to Numerical Relativity" paper found [here](http://cgwa.phys.utb.edu/Files/Events/29_610_Alcubierre_numerical.pdf)


`euler1d.lua` is the Euler hydrodynamic equations I used as a test-case to verify the Roe solver was working.  Good ol' Euler fluid equations.


`maxwell.lua` is Maxwell equations described as hyperbolic equations in section 4.3 of Trangenstein's "Numerical Solutions of Hyperbolic Partial Differential Equations"


`bssnok1d.lua` is a hodge-podge of the BSSNOK partial wrt time and the hyperbolic formalism of BSSNOK described in Alcubierre's Hyperbolicity chapter.


`mhd.lua` is based on "ATHENA: A NEW CODE FOR ASTROPHYSICAL MHD" 2008 by Stone, Gardiner, Teuben, Hawley, and Simon


and then I started adding attempts at implicit solvers.

the backward-euler conj.res. implicit and the backward-euler newton+conj.res. implicit don't work so great at the moment.

the roe+backward-euler conj.res. implicit is stable for a while, but with some errors at boundary conditions (which probably what makes it explode... I hope)

