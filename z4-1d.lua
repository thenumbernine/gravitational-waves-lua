--[[

http://arxiv.org/pdf/1106.2254v2.pdf

alpha = lapse
beta^i = shift, ignored for now
B^i is something ...
gammaTilde_ij = conformal spatial metric ... det = 1


alpha,t = -2 alpha (K - 2 Theta)
-- beta is ignored, no shift
B^x,t = connTilde^x,t - eta B^x
gammaTilde_xx,t = -2 alpha ATilde_ij^TF
ATilde_ij,t = 

ATilde_ij^TF = trace-free ATilde_ij
... = ATilde_ij - 1/3 gammaTilde_ij ATilde^k_k 
	(in a 1D sim do we use 1/3 for simulating 1D of a 3D sim, or do we use 1/1 to simulate a 1D spatial dimension?)
1D of 3D option: ATilde_xx - 1/3 ATilde_xx / gammaTilde_xx = (1 - 1/(3 gammaTilde_xx)) ATilde_xx
1D of 1D option: ATilde_xx - 1/1 ATilde_xx / gammaTilde_xx = (1 - 1 / gammaTilde_xx) ATilde_xx
(but isn't A_ij defined as the trace-free portion of the extrinsic curvature?)
(and ATilde_ij is the conformal metric version of A_ij)
I'm just going to ignore that ^TF, because ATilde_ij should be ^TF already, and because in BSSN it's not there (as it is missing in all other ATilde_ij references in th Z4)

... and this isn't deconstructed into flux matrix, so I'm going to pass for now

--]]
