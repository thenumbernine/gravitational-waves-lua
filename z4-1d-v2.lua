--[[
same as Z4 except using the 1D ADM's 
	D_g = (ln gamma_xx),x = gamma_xx,x / gamma_xx
	and KTilde_xx = sqrt(gamma_xx) K_xx

alpha,t = -f alpha^2 (tr K - m Theta)
	= -f alpha^2 (K_xx / gamma_xx - m Theta)
gamma_xx,t = -2 alpha K_xx
a_x,t = - f alpha K_xx,x / gamma_xx 
		+ f alpha m Theta,x 
		- alpha a_x (f' alpha + f) (K_xx / gamma_xx - m Theta)
		+ 2 f alpha K_xx d_xxx / gamma_xx^2
d_xxx,t = -alpha K_xx,x - alpha a_x K_xx
K_xx,t = - alpha a_x,x 
		+ 2 alpha Z_x,x 
		- alpha a_x^2 
		+ alpha d_xxx a_x / gamma_xx
		- 2 alpha d_xxx Z_x / gamma_xx
		- alpha K_xx^2 / gamma_xx
		- 2 alpha Theta K_xx
		- alpha S_xx
		+ 1/2 alpha (S_xx / gamma_xx - tau) gamma_xx
Theta,t = alpha Z_x,x / gamma_xx
		- alpha Z_x d_xxx / gamma_xx^2
		- alpha Theta K_xx / gamma_xx
		- alpha tau
		- alpha a_x Z_x / gamma_xx
Z_x,t = alpha Theta,x
		- 2 alpha Z_x K_xx / gamma_xx
		- alpha S_x
		- alpha Theta a_x

so D_g,t = (ln gamma_xx),xt
	= (gamma_xx,x / gamma_xx),t
	= (2 d_xxx / gamma_xx),t
	= 2 d_xxx,t / gamma_xx - 2 d_xxx gamma_xx,t / gamma_xx^2
	= 2 (-alpha K_xx,x - alpha a_x K_xx) / gamma_xx
		- 2 d_xxx (-2 alpha K_xx) / gamma_xx^2
	= -2 alpha K_xx,x / gamma_xx 
		- 2 alpha a_x K_xx / gamma_xx
		+ 4 alpha d_xxx K_xx / gamma_xx^2

and KTilde_xx,t = (sqrt(gamma_xx) K_xx),t
	= gamma_xx,t K_xx / (2 sqrt(gamma_xx)) + sqrt(gamma_xx) K_xx,t
	= (-2 alpha K_xx) K_xx / (2 sqrt(gamma_xx)) + sqrt(gamma_xx) (
	- alpha a_x,x 
		+ 2 alpha Z_x,x 
		- alpha a_x^2 
		+ alpha d_xxx a_x / gamma_xx
		- 2 alpha d_xxx Z_x / gamma_xx
		- alpha K_xx^2 / gamma_xx
		- 2 alpha Theta K_xx
		- alpha S_xx
		+ 1/2 alpha (S_xx / gamma_xx - tau) gamma_xx
	)
--]]
