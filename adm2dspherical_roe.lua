local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'roe'
local ADM2DSpherical = require 'adm2dspherical'
local ADM2DSphericalRoe = class(Roe)
function ADM2DSphericalRoe:init(args)
	args = table(args)
	args.equation = ADM2DSpherical(args)
	ADM2DSphericalRoe.super.init(self, args)
end
return ADM2DSphericalRoe
