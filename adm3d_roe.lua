local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'roe'
local ADM3D = require 'adm3d'
local ADM3DRoe = class(Roe)
ADM3DRoe.name = 'ADM 3D Roe'
function ADM3DRoe:init(args)
	args = table(args)
	args.equation = ADM3D(args)
	ADM3DRoe.super.init(self, args)
end
return ADM3DRoe
