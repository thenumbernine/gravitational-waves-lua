local class = require 'ext.class'
local table = require 'ext.table'
local ADM1D5Var = require 'adm1d5var'
local Roe = require 'roe'

local ADM1D5VarRoe = class(Roe)

ADM1D5VarRoe.name = 'ADM 1D 5-Var Roe'

function ADM1D5VarRoe:init(args)
	args = table(args)
	args.equation = ADM1D5Var(args)
	ADM1D5VarRoe.super.init(self, args)
end

return ADM1D5VarRoe
