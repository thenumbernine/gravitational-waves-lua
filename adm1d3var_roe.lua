local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'roe'
local ADM1D3VarRoe = class(Roe)
local ADM1D3Var = require 'adm1d3var'
function ADM1D3VarRoe:init(args)
	args = table(args)
	args.equation = ADM1D3Var(args)
	ADM1D3VarRoe.super.init(self, args)
end
return ADM1D3VarRoe
