local class = require 'ext.class'
local Roe = require 'roe'
local Maxwell = require 'maxwell'
local MaxwellRoe = class(Roe)
MaxwellRoe.equation = Maxwell()
return MaxwellRoe
