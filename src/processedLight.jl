module processedLight

using DifferentialEquations
using ParameterizedFunctions
using ODEInterfaceDiffEq

function setLight(L)
	# Pass in a function of time L(t) to set this for the system
	global Light
	Light=L
end

lightTransformSystem=@ode_def_bare begin
  α0(t)=αconst*abs(Light(t))^p/(abs(Light(t))^p+I0)
  dn=60.0*(α0(t)*(1.0-n)-δ*n)
end αconst p I0 δ

function getBfunc(;repeatDays=500)
	#=
		Solves the nonlinear ode using the function L(t) to give the
		transformed light B(t) which is actually presented to the SCN

		getBfunc()

		Assumes that the light function is given in a way that t in 0 Infinity is allowed....
		For data you will need to periodically extend the data

		Returns the interpolation function for B(t) the processed light input
		presented to the SCN
	=#
	lightmax=24.0*repeatDays
	lightParams=[0.05,1.5,9325.0, 0.0075]
	#Estimate the initial conditions
	problight = ODEProblem(lightTransformSystem,[0.01],(0.0,lightmax*2),lightParams)
	solTrans=solve(problight, RadauIIA5())[end]
	problight = ODEProblem(lightTransformSystem,solTrans,(0.0,lightmax),lightParams)
	solLight=solve(problight, RadauIIA5())
	#Changed this to not do rem in solLight so no automatic wrapping
	B(t)=33.75*(1.0-solLight(t)[1])*0.05*abs(Light(t))^1.5/(abs(Light(t))^1.5+9325.0)
	return B
end
end #module processedLight
