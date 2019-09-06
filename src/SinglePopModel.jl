#=
Single Population Model Julia Code for Simulation and Analysis

Some Useful translations for the single population model to biological timepoints
CBT=DLMO+7hrs
CBT=DLMO_mid+2hrs
CBT=circadian phase pi in the model
DLMO=circadian phase 5pi/12=1.309 in the model
=#

module SinglePopModel

using DifferentialEquations
using ParameterizedFunctions
using DataFrames
using ODEInterfaceDiffEq
using Sundials
include("./processedLight.jl")
#solver_method=RadauIIA5() #safe but slow option for integrating the odes
#solver_method=ARKODE(Sundials.Explicit(),etable = Sundials.CASH_KARP_6_4_5)
solver_method=CVODE_BDF()

spmodel = @ode_def_bare begin
    dR=-1.0*γ*R+K*cos(β1)/2.0*R*(1.0-R^4)+B(t)*(A1*0.5*(1.0-R^4)*cos(Ψ+βL1)+A2*0.5*R*(1.0-R^8)*cos(2.0*Ψ+βL2))
    dΨ=ω0+K/2.0*sin(β1)*(1+R^4)+B(t)*(σ-A1*0.5*(R^3+1.0/R)*sin(Ψ+βL1)-A2*0.5*(1.0+R^8)*sin(2.0*Ψ+βL2))
end ω0 K γ β1 A1 A2 βL1 βL2 σ


function setParameters(;paramsIn=nothing)
#=
Set the parameters for the system
Using the order as given in the spmodel parameterixed function above

=#
	if paramsIn==nothing
    	pfilename="optimalParams.dat"
    	global pvalues=[]
    	open(pfilename) do filep
        	s=read(filep, String)
        	push!(pvalues, map(x->parse(Float64,x),split(s)))
    	end
    	pvalues=pvalues[1]
    	pvalues=pvalues[1:9] #non light parameters only
	else
		global pvalues=paramsIn
	end
end #setParameters

function setB(LightIn)
	#=
	 Pass in a function of time L(t) to set this for the system. Will find B(t) and
	 save Light(t) as well.
	 =#
	global B
	global Light
	processedLight.setLight(LightIn)
	B=processedLight.getBfunc()
	Light=LightIn
end

function integrateModel(tend, initial)
    prob=ODEProblem(spmodel,initial,(0.0,tend), pvalues)
    sol=solve(prob, solver_method)
    return sol
end

function integrateSegment(initial, tstart, tend)
	#=
	Integrate a segment of the model given a start condition and tstart and tend
	=#
    prob=ODEProblem(spmodel,initial,(tstart,tend), pvalues)
    sol=solve(prob, solver_method)
    return sol
end

function integrateTransients(;numdays=50)
    tend=numdays*24.0
    prob=ODEProblem(spmodel,[0.75,0.0],(0.0,tend), pvalues)
    last=solve(prob,solver_method).u[end]
    return last
end

function getTS(tend, initial; dt=0.1)
	#Return a time series data frame for the system, used for creating actogram plots
	sol=integrateModel(tend, initial)
	ts=collect(0.0:dt:tend)
	Light_L=map(Light, ts)
	PhaseA=map(t->sol(t)[2], ts)
	RA=map(t->sol(t)[1], ts)
	df=DataFrame(Time=ts, Light_Level=Light_L, Phase=PhaseA, R=RA )
	return(df)
end

end #module single pop model
