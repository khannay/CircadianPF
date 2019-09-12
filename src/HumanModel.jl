#=
Single Population Model Julia Code for Simulation and Analysis

Some Useful translations for the single population model to biological timepoints
CBT=DLMO+7hrs
CBT=DLMO_mid+2hrs
CBT=circadian phase pi in the model
DLMO=circadian phase 5pi/12=1.309 in the model
=#
module HumanModel

using Parameters
using DifferentialEquations
using ParameterizedFunctions
using DataFrames
using ODEInterfaceDiffEq
using Sundials
#solver_method=RadauIIA5() #safe but slow option for integrating the odes
#solver_method=ARKODE(Sundials.Explicit(),etable = Sundials.CASH_KARP_6_4_5)
solver_method=CVODE_BDF()

@with_kw struct sp_parameters
	ω0::Float64=0.263524;
	K::Float64=0.0635842;
	γ::Float64=0.024;
	β1::Float64=-0.0931817;
	A1::Float64=0.3855;
	A2::Float64=0.195123;
	βL1::Float64=-0.0026;
	βL2::Float64=-0.957756;
	σ::Float64=0.0400692;
	δ::Float64=0.0075;
	I0::Float64=9325.0;
	pp::Float64=1.5;
	α0const::Float64=0.05;
	G::Float64=33.75
end

export sp_parameters

function spmodelRHS(du,u,p,t,L)
	@unpack ω0,K,γ, β1, A1, A2, βL1, βL2, σ, δ, I0, pp, α0const,G=p
	α0var=α0const*abs(L(t))^pp/(abs(L(t))^pp+I0)
	R=u[1]
	Ψ=u[2]
	n=u[3]
	B=G*(1.0-n)*α0var
	du[1]=-1.0*γ*R+K*cos(β1)/2.0*R*(1.0-R^4)+B*(A1*0.5*(1.0-R^4)*cos(Ψ+βL1)+A2*0.5*R*(1.0-R^8)*cos(2.0*Ψ+βL2))
    du[2]=ω0+K/2.0*sin(β1)*(1+R^4)+B*(σ-A1*0.5*(R^3+1.0/R)*sin(Ψ+βL1)-A2*0.5*(1.0+R^8)*sin(2.0*Ψ+βL2))
    du[3]=60.0*(α0var*(1.0-n)-δ*n)
end


function integrateModel(tend, initial, L, p)
	spmodelThis(du,u,p,t)=spmodelRHS(du,u,p,t,L)
    prob=ODEProblem(spmodelThis,initial,(0.0,tend), p)
    sol=solve(prob, solver_method)
    return sol
end

function integrateSegment(initial, tstart, tend, L, p)
	#=
	Integrate a segment of the model given a start condition and tstart and tend
	=#
	spmodelThis(du,u,p,t)=spmodelRHS(du,u,p,t,L)
    prob=ODEProblem(spmodelThis,initial,(tstart,tend), p)
    sol=solve(prob, solver_method)
    return sol
end

function integrateTransients(;numdays=50, L, p)
	spmodelThis(du,u,p,t)=spmodelRHS(du,u,p,t,L)
    tend=numdays*24.0
    prob=ODEProblem(spmodelThis,[0.75,0.0],(0.0,tend), p)
    last=solve(prob,solver_method).u[end]
    return last
end

function getTS(tend, initial, L, p; dt=0.1)
	#Return a time series data frame for the system, used for creating actogram plots
	spmodelThis(du,u,p,t)=spmodelRHS(du,u,p,t,L)
	sol=integrateModel(tend, initial,L,p)
	ts=collect(0.0:dt:tend)
	Light_L=map(Light, ts)
	PhaseA=map(t->sol(t)[2], ts)
	RA=map(t->sol(t)[1], ts)
	df=DataFrame(Time=ts, Light_Level=Light_L, Phase=PhaseA, R=RA )
	return(df)
end

end #module HumanModel
