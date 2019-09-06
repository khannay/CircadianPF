module CircadianPF

using Distributions, StatsBase, DifferentialEquations
using CSV, DataFrames, LinearAlgebra, SharedArrays, Plots
include("./SinglePopModel.jl")
include("./HCHSLight.jl")
include("./LightSchedules.jl")

greet() = print("Hello World!")

function initModel(LightIn)
    SinglePopModel.setB(LightIn)
end

function makeFakeData(params, σp, σy; num_data=40)
    #=
        Generate some fake data using the specified LightSchedule
    =#
    SinglePopModel.setParameters(paramsIn=params)
    currentVal=SinglePopModel.integrateTransients()
    Psi=[]
    push!(Psi, angle(exp(im*currentVal[2])))
    tcurrent=0.0
    while tcurrent < num_data*24.0
        currentVal=systemDynamics(currentVal, params, σp, tcurrent, tcurrent+24.0)
        observedState=angle(exp(im*currentVal[2]))
        push!(Psi, observedState)
        tcurrent+=24.0
    end
    #Add the measurement noise
    PsiMeasured=Psi .+ rand(VonMises(0.0, 1.0/σy), length(Psi))
    return(Psi, PsiMeasured)
end

function systemDynamics(ustart, params, σp, tstart, tend)
    #=
    This function implements the sp model predictions with intrinsic noise
    in the dynamics. You need to give the time start and tend because the
    system is not autonomous. Also σp has two elements for the noise in each
    component of the model (R, ψ)
    =#
    SinglePopModel.setParameters(paramsIn=params)
    sol=SinglePopModel.integrateSegment(ustart, tstart, tend)
    vp=[Normal(0.0, σp[1]), VonMises(1.0/σp[2],1)]
    processNoise=[rand(vp[1],1)[1], rand(vp[2],1)[1]]
    state=sol[end] .+ processNoise #add the process noise to the system
end

function runParticleFilter(params, dataObs, σp, σy; N=1000, init=[0.70, 0.0])
    #=
        Implement a basic bootstrap particle filter
    =#

    tau=length(dataObs)
    xpf=zeros(tau+1,N,2) #Store the results for each particle here
    wloglik=log(1.0/N)

    #Allow MH to find initial conditions as well
    if init==[0.70,0.0]
        xpf[1,:,:] = rand(N,2) #choose random initial states for particles
        xpf[1,:,2]*=2π #Put the initial phases between zero and 2π
    else
        xpf[1,:,:]=init
    end
    tstart=0.0

    #Now interate through each time point
    for t in 2:(tau+1)
        tend=tstart+24.0

        w_tilde=ones(N)*1/N #Likelihood (weights unnormalized)

        # Parallel Part is not working
        Threads.@threads for p in 1:N
            xpf[t,p,:]= systemDynamics(xpf[t-1,p,:],params,σp,tstart,tend ) # 1. Importance sampling step, project each particle forward using the underlying dynamics
            observedState=angle(exp(im*xpf[t,p,2])) #Find the state given by the simulation
            w_tilde[p]= pdf(VonMises(observedState,1.0/σy), dataObs[t-1]) # Simulation prob given the data observed at that point
        end

        w=w_tilde./sum(w_tilde) #Normalize the weights

        #=
        Find the mean and std of the phase variable at each time step
        Z=sum(w .* exp.(im*xpf[t,:,2]))
        m1[t-1]=angle(Z) #mean phase estimate
        stdpf1[t-1]=abs(Z) #coherence of angle estimates
        =#


        #Add a check to only resample if the effective sample size is below N/10
        ESS=1.0/sum(w.^2)
        if ESS< N/10.0
            s= sample(1:N, Weights(w), N) #sample with replacement according to the weights
            xpf= xpf[:,s,:] #Note: resample WHOLE path, not just x.pf[t, ]
        end
        wloglik+=log(mean(w_tilde)) #caculate the log marginal likelihood of the data
        tstart+=24.0
    end
    return xpf, wloglik
end

function pMMH_likelihood(paramsProposal, dataObs, σp, σy; Nin=1000)
    xpf, wloglik=runParticleFilter(paramsProposal, dataObs, σp, σy, N=Nin)
    #s=sample(1:N, Weights(exp.(wloglik))) #pick a trajectory
    return wloglik
end

function pMMH(dataObs, M, σp, σy)
    #=
        Implement of Particle filter metropolis hastings walk through the
        parameter space.
            dataObs: observed phase every 24 hours with noise σy
            M: the number of steps to take in the walk
            σp: this is the length 2 vector of instrinsic noise in the system dynamics
            σy: This is is observational noise in the experimental measurements
        Returns the parameter values θ for each accepted step in the walk
    =#
    numParams=9
    priorCov=getCovarianceParameters(inflatedFactor=1.0)
    priorDist=MultivariateNormal(SinglePopModel.pvalues, priorCov)
    #Use the prior covariance matrix in the proposals to increase the rate of acceptances?
    # I think this is okay to do........
    proposalDist=MultivariateNormal(zeros(numParams), 0.5 .* priorCov)
    θlast=rand(priorDist)

    likθlast=pMMH_likelihood(θlast, dataObs, σp, σy)

    θ=zeros(M,numParams)

    for i in 1:M
        θstar=rand(proposalDist)+θlast
        likθstar=pMMH_likelihood(θlast, dataObs, σp, σy)
        logα=(likθstar-likθlast)+log(pdf(priorDist, θstar)/pdf(priorDist, θlast))
        if (logα>0.0)
            r=1.0
        else
            r=exp(logα)
        end
        if (rand()<=r)
            #Switch parameter values
            likθlast=likθstar
            θlast=θstar
        end
        θ[i,:]=θlast
        println(θlast)
    end
    return(θ)
end

function getCovarianceParameters(;filename="./mcmc_run_params.dat", inflatedFactor=1.0)
    #=
    Use the PRC fitting MCMC to form a prior for the parameters
    =#
    covM=diagm(0=>ones(9))
    println(size(covM))
    mcmcData=CSV.File("mcmc_run_params.dat", delim="\t") |> DataFrame
    mcmcData[:,1]= 2π ./ mcmcData[:,1] #convert period to frequency
    mcmcData=mcmcData[1:7] #take the first 7 columns only
    mcmcData[8]=0.024*ones(length(mcmcData[1])) #add fixed γ column
    mcmcData[9]=-0.0931817*ones(length(mcmcData[1])) #add a fixed β1 column


    permutecols!(mcmcData, [1,2,8,9,3,4,5,6,7]) #Arrange the order of the parameters so that it matches SinglePopModel
    mcmcData=Matrix(mcmcData) #convert to a matrix to get the covariance parameters

    covM=inflatedFactor*cov(mcmcData)
    #Set small entries to zero
    ϵcut=1.0e-8
    for i in 1:9
        for j in 1:9
            if abs(covM[i,j])<=ϵcut
                covM[i,j]=0.0
            end
        end
    end
    covM[3,3]=1e-6
    covM[4,4]=1e-6

    return covM
end

function runHCHS_PF(filename)
    SinglePopModel.setParameters()
    L, hchsDataFrame=HCHSLight.readData(filename)
    initModel(L)
    #@time ParticleFilterSP.systemDynamics([0.7,0.0], SinglePopModel.pvalues, [0.01,10.0], 0.0,24.0)


    trueStates, dataObs=makeFakeData(SinglePopModel.pvalues,[0.01,10.0], 10.0)
    xpf, wloglik=runParticleFilter(SinglePopModel.pvalues, dataObs, [0.01, 10.0], 10.0)
    p=plot(xpf[:,:,2], label="")
    scatter!(p,unwrap(trueStates), label="true states", color=:red)
    scatter!(p,)
    display(p)
    return trueStates,xpf

end

function PMCMC_HCHS(filename)
    SinglePopModel.setParameters()
    L, hchsDataFrame=HCHSLight.readData(filename)
    initModel(L)

    trueStates, dataObs=makeFakeData(SinglePopModel.pvalues,[0.01,10.0], 10.0)
    θ=pMMH(dataObs, 100, [0.01,10.0], 10.0)
    return(θ)
end



#trueStates,xpf=runHCHS_PF("hchs-sol-sueno-00579338.csv")
#PMCMC_HCHS("hchs-sol-sueno-00579338.csv")


end
