# solve.jl
# Solve ODEs in three-stage chemostat experiment
# Callbacks determine when to dilute and which species to refresh
# Uses CVODE_BDF solver from Sundials
# July 2020, Nadanai Laohakunakorn (nadanai.laohakunakorn@ed.ac.uk)

using Sundials

function marktime!(integrator)   
end

function condition(u,t,integrator) 
    t%INTERVAL_DIL 
end

###############################################
######### 1. Resource dependent model #########
###############################################

function solvemodel(grads,u0,params,TMAX,INTERVAL_DIL,tsave,DIL_FRAC,NSPECIES,SWITCHTIMES,FLAG)
	tspan = (0.0,TMAX); 
	prob = ODEProblem(grads,u0,tspan,params);

	# Callbacks
	saved_values = SavedValues(Float64, Array{Float64});
	savecb = SavingCallback((u,t,integrator)->integrator(t,Val{1}),  # Save derivatives as well
	    saved_values,saveat=tsave);
    if FLAG=="SR"
	   contcb = ContinuousCallback(condition,diluteSR!;save_positions=(false,false))
    elseif FLAG=="PC"
        contcb = ContinuousCallback(condition,dilutePC!;save_positions=(false,false))
    elseif FLAG=="NC"   
        contcb = ContinuousCallback(condition,diluteNC!;save_positions=(false,false))
    end
	periodcb = PeriodicCallback(marktime!,INTERVAL_DIL;save_positions=(false,false)) # hack to mark time points for dilution
	cb = CallbackSet(periodcb,contcb,savecb)

	# Solve
	sol = solve(prob,CVODE_BDF(),abstol=1e-10,reltol=1e-10, callback=cb,saveat=tsave);
	return(sol,saved_values)
end

###############################################
######## 2. Resource independent model ########
###############################################

function solvemodel_RI(grads,u0,params,TMAX,INTERVAL_DIL,tsave,DIL_FRAC,NSPECIES,SWITCHTIMES,FLAG)
    tspan = (0.0,TMAX); 
    prob = ODEProblem(grads,u0,tspan,params);

    # Callbacks
    saved_values = SavedValues(Float64, Array{Float64});
    savecb = SavingCallback((u,t,integrator)->integrator(t,Val{1}),  # Save derivatives as well
        saved_values,saveat=tsave);
    if FLAG=="SR"
       contcb = ContinuousCallback(condition,diluteSR_RI!;save_positions=(false,false))
    elseif FLAG=="PC"
        contcb = ContinuousCallback(condition,dilutePC_RI!;save_positions=(false,false))
    elseif FLAG=="NC"   
        contcb = ContinuousCallback(condition,diluteNC_RI!;save_positions=(false,false))
    end
    periodcb = PeriodicCallback(marktime!,INTERVAL_DIL;save_positions=(false,false)) # hack to mark time points for dilution
    cb = CallbackSet(periodcb,contcb,savecb)

    # Solve
    sol = solve(prob,CVODE_BDF(),abstol=1e-10,reltol=1e-10, callback=cb,saveat=tsave);
    return(sol,saved_values)
end


###############################################
## 3. Resource independent model  no loading ##
###############################################
function solvemodel_RInl(grads,u0,params,TMAX,INTERVAL_DIL,tsave,DIL_FRAC,NSPECIES,SWITCHTIMES,FLAG)
    tspan = (0.0,TMAX); 
    prob = ODEProblem(grads,u0,tspan,params);

    # Callbacks
    saved_values = SavedValues(Float64, Array{Float64});
    savecb = SavingCallback((u,t,integrator)->integrator(t,Val{1}),  # Save derivatives as well
        saved_values,saveat=tsave);
    if FLAG=="SR"
       contcb = ContinuousCallback(condition,diluteSR_RInl!;save_positions=(false,false))
    elseif FLAG=="PC"
        contcb = ContinuousCallback(condition,dilutePC_RInl!;save_positions=(false,false))
    elseif FLAG=="NC"   
        contcb = ContinuousCallback(condition,diluteNC_RInl!;save_positions=(false,false))
    end
    periodcb = PeriodicCallback(marktime!,INTERVAL_DIL;save_positions=(false,false)) # hack to mark time points for dilution
    cb = CallbackSet(periodcb,contcb,savecb)

    # Solve
    sol = solve(prob,CVODE_BDF(),abstol=1e-10,reltol=1e-10, callback=cb,saveat=tsave);
    return(sol,saved_values)
end

###############################################
######## 4. Coupled model TL saturation #######
###############################################

function solvemodel_1RI(grads,u0,params,TMAX,INTERVAL_DIL,tsave,DIL_FRAC,NSPECIES,SWITCHTIMES,FLAG)
    tspan = (0.0,TMAX); 
    prob = ODEProblem(grads,u0,tspan,params);

    # Callbacks
    saved_values = SavedValues(Float64, Array{Float64});
    savecb = SavingCallback((u,t,integrator)->integrator(t,Val{1}),  # Save derivatives as well
        saved_values,saveat=tsave);
    if FLAG=="SR"
       contcb = ContinuousCallback(condition,diluteSR_1RI!;save_positions=(false,false))
    elseif FLAG=="PC"
        contcb = ContinuousCallback(condition,dilutePC_1RI!;save_positions=(false,false))
    elseif FLAG=="NC"   
        contcb = ContinuousCallback(condition,diluteNC_1RI!;save_positions=(false,false))
    end
    periodcb = PeriodicCallback(marktime!,INTERVAL_DIL;save_positions=(false,false)) # hack to mark time points for dilution
    cb = CallbackSet(periodcb,contcb,savecb)

    # Solve
    sol = solve(prob,CVODE_BDF(),abstol=1e-10,reltol=1e-10, callback=cb,saveat=tsave);
    return(sol,saved_values)
end

###############################################
####### 5. Coupled model TXTL saturation ######
###############################################

function solvemodel_1RIf(grads,u0,params,TMAX,INTERVAL_DIL,tsave,DIL_FRAC,NSPECIES,SWITCHTIMES,FLAG)
    tspan = (0.0,TMAX); 
    prob = ODEProblem(grads,u0,tspan,params);

    # Callbacks
    saved_values = SavedValues(Float64, Array{Float64});
    savecb = SavingCallback((u,t,integrator)->integrator(t,Val{1}),  # Save derivatives as well
        saved_values,saveat=tsave);
    if FLAG=="SR"
       contcb = ContinuousCallback(condition,diluteSR_1RIf!;save_positions=(false,false))
    elseif FLAG=="PC"
        contcb = ContinuousCallback(condition,dilutePC_1RIf!;save_positions=(false,false))
    elseif FLAG=="NC"   
        contcb = ContinuousCallback(condition,diluteNC_1RIf!;save_positions=(false,false))
    end
    periodcb = PeriodicCallback(marktime!,INTERVAL_DIL;save_positions=(false,false)) # hack to mark time points for dilution
    cb = CallbackSet(periodcb,contcb,savecb)

    # Solve
    sol = solve(prob,CVODE_BDF(),abstol=1e-10,reltol=1e-10, callback=cb,saveat=tsave);
    return(sol,saved_values)
end