# callbacks.jl
# Define callbacks for periodic dilution
# Separate definitions for different models
# July 2020, Nadanai Laohakunakorn (nadanai.laohakunakorn@ed.ac.uk)


###############################################
######### 1. Resource dependent model #########
###############################################

function diluteSR!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_R,idx_dG,idx_pT,idx_dT]
        CONC_REFRESH = [R0,dG0,pT0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_R,idx_dG,idx_dT]
        CONC_REFRESH = [R0,dG0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_R,idx_dG]
        CONC_REFRESH = [R0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function dilutePC!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_R,idx_dG,idx_pT]
        CONC_REFRESH = [R0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_R,idx_dG,idx_pT]
        CONC_REFRESH = [R0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_R,idx_dG,idx_pT]
        CONC_REFRESH = [R0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function diluteNC!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_R,idx_dG,idx_pT]
        CONC_REFRESH = [R0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_R,idx_dG]
        CONC_REFRESH = [R0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_R,idx_dG]
        CONC_REFRESH = [R0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end


###############################################
######## 2. Resource independent model ########
###############################################

function diluteSR_RI!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT,idx_dT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_dT]
        CONC_REFRESH = [TX0,TL0,dG0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function dilutePC_RI!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function diluteNC_RI!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end


###############################################
## 3. Resource independent model  no loading ##
###############################################

function diluteSR_RInl!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT,idx_dT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_dT]
        CONC_REFRESH = [TX0,TL0,dG0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function dilutePC_RInl!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function diluteNC_RInl!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end


###############################################
######## 4. Coupled model TL saturation #######
###############################################

function diluteSR_1RI!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT,idx_dT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_dT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [R0,TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function dilutePC_1RI!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function diluteNC_1RI!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [R0,TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [R0,TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

###############################################
####### 5. Coupled model TXTL saturation ######
###############################################

function diluteSR_1RIf!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT,idx_dT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_dT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,dT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [R0,TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function dilutePC_1RIf!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end

function diluteNC_1RIf!(integrator)
    if integrator.t<SWITCHTIMES[1]*60
        # Stage 1
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG,idx_pT]
        CONC_REFRESH = [R0,TX0,TL0,dG0,pT0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[1]*60<=integrator.t<SWITCHTIMES[2]*60
        # Stage 2
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [R0,TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    elseif SWITCHTIMES[2]*60<=integrator.t
        # Stage 3
        INDEX_REFRESH = [idx_R,idx_TX,idx_TL,idx_dG]
        CONC_REFRESH = [R0,TX0,TL0,dG0]
        for j in 1:NSPECIES
            integrator.u[j] = integrator.u[j]*(1-DIL_FRAC)
        end
        for j in 1:size(INDEX_REFRESH)[1]
            integrator.u[INDEX_REFRESH[j]] = integrator.u[INDEX_REFRESH[j]] + DIL_FRAC*CONC_REFRESH[j]
        end
        
    end
end