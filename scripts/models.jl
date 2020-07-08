# models.jl
# define ODEs
# 5 models: resource dependent, resource independent
# no TL coupling, with TL coupling, with TX and TL coupling
# July 2020, Nadanai Laohakunakorn (nadanai.laohakunakorn@ed.ac.uk)

###############################################
######### 1. Resource dependent model #########
###############################################

function grads(du,u,params,t)
	
    R = u[1]
    dT = u[2]
    dG = u[3]
    mT = u[4]
    mG = u[5]
    pT = u[6]
    pG = u[7]

    alpha = params[1]
    beta = params[2]
    K = params[3]

    du[1] = -alpha*R/(R+K)*dT*pT-alpha*R/(R+K)*dG*pT-beta*R/(R+K)*mT-beta*R/(R+K)*mG
    du[2] = 0
    du[3] = 0
    du[4] = alpha*R/(R+K)*dT*pT
    du[5] = alpha*R/(R+K)*dG*pT
    du[6] = beta*R/(R+K)*mT
    du[7] = beta*R/(R+K)*mG

end

###############################################
######## 2. Resource independent model ########
###############################################

function grads_RI(du,u,params,t)
    
    R = u[1]
    dT = u[2]
    dG = u[3]
    mT = u[4]
    mG = u[5]
    pT = u[6]
    pG = u[7]
    TX = u[8]
    TL = u[9]

    alpha = params[1]
    beta = params[2]
    lamb1 = params[3]
    lamb2 = params[4]
    KTL = params[5]

    du[1] = 0
    du[2] = 0
    du[3] = 0
    du[4] = alpha*TX*dT*pT
    du[5] = alpha*TX*dG*pT
    du[6] = beta*TL*mT/(mT+mG+KTL)
    du[7] = beta*TL*mG/(mT+mG+KTL)
    du[8] = -lamb1*TX
    du[9] = -lamb2*TL

end

###############################################
## 3. Resource independent model  no loading ##
###############################################

function grads_RInl(du,u,params,t)
    
    R = u[1]
    dT = u[2]
    dG = u[3]
    mT = u[4]
    mG = u[5]
    pT = u[6]
    pG = u[7]
    TX = u[8]
    TL = u[9]

    alpha = params[1]
    beta = params[2]
    lamb1 = params[3]
    lamb2 = params[4]
    KTL = params[5]

    du[1] = 0
    du[2] = 0
    du[3] = 0
    du[4] = alpha*TX*dT*pT
    du[5] = alpha*TX*dG*pT
    du[6] = beta*TL*mT/(mT+KTL)
    du[7] = beta*TL*mG/(mG+KTL)
    du[8] = -lamb1*TX
    du[9] = -lamb2*TL

end

###############################################
######## 4. Coupled model TL saturation #######
###############################################

function grads_1RI(du,u,params,t)
    
    R = u[1]
    dT = u[2]
    dG = u[3]
    mT = u[4]
    mG = u[5]
    pT = u[6]
    pG = u[7]
    TX = u[8]
    TL = u[9]

    alpha = params[1]
    beta = params[2]
    lamb1 = params[3]
    lamb2 = params[4]
    KTL = params[5]
    K = params[6]

    du[1] = -alpha*TX*dT*pT*R/(R+K)-alpha*TX*dG*pT*R/(R+K)-beta*TL*mT/(mT+mG+KTL)*R/(R+K)-beta*TL*mG/(mT+mG+KTL)*R/(R+K)
    du[2] = 0
    du[3] = 0
    du[4] = alpha*TX*dT*pT*R/(R+K)
    du[5] = alpha*TX*dG*pT*R/(R+K)
    du[6] = beta*TL*mT/(mT+mG+KTL)*R/(R+K)
    du[7] = beta*TL*mG/(mT+mG+KTL)*R/(R+K)
    du[8] = -lamb1*TX
    du[9] = -lamb2*TL

end

###############################################
####### 5. Coupled model TXTL saturation ######
###############################################

function grads_1RIf(du,u,params,t)
    
    R = u[1]
    dT = u[2]
    dG = u[3]
    mT = u[4]
    mG = u[5]
    pT = u[6]
    pG = u[7]
    TX = u[8]
    TL = u[9]

    alpha = params[1]
    beta = params[2]
    lamb1 = params[3]
    lamb2 = params[4]
    KTL = params[5]
    K = params[6]
    KTX = params[7]

    du[1] = -alpha*TX*dT/(dT+dG+KTX)*pT*R/(R+K)-alpha*TX*dG/(dT+dG+KTX)*pT*R/(R+K)-beta*TL*mT/(mT+mG+KTL)*R/(R+K)-beta*TL*mG/(mT+mG+KTL)*R/(R+K)
    du[2] = 0
    du[3] = 0
    du[4] = alpha*TX*dT/(dT+dG+KTX)*pT*R/(R+K)
    du[5] = alpha*TX*dG/(dT+dG+KTX)*pT*R/(R+K)
    du[6] = beta*TL*mT/(mT+mG+KTL)*R/(R+K)
    du[7] = beta*TL*mG/(mT+mG+KTL)*R/(R+K)
    du[8] = -lamb1*TX
    du[9] = -lamb2*TL

end