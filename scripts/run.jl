# run.jl
# Example script to run three-stage chemostat experiment
# July 2020, Nadanai Laohakunakorn (nadanai.laohakunakorn@ed.ac.uk)

using DelimitedFiles
using Plots
using DifferentialEquations
using ProgressBars
using Statistics
using PyCall

include("./models.jl");
include("./callbacks.jl");
include("./solve.jl");

# Global settings: paths and plot
Plots.pyplot()
PATH_OUT = "../output/"
FN="script_titration_plot.pdf" # Filename
fntsm = Plots.font("sans-serif", pointsize=round(14.0))
fntlg = Plots.font("sans-serif", pointsize=round(18.0))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntlg, legendfont=fntsm)

# Global simulation settings
TMAX = 30.0*60 # in minutes
INTERVAL_DIL = 15.0 # in minutes
DIL_FRAC = 0.2;
NSPECIES = 7
SWITCHTIMES = [4.0,16.0] # Time to swap between experimental stages

# 1. Define names for species indices
idx_R=1;
idx_dT=2;
idx_dG=3;
idx_mT=4;
idx_mG=5;
idx_pT=6;
idx_pG=7;
idx_TX = 8;
idx_TL = 9;

# 2. Set up simulation and plot
CONCS = [0.001,0.01,1]
LABELS = ["SR low", "SR medium", "SR high"]
COLOURS = ["#afb1b6","#aec7e8","#1f77b4"]
p1 = plot(grid=:false,legend=:false)#,legend=:outertopright)
vspan!([4,16],fill=:black,alpha=:0.1); vline!([4,16],color=:black,linestyle=:dot,linewidth=2);

# 3. Set initial conditions and parameters
R0=100.0;
dT0=0.0;
dG0=2.0;
mT0=0.0;
mG0=0.0;
pT0=1.0;
pG0=0.0;
alpha=0.7;
beta=0.07;
K=1.0;
u0 = [R0,dT0,dG0,mT0,mG0,pT0,pG0];
params = [alpha,beta,K];

# 4. Solve 
#TSAVE = 0.0:5.0:TMAX # either array or range can be passed
TSAVE = collect(0:15:TMAX).+1 

# Positive control
solU,solDU=solvemodel(grads,u0,params,TMAX,INTERVAL_DIL,TSAVE,DIL_FRAC,NSPECIES,SWITCHTIMES,"PC");
t = solU.t/60;
pGPC = [datum for subarr in solU.u for datum in subarr[idx_pG]];

# Negative control
solU,solDU=solvemodel(grads,u0,params,TMAX,INTERVAL_DIL,TSAVE,DIL_FRAC,NSPECIES,SWITCHTIMES,"NC");
pGNC = [datum for subarr in solU.u for datum in subarr[idx_pG]];

plot!(t,pGPC,label="PC",color="#2ca02c",xaxis="time (h)",yaxis="GFP",lw=4);
plot!(t,pGNC,label="NC",color="#ffbb78",lw=4)

# T7RNAP DNA titration
for j in (1:size(CONCS)[1])
    dT0=CONCS[j]
    u0 = [R0,dT0,dG0,mT0,mG0,pT0,pG0];
    params = [alpha,beta,K];
    solU,solDU=solvemodel(grads,u0,params,TMAX,INTERVAL_DIL,TSAVE,DIL_FRAC,NSPECIES,SWITCHTIMES,"SR");
    pG = [datum for subarr in solU.u for datum in subarr[idx_pG]];
    plot!(t,pG,label=LABELS[j],color=COLOURS[j],lw=4)
end

p1 # Show plot
savefig(PATH_OUT*FN) # Save plot


