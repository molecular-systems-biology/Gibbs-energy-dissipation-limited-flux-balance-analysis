*** Flux balance analysis with limited glucose uptake rate (GUR)
*** required input: --GUR = x (if not set GUR = 1)

$If NOT set GUR $set GUR 1

$If not dexist results $call "mkdir results"

option minlp=antigone, mip=cplex, nlp=conopt;

$onecho > cplex.opt
threads 1
parallelmode -1
eprhs 1e-6
epopt 1e-6
epagap 1e-6
epgap 1e-7
epint 1e-6
$offecho

$onecho > antigone.opt
abs_opt_tol=1e-6
cplex_optfile=cplex.opt
feas_tolerance=1E-6
logging_freq=10
logging_level=-1
max_rlt_cuts=100
max_time=36000
nominal_time_limit=60
print_options=1
rel_opt_tol=1e-7
$offecho

** load sets and parameters needed to construct the model
sets mets,comps,rxns,grp,exchanges,sndlaw,gav,bndID,nsID;
scalar temp,gascons,epsilon;
parameters S(mets,comps,rxns),Sg(mets,comps,grp),dfG0(mets,comps),drG0(rxns),ratio(grp,rxns),sigmacLimits(bndID),lncLimits(mets,comps,bndID)
vgLimits(grp,bndID),drGLimits(rxns,bndID),dfGLimits(rxns,bndID),drGgLimits(grp,bndID),dfGgLimits(grp,bndID);

$GDXIN data/model_Scerevisiae.gdx
** sets
$LOAD mets
$LOAD comps
$LOAD rxns
$LOAD grp
$LOAD exchanges
$LOAD sndlaw
$LOAD gav
$LOAD bndID
** scalars
$LOAD Temp
$LOAD GasCons
$LOAD epsilon
$LOAD S
** parameter
$LOAD Sg
$LOAD dfG0
$LOAD drG0
$LOAD ratio
** limits
$LOAD sigmacLimits
$LOAD lncLimits
$LOAD vgLimits
$LOAD drGLimits
$LOAD dfGLimits
$LOAD drGgLimits
$LOAD dfGgLimits
$GDXIN

** include model equations from tcbm.inc
$include data/tcbm.inc

** define growth as objective function
variable mu;
equation defmu;
defmu.. mu =e= sum(grp,ratio(grp,'biomass')*vg(grp));

** define model
model fba_mip /defmu,tcbm_mip,tcbm_McCormick/;
fba_mip.optfile = 1;
fba_mip.holdfixed = 1;

model fba_minlp /defmu,tcbm_minlp,tcbm_McCormick/;
fba_minlp.optfile = 1;
fba_minlp.holdfixed = 1;

** solve MIP linear approximation to generate start value
solve fba_mip us mip max mu;

** solve model with start value from MIP solution
solve fba_minlp us minlp max mu;

execute_unload'results/FBA_%GUR%.gdx';
