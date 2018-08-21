*** Flux balance analysis with unlimited uptake and amino acids as carbon sources (scond)
*** required input: --scond = x (if not set scond = glucose+AA)
***                 --rep = x (if not set rep = 1) ... iteration

$If NOT set scond $set scond glucoseAA
$If NOT set rep $set rep 1   

$If not dexist results $call "mkdir results"
$If not dexist results/unlimCsource $call "mkdir results/unlimCsource"
$If not dexist results/unlimCsource/%scond%_%rep% $call "mkdir results/unlimCsource/%scond%_%rep%"

$onecho > ./results/unlimCsource/%scond%_%rep%/cplex.op3
interactive 1
iatriggerfile results/unlimCsource/%scond%_%rep%/stopme
iatriggertime 43200
iafile results/unlimCsource/%scond%_%rep%/cplex.op4
threads 4
parallelmode -1
eprhs 1e-9
epopt 1e-9
epagap 1e-9
epgap 1e-9
epint 1e-9
mipstart 1
memoryemphasis 1
solnpool results/unlimCsource/%scond%_%rep%/solnpool.gdx
solnpoolreplace 2
solnpoolcapacity 1000
populatelim 1000
solnpoolprefix soln_unlimCsource_%scond%_%rep%
solnpoolpop 2
$offecho

$onecho > results/unlimCsource/%scond%_%rep%/cplex.op4
epagap 1e9
$offecho

$echo stopme > results/unlimCsource/%scond%_%rep%/stopme

$onecho > ./results/unlimCsource/%scond%_%rep%/conopt.opt
LFNICR 100
LS2NDI 1
RTBND1 1e-9
RTBNDT 1e-10
RTMAXJ 1e15
RTMAXV 1e15
RTNWMA 1e-9
RTNWMI 5e-11
RTNWTR 1e-9
RTREDG 1e-9
RVHESS 0
RVTIME 43200
$offecho

$onecho > ./results/unlimCsource/%scond%_%rep%/conopt4.opt
LFNICR 100
RTNWMA 1e-9
RTNWMI 5e-11
RTNWTR 1e-9
RTREDG 1e-9
RVTIME 43200
Num_Rounds 4
$offecho

$onecho > ./results/unlimCsource/%scond%_%rep%/antigone.opt
abs_opt_tol=1e-6
cplex_optfile=optfiles/cplex.op2
conopt_optfile conopt.opt
feas_tolerance=1E-9
logging_freq 50
logging_level=2
max_rlt_cuts=100
max_time=50
nominal_time_limit=1000
feas_soln_time_limit=4000
print_options=1
rel_opt_tol=1e-7
piecewise_linear_partitions=1
number_of_partitions=2
eigenvector_projection_partitioning=1
max_partitioned_quantities=50
$offecho

$onecho > ./results/unlimCsource/%scond%_%rep%/cplex.op2
threads 4
epint 1e-9
epopt 1e-9
eprhs 1e-9
epagap 0.1
epgap 0.1
memoryemphasis 1
$offecho

option minlp=antigone, mip=cplex, nlp=conopt, optcr = 1e-9;

** load sets and parameters needed to construct the model
sets mets,comps,rxns,grp,exchanges,sndlaw,gav,bndID,nsID,consdrGerr(rxns),consNS(nsID),blocked(rxns);
scalar temp,gascons,epsilon,CI;
parameters S(mets,comps,rxns),Sg(mets,comps,grp),drGt0tr(rxns),drGSE(rxns),dfG0(mets,comps),ratio(grp,rxns),sigmacLimits(bndID),lncLimits(mets,*,bndID),dfGt0(mets,comps)
vgLimits(grp,bndID),drGLimits(rxns,bndID),dfGLimits(rxns,bndID),drGgLimits(grp,bndID),dfGgLimits(grp,bndID),K(rxns,nsID);

$GDXIN data/model_Ecoli.gdx
** sets
$LOAD mets
$LOAD comps
$LOAD rxns
$LOAD grp
$LOAD exchanges
$LOAD sndlaw
$LOAD gav
$LOAD nsID
$LOAD consdrGerr
$LOAD bndID
$LOAD consNS
$LOAD blocked
** scalars
$LOAD Temp
$LOAD GasCons
$LOAD epsilon
$LOAD CI
$LOAD S
** parameter
$LOAD Sg
$LOAD dfG0
$LOAD drGt0tr
$LOAD drGSE
$LOAD dfGt0
$LOAD ratio
$LOAD K
** limits
$LOAD sigmacLimits
$LOAD lncLimits
$LOAD vgLimits
$LOAD drGLimits
$LOAD dfGLimits
$LOAD drGgLimits
$LOAD dfGgLimits
$GDXIN

** load literature variable bounds (rMIPbnds were determined using literature concentration bounds and uptake constraints as specified below)
parameter V_vgLimits(grp,bndID),V_drGgLimits(grp,bndID),V_dfGgLimits(grp,bndID),V_sigmagLimits(grp,bndID),V_gexgLimits(grp,bndID),V_drGLimits(rxns,bndID),V_lncLimits(mets,comps,bndID);

$GDXIN data/bnds/literature_bnds_rMIP_%scond%500.gdx
$LOAD V_vgLimits
$LOAD V_drGgLimits
$LOAD V_dfGgLimits
$LOAD V_sigmagLimits
$LOAD V_gexgLimits
$LOAD V_drGLimits
$LOAD V_lncLimits
$GDXIN

lncLimits(mets,comps,'lo')$((lncLimits(mets,comps,'lo') ne lncLimits(mets,comps,'up')) and (V_lncLimits(mets,comps,'lo') ne 0) and (lncLimits(mets,comps,'lo') le V_lncLimits(mets,comps,'lo')) and (V_lncLimits(mets,comps,'lo') le V_lncLimits(mets,comps,'up'))) = V_lncLimits(mets,comps,'lo');
lncLimits(mets,comps,'up')$((lncLimits(mets,comps,'lo') ne lncLimits(mets,comps,'up')) and (V_lncLimits(mets,comps,'up') ne 0) and  (lncLimits(mets,comps,'up') ge V_lncLimits(mets,comps,'up')) and (V_lncLimits(mets,comps,'lo') le V_lncLimits(mets,comps,'up'))) = V_lncLimits(mets,comps,'up');
vgLimits(grp,'lo')$((vgLimits(grp,'lo') ne vgLimits(grp,'up')) and (V_vgLimits(grp,'lo') ne 0) and  (vgLimits(grp,'lo') le V_vgLimits(grp,'lo')) and (V_vgLimits(grp,'lo') le V_vgLimits(grp,'up'))) = V_vgLimits(grp,'lo');
vgLimits(grp,'up')$((vgLimits(grp,'lo') ne vgLimits(grp,'up')) and (V_vgLimits(grp,'up') ne 0) and  (vgLimits(grp,'up') ge V_vgLimits(grp,'up')) and (V_vgLimits(grp,'lo') le V_vgLimits(grp,'up'))) = V_vgLimits(grp,'up');
drGLimits(rxns,'lo')$((drGLimits(rxns,'lo') ne drGLimits(rxns,'up')) and (V_drGLimits(rxns,'lo') ne 0) and  (drGLimits(rxns,'lo') le V_drGLimits(rxns,'lo')) and (V_drGLimits(rxns,'lo') le V_drGLimits(rxns,'up'))) = V_drGLimits(rxns,'lo');
drGLimits(rxns,'up')$((drGLimits(rxns,'lo') ne drGLimits(rxns,'up')) and (V_drGLimits(rxns,'up') ne 0) and  (drGLimits(rxns,'up') ge V_drGLimits(rxns,'up')) and (V_drGLimits(rxns,'lo') le V_drGLimits(rxns,'up'))) = V_drGLimits(rxns,'up');
drGgLimits(grp,'lo')$((drGgLimits(grp,'lo') ne drGgLimits(grp,'up')) and (V_drGgLimits(grp,'lo') ne 0) and  (drGgLimits(grp,'lo') le V_drGgLimits(grp,'lo')) and (V_drGgLimits(grp,'lo') le V_drGgLimits(grp,'up'))) = V_drGgLimits(grp,'lo');
drGgLimits(grp,'up')$((drGgLimits(grp,'lo') ne drGgLimits(grp,'up')) and (V_drGgLimits(grp,'up') ne 0) and  (drGgLimits(grp,'up') ge V_drGgLimits(grp,'up')) and (V_drGgLimits(grp,'lo') le V_drGgLimits(grp,'up'))) = V_drGgLimits(grp,'up');
dfGgLimits(grp,'lo')$((dfGgLimits(grp,'lo') ne dfGgLimits(grp,'up')) and (V_dfGgLimits(grp,'lo') ne 0) and  (dfGgLimits(grp,'lo') le V_dfGgLimits(grp,'lo')) and (V_dfGgLimits(grp,'lo') le V_dfGgLimits(grp,'up'))) = V_dfGgLimits(grp,'lo');
dfGgLimits(grp,'up')$((dfGgLimits(grp,'lo') ne dfGgLimits(grp,'up')) and (V_dfGgLimits(grp,'up') ne 0) and  (dfGgLimits(grp,'up') ge V_dfGgLimits(grp,'up')) and (V_dfGgLimits(grp,'lo') le V_dfGgLimits(grp,'up'))) = V_dfGgLimits(grp,'up');

** set all uptake rates to 0 but allow production of acetate, succinate, ethanol, formate
vgLimits(grp,bndID)$(ratio(grp,'EX_fru') ne 0) = 0;
vgLimits(grp,bndID)$(ratio(grp,'EX_gal') ne 0) = 0;
vgLimits(grp,bndID)$(ratio(grp,'EX_glcn') ne 0) = 0;
vgLimits(grp,bndID)$(ratio(grp,'EX_glyc') ne 0) = 0;
vgLimits(grp,'up')$(ratio(grp,'EX_ac') ne 0) = 0;
vgLimits(grp,bndID)$(ratio(grp,'EX_pyr') ne 0) = 0;
vgLimits(grp,'lo')$(ratio(grp,'EX_succ') ne 0) = 0;
vgLimits(grp,'lo')$(ratio(grp,'EX_etoh') ne 0) = 0;
vgLimits(grp,'lo')$(ratio(grp,'EX_for') ne 0) = 0;
vgLimits(grp,bndID)$(ratio(grp,'EX_glc') ne 0) = 0;
vgLimits(grp,'lo')$(ratio(grp,'EX_co2') ne 0) = 0;

** allow the uptake of the specified additional carbon source
set cond/glucoseAA,glycerolAA/,runcond(cond);
runcond(cond)=no;
runcond('%scond%')=yes;
loop(cond$runcond(cond),
         if(sameas(cond,'glucoseAA'),vgLimits(grp,'lo')$(ratio(grp,'EX_glc') eq 1)=-500);
         if(sameas(cond,'glycerolAA'),vgLimits(grp,'lo')$(ratio(grp,'EX_glyc') eq 1)=-500);
);

** include model equations from tcbm.inc
$include data/tcbm_v2bnd.inc

** define growth as objective function
variable mu;
equation defmu;
defmu.. mu =e= sum(grp,ratio(grp,'biomass')*vg(grp));

** define model
model FBA_mip /defmu,tcbm_mip,nullspace,tcbm_McCormick/;
FBA_mip.optfile = 3;
FBA_mip.holdfixed = 1;
FBA_mip.solprint = 2;

model FBA_nlp /defmu,tcbm_nlp,nullspace,tcbm_McCormick/;
FBA_nlp.optfile = 1;
FBA_nlp.holdfixed = 1;
FBA_nlp.solprint = 2;

option reslim = 54000;

** solve MIP linear approximation to generate start value
solve FBA_mip us mip max mu;

set soln           possible solutions in the solution pool /file1*file1000/
    solnpool(soln) actual solutions;
file fsol;

execute_load 'results/unlimCsource/%scond%_%rep%/solnpool.gdx', solnpool=Index;

loop(solnpool(soln),
  put_utility fsol 'shell' / 'mv ./'solnpool.te(soln):0:0' ./results/unlimCsource/%scond%_%rep%';
);

loop(solnpool(soln),
  put_utility fsol 'gdxin' / 'results/unlimCsource/%scond%_%rep%/'solnpool.te(soln):0:0;
  execute_loadpoint;
  solve FBA_nlp us nlp max mu;
  if (((FBA_nlp.modelstat eq 2) or (FBA_nlp.modelstat eq 8) or (FBA_nlp.modelstat eq 1) or (FBA_nlp.modelstat eq 7)),
        put_utility fsol 'gdxout' / 'results/unlimCsource/%scond%_%rep%/NLP_'solnpool.te(soln):0:0;
        execute_unload drGerror, sigmac, vg, drGg, dfGg, sigmag, gexg, drG, dfG, lnc, v, g, b, mu,
                       active_drG, active_drGg, active_dfG, active_dfGg, active_snd, active_sndg, active_gexg, active_sigmag,
                       rxns,mets,comps,grp;
  );
);

loop(solnpool(soln),
  put_utility fsol 'shell' / 'rm ./results/unlimCsource/%scond%_%rep%/'solnpool.te(soln):0:0'';
);

put_utility fsol 'shell' / 'rm ./results/unlimCsource/%scond%_%rep%/antigone*';
put_utility fsol 'shell' / 'rm ./results/unlimCsource/%scond%_%rep%/cplex*';
put_utility fsol 'shell' / 'rm ./results/unlimCsource/%scond%_%rep%/conopt*';
