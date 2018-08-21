*** Themodynamic constrained based model
* note that in this model the entropy production rate (sigma) is models instead of the energy dissipation rate (g)
* the energy dissipation rate (g) is entropy production rate (sigma) * T

** set variable bounds
variables obj,sigmac,vg(grp),drGg(grp),dfGg(grp),sigmag(grp),gexg(grp),drG(rxns),dfG(rxns),lnc(mets,*),v(rxns),g(rxns);
positive variables vp(rxns),vn(rxns),drGp(rxns),drGn(rxns),sqrdrG(rxns),conc(mets,*),conccomb(mets);
binary variables b(rxns);

sigmac.lo = sigmacLimits('lo');
sigmac.up = 15.9;
vg.lo(grp) = vgLimits(grp,'lo');
vg.up(grp) = vgLimits(grp,'up');
v.lo(rxns) = sum(grp,ratio(grp,rxns)*vg.lo(grp)$(ratio(grp,rxns) gt 0) + ratio(grp,rxns)*vg.up(grp)$(ratio(grp,rxns) lt 0))-1e-4;
v.up(rxns) = sum(grp,ratio(grp,rxns)*vg.lo(grp)$(ratio(grp,rxns) lt 0) + ratio(grp,rxns)*vg.up(grp)$(ratio(grp,rxns) gt 0))+1e-4;
drG.lo(rxns) = drGLimits(rxns,'lo');
drG.up(rxns) = drGLimits(rxns,'up');
drGg.lo(grp) = drGgLimits(grp,'lo');
drGg.up(grp) = drGgLimits(grp,'up');
dfG.lo(rxns) = dfGLimits(rxns,'lo');
dfG.up(rxns) = dfGLimits(rxns,'up');
dfGg.lo(grp) = dfGgLimits(grp,'lo');
dfGg.up(grp) = dfGgLimits(grp,'up');
lnc.lo(mets,comps) = lncLimits(mets,comps,'lo');
lnc.up(mets,comps) = lncLimits(mets,comps,'up');
sigmag.lo(grp)$(sum(rxns$sndlaw(rxns),1$(ratio(grp,rxns) ne 0)) eq sum(rxns$(gav(rxns)),1$(ratio(grp,rxns) ne 0))) = 0;
sigmag.lo(grp)$(sum(rxns$sndlaw(rxns),1$(ratio(grp,rxns) ne 0)) ne sum(rxns$(gav(rxns)),1$(ratio(grp,rxns) ne 0))) = -sigmacLimits('up');
sigmag.up(grp) = sigmacLimits('up');
sigmag.lo(grp)$((sigmag.lo(grp) ne sigmag.up(grp)) and (sigmag.lo(grp) le V_sigmagLimits(grp,'lo')) and (V_sigmagLimits(grp,'lo') ne 0) and (V_sigmagLimits(grp,'lo') le V_sigmagLimits(grp,'up'))) = V_sigmagLimits(grp,'lo');
sigmag.up(grp)$((sigmag.lo(grp) ne sigmag.up(grp)) and (sigmag.up(grp) ge V_sigmagLimits(grp,'up')) and (V_sigmagLimits(grp,'up') ne 0) and (V_sigmagLimits(grp,'lo') le V_sigmagLimits(grp,'up'))) = V_sigmagLimits(grp,'up');
g.lo(rxns)$sndlaw(rxns) = -Temp*sigmacLimits('up');
g.up(rxns)$sndlaw(rxns) = 0;
g.up(rxns)$(gav(rxns) and not sndlaw(rxns)) = Temp*sigmacLimits('up');
gexg.lo(grp)$(sum(rxns$exchanges(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = -1e5;
gexg.up(grp)$(sum(rxns$exchanges(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = +1e5;
gexg.lo(grp)$((gexg.lo(grp) ne gexg.up(grp)) and (gexg.lo(grp) le V_gexgLimits(grp,'lo')) and (V_gexgLimits(grp,'lo') ne 0) and (V_gexgLimits(grp,'lo') le V_gexgLimits(grp,'up'))) = V_gexgLimits(grp,'lo');
gexg.up(grp)$((gexg.lo(grp) ne gexg.up(grp)) and (gexg.lo(grp) ge V_gexgLimits(grp,'up')) and (V_gexgLimits(grp,'up') ne 0) and (V_gexgLimits(grp,'lo') le V_gexgLimits(grp,'up'))) = V_gexgLimits(grp,'up');
sqrdrG.lo(rxns) = sqr(epsilon);

** set drG0 as determined in the regression
* drG0_CCM + drGerror * drGSE_CCM = drG0_reg
* drG0_CCM is the drG0 determined from the CCM
* drGSE is the std deviation determined from the CCM

variable drGerror_reg(rxns);
parameter drGerr_max(rxns),drGerr_min(rxns);
$GDXIN data/bnds/drGerror_reg.gdx
$LOAD drGerr_max
$LOAD drGerr_min
$GDXIN

variable drGerror(rxns);
drGerror.lo(rxns)$consdrGerr(rxns) = drGerr_min(rxns) - 0.15 * abs(drGerr_min(rxns));
drGerror.up(rxns)$consdrGerr(rxns) = drGerr_max(rxns) + 0.15 * abs(drGerr_min(rxns));

** identify active variables (e.g. not blocked fluxes) & adapt bounds accordingly
set  map(mets,comps),active_drG(rxns), active_drGg(grp),active_dfG(rxns),active_dfGg(grp),active_snd(rxns),active_sndg(grp),active_gexg(grp),active_sigmag(grp);
map(mets,comps)$(not sameas(mets,'h') and not sameas(mets,'charge') and sum(rxns,1$(S(mets,comps,rxns) ne 0) gt 0)) = yes;
active_gexg(grp)$(sum(rxns$exchanges(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = yes;
active_sigmag(grp)$(sum(rxns$gav(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = yes;
active_snd(rxns)$(sndlaw(rxns) and sign(drG.lo(rxns)) ne sign(drG.up(rxns)) and not blocked(rxns)) = yes;
active_sndg(grp)$(not sameas(grp,'g1') and sum(rxns$sndlaw(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = yes;
vg.lo(grp)$(active_sndg(grp) and drGg.up(grp) lt 0 and vg.lo(grp) lt 0) = 0;
vg.up(grp)$(active_sndg(grp) and drGg.lo(grp) gt 0 and vg.up(grp) gt 0) = 0;
v.lo(rxns) = max(sum(grp,ratio(grp,rxns)*vg.lo(grp)$(ratio(grp,rxns) gt 0) + ratio(grp,rxns)*vg.up(grp)$(ratio(grp,rxns) lt 0)),-500)-1e-4;
v.up(rxns) = min(sum(grp,ratio(grp,rxns)*vg.lo(grp)$(ratio(grp,rxns) lt 0) + ratio(grp,rxns)*vg.up(grp)$(ratio(grp,rxns) gt 0)),500)+1e-4;
active_sndg(grp)$(not sameas(grp,'g1') and sum(rxns$sndlaw(rxns),1$(ratio(grp,rxns) ne 0)) eq sum(rxns$(sndlaw(rxns) and not active_snd(rxns)),1$(ratio(grp,rxns) ne 0))) = no;
active_drG(rxns)$(gav(rxns) and (drG.lo(rxns) ne drG.up(rxns) and not blocked(rxns)) and sum(map,1$(S(map,rxns) ne 0) - 1$(S(map,rxns) ne 0 and lnc.lo(map) eq lnc.up(map))) gt 0) = yes;
active_dfG(rxns)$(exchanges(rxns) and (dfG.lo(rxns) ne dfG.up(rxns))) = yes;
active_drGg(grp)$(active_sigmag(grp) and (drGg.lo(grp) ne drGg.up(grp)) and sum(rxns$(ratio(grp,rxns) ne 0),1$active_drG(rxns)) gt 0) = yes;
active_dfGg(grp)$(active_gexg(grp) and (dfGg.lo(grp) ne dfGg.up(grp)) and sum(rxns$(ratio(grp,rxns) ne 0),1$active_dfG(rxns)) gt 0) = yes;

********** EQUATIONS **********

********** LP **********
equations massbalances,vgtov,gibbsenergybalance,defsigmac,defdrG,defdrGg,defdfG,defdfGg,ATPM;

* massbalance
massbalances(mets,comps).. sum(grp,Sg(mets,comps,grp)*vg(grp)) =e= 0;

* v, vg transformation, GEB & sigmac
vgtov(rxns).. v(rxns) =e= sum(grp,ratio(grp,rxns)*vg(grp));
gibbsenergybalance.. -Temp*sigmac =e= sum(grp$active_gexg(grp),gexg(grp));
defsigmac.. sigmac =e= sum(grp$active_sigmag(grp),sigmag(grp));

* drG, dfG + drGg, dfGg transformation
defdrG(rxns)$active_drG(rxns)..
drG(rxns) =e= drGt0tr(rxns) + drGerror(rxns)*drGSE(rxns)$consdrGerr(rxns) + GasCons*Temp*sum((mets,comps)$(map(mets,comps) and not sameas(mets,'h') and not sameas(mets,'charge')),S(mets,comps,rxns)*lnc(mets,comps));
defdrGg(grp)$active_drGg(grp).. drGg(grp) =e= sum(rxns$gav(rxns),ratio(grp,rxns)*drG(rxns));
defdfG(rxns)$active_dfG(rxns).. dfG(rxns) =e= -sum(mets$map(mets,'e'),S(mets,'e',rxns)*(dfG0(mets,'e')+GasCons*Temp*lnc(mets,'e')));
defdfGg(grp)$active_dfGg(grp).. dfGg(grp) =e= sum(rxns$exchanges(rxns),ratio(grp,rxns)*dfG(rxns));

* ATP maintenance
ATPM.. vg('g62') =g= 3.15 + 31.2622 * vg('g1');
model tcbm_lp /all/;

********** MIP **********
* MIP 2ed law
equation defvgn,defvp,defvn,defdrGpn,defdrGp,defdrGn,defstrict;
defvgn(rxns)$active_snd(rxns).. v(rxns) =e= vp(rxns) - vn(rxns);
defvp(rxns)$active_snd(rxns).. vp(rxns) =l= max(abs(v.lo(rxns)),abs(v.up(rxns)))*b(rxns);
defvn(rxns)$active_snd(rxns).. vn(rxns) =l= max(abs(v.lo(rxns)),abs(v.up(rxns)))*(1-b(rxns));
defdrGpn(rxns)$active_snd(rxns).. drG(rxns) =e= drGp(rxns)-drGn(rxns);
defdrGp(rxns)$active_snd(rxns).. drGp(rxns) =l= max(abs(drG.lo(rxns)),abs(drG.up(rxns)))*(1-b(rxns));
defdrGn(rxns)$active_snd(rxns).. drGn(rxns) =l= max(abs(drG.lo(rxns)),abs(drG.up(rxns)))*b(rxns);
defstrict(rxns)$active_snd(rxns).. drGp(rxns) + drGn(rxns) =g= epsilon;
model tcbm_mip /all/;

********** MINLP **********
* GEB
equation defsigma,defgex;
defsigma(grp)$active_sigmag(grp).. -Temp*sigmag(grp) =e= drGg(grp)*vg(grp);
defgex(grp)$active_gexg(grp).. gexg(grp) =e= dfGg(grp)*vg(grp);
model tcbm_minlp /all/;

********** NLP **********
equation defgb,strictnlp,strictnlp_lp,defg_mccormick1b,defg_mccormick2b,defg_mccormick3b,defg_mccormick4b;
defgb(rxns)$active_snd(rxns).. g(rxns) =e= v(rxns)*drG(rxns);
strictnlp(rxns)$active_snd(rxns).. sqr(drG(rxns)) =e= sqrdrG(rxns);
strictnlp_lp(rxns)$active_snd(rxns).. (drG.lo(rxns)+drG.up(rxns))*drG(rxns)-drG.lo(rxns)*drG.up(rxns) =g= sqrdrG(rxns);
defg_mccormick1b(rxns)$active_snd(rxns).. g(rxns) =g= v(rxns)*drG.up(rxns)+v.up(rxns)*drG(rxns)-v.up(rxns)*drG.up(rxns);
defg_mccormick2b(rxns)$active_snd(rxns).. g(rxns) =g= v(rxns)*drG.lo(rxns)+v.lo(rxns)*drG(rxns)-v.lo(rxns)*drG.lo(rxns);
defg_mccormick3b(rxns)$active_snd(rxns).. g(rxns) =l= v(rxns)*drG.up(rxns)+v.lo(rxns)*drG(rxns)-v.lo(rxns)*drG.up(rxns);
defg_mccormick4b(rxns)$active_snd(rxns).. g(rxns) =l= v(rxns)*drG.lo(rxns)+v.up(rxns)*drG(rxns)-v.up(rxns)*drG.lo(rxns);
model tcbm_nlp /tcbm_lp,strictnlp,strictnlp_lp,defsigma,defgex,defgb,defg_mccormick1b,defg_mccormick2b,defg_mccormick3b,defg_mccormick4b/;

********** nullspace **********
equation nullspace;
nullspace(nsID)$consNS(nsID)..
sum(rxns,K(rxns,nsID)*(sum((mets,comps)$(map(mets,comps) and not sameas(mets,'h') and not sameas(mets,'charge')),S(mets,comps,rxns)*dfGt0(mets,'c'))$(drGSE(rxns) lt 1000)+drGerror(rxns)*drGSE(rxns)$consdrGerr(rxns))) =e= 0;

********** McCormicks envelopes (used for linear approximation) **********
equations defsigma_mccormick1,defsigma_mccormick2,defsigma_mccormick3,defsigma_mccormick4,defgex_mccormick1,defgex_mccormick2,defgex_mccormick3,defgex_mccormick4;

defsigma_mccormick1(grp)$active_sigmag(grp).. -Temp*sigmag(grp) =g= vg(grp)*drGg.up(grp)+vg.up(grp)*drGg(grp)-vg.up(grp)*drGg.up(grp);
defsigma_mccormick2(grp)$active_sigmag(grp).. -Temp*sigmag(grp) =g= vg(grp)*drGg.lo(grp)+vg.lo(grp)*drGg(grp)-vg.lo(grp)*drGg.lo(grp);
defsigma_mccormick3(grp)$active_sigmag(grp).. -Temp*sigmag(grp) =l= vg(grp)*drGg.up(grp)+vg.lo(grp)*drGg(grp)-vg.lo(grp)*drGg.up(grp);
defsigma_mccormick4(grp)$active_sigmag(grp).. -Temp*sigmag(grp) =l= vg(grp)*drGg.lo(grp)+vg.up(grp)*drGg(grp)-vg.up(grp)*drGg.lo(grp);
defgex_mccormick1(grp)$active_gexg(grp).. gexg(grp) =g= vg(grp)*dfGg.up(grp)+vg.up(grp)*dfGg(grp)-vg.up(grp)*dfGg.up(grp);
defgex_mccormick2(grp)$active_gexg(grp).. gexg(grp) =g= vg(grp)*dfGg.lo(grp)+vg.lo(grp)*dfGg(grp)-vg.lo(grp)*dfGg.lo(grp);
defgex_mccormick3(grp)$active_gexg(grp).. gexg(grp) =l= vg(grp)*dfGg.up(grp)+vg.lo(grp)*dfGg(grp)-vg.lo(grp)*dfGg.up(grp);
defgex_mccormick4(grp)$active_gexg(grp).. gexg(grp) =l= vg(grp)*dfGg.lo(grp)+vg.up(grp)*dfGg(grp)-vg.up(grp)*dfGg.lo(grp);
model tcbm_McCormick /defsigma_mccormick1,defsigma_mccormick2,defsigma_mccormick3,defsigma_mccormick4,defgex_mccormick1,defgex_mccormick2,defgex_mccormick3,defgex_mccormick4/;

