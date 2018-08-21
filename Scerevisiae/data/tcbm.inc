*** Themodynamic constrained based model
* note that in this model the entropy production rate (sigma) is models instead of the energy dissipation rate (g)
* the energy dissipation rate (g) is entropy production rate (sigma) * T

** set variable bounds
variables obj,sigmac,vg(grp),drGg(grp),dfGg(grp),sigmag(grp),gexg(grp),drG(rxns),dfG(rxns),lnc(mets,comps),v(rxns),g(rxns);
positive variables vp(rxns),vn(rxns),drGp(rxns),drGn(rxns),sqrdrG(rxns),conc(mets,*),conccomb(mets);
binary variables b(rxns);

sigmac.lo = sigmacLimits('lo');
sigmac.up = sigmacLimits('up');

* vg bounds are scaled to GUR = 1
vg.lo(grp) = max(vgLimits(grp,'lo')*%gur%,-500);
vg.up(grp) = min(vgLimits(grp,'up')*%gur%,500);

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
g.lo(rxns)$sndlaw(rxns) = -Temp*sigmacLimits('up');
g.up(rxns)$sndlaw(rxns) = 0;
g.up(rxns)$(gav(rxns) and not sndlaw(rxns)) = Temp*sigmacLimits('up');
gexg.lo(grp)$(sum(rxns$exchanges(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = -1e5;
gexg.up(grp)$(sum(rxns$exchanges(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = +1e5;
sqrdrG.lo(rxns) = sqr(epsilon);

** identify active variables (e.g. not blocked fluxes) & adapt bounds accordingly
set  map(mets,comps),active_drG(rxns), active_drGg(grp),active_dfG(rxns),active_dfGg(grp),active_snd(rxns),active_sndg(grp),active_gexg(grp),active_sigmag(grp);
map(mets,comps)$(not sameas(mets,'h') and not sameas(mets,'charge') and sum(rxns,1$(S(mets,comps,rxns) ne 0) gt 0)) = yes;
active_gexg(grp)$(sum(rxns$exchanges(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = yes;
active_sigmag(grp)$(sum(rxns$gav(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = yes;
active_snd(rxns)$(sndlaw(rxns) and sign(drG.lo(rxns)) ne sign(drG.up(rxns))) = yes;
active_sndg(grp)$(sum(rxns$sndlaw(rxns),1$(ratio(grp,rxns) ne 0)) gt 0) = yes;
vg.lo(grp)$(active_sndg(grp) and drGg.up(grp) lt 0 and vg.lo(grp) lt 0) = 0;
vg.up(grp)$(active_sndg(grp) and drGg.lo(grp) gt 0 and vg.up(grp) gt 0) = 0;
v.lo(rxns) = max(sum(grp,ratio(grp,rxns)*vg.lo(grp)$(ratio(grp,rxns) gt 0) + ratio(grp,rxns)*vg.up(grp)$(ratio(grp,rxns) lt 0)),-500)-1e-4;
v.up(rxns) = min(sum(grp,ratio(grp,rxns)*vg.lo(grp)$(ratio(grp,rxns) lt 0) + ratio(grp,rxns)*vg.up(grp)$(ratio(grp,rxns) gt 0)),500)+1e-4;
active_sndg(grp)$(sum(rxns$sndlaw(rxns),1$(ratio(grp,rxns) ne 0)) eq sum(rxns$(sndlaw(rxns) and not active_snd(rxns)),1$(ratio(grp,rxns) ne 0))) = no;
active_drG(rxns)$(gav(rxns) and (drG.lo(rxns) ne drG.up(rxns)) and sum(map,1$(S(map,rxns) ne 0) - 1$(S(map,rxns) ne 0 and lnc.lo(map) eq lnc.up(map))) gt 0) = yes;
active_dfG(rxns)$(exchanges(rxns) and (dfG.lo(rxns) ne dfG.up(rxns))) = yes;
active_drGg(grp)$(active_sigmag(grp) and (drGg.lo(grp) ne drGg.up(grp)) and sum(rxns$(ratio(grp,rxns) ne 0),1$active_drG(rxns)) gt 0) = yes;
active_dfGg(grp)$(active_gexg(grp) and (dfGg.lo(grp) ne dfGg.up(grp)) and sum(rxns$(ratio(grp,rxns) ne 0),1$active_dfG(rxns)) gt 0) = yes;

********** EQUATIONS **********

********** LP **********
equations massbalances,vgtov,gibbsenergybalance,defsigmac,defdrG,defdrGg,defdfG,defdfGg;

* massbalance
massbalances(mets,comps).. sum(grp,Sg(mets,comps,grp)*vg(grp)) =e= 0;

* v, vg transformation, GEB & sigmac
vgtov(rxns).. v(rxns) =e= sum(grp,ratio(grp,rxns)*vg(grp));
gibbsenergybalance.. -Temp*sigmac =e= sum(grp$active_gexg(grp),gexg(grp));
defsigmac.. sigmac =e= sum(grp$active_sigmag(grp),sigmag(grp));

* drG, dfG and drGg, dfGg transformation
defdrG(rxns)$active_drG(rxns).. drG(rxns) =e= drG0(rxns)+GasCons*Temp*sum((mets,comps)$map(mets,comps),S(mets,comps,rxns)*lnc(mets,comps));
defdrGg(grp)$active_drGg(grp).. drGg(grp) =e= sum(rxns$gav(rxns),ratio(grp,rxns)*drG(rxns));
defdfG(rxns)$active_dfG(rxns).. dfG(rxns) =e= -sum(mets$map(mets,'e'),S(mets,'e',rxns)*(dfG0(mets,'e')+GasCons*Temp*lnc(mets,'e')));
defdfGg(grp)$active_dfGg(grp).. dfGg(grp) =e= sum(rxns$exchanges(rxns),ratio(grp,rxns)*dfG(rxns));
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

********** McCormicks envelopes (used for linear approximation) **********
equations defsigma_mccormick1,defsigma_mccormick2,defsigma_mccormick3,defsigma_mccormick4,defgex_mccormick1,defgex_mccormick2,defgex_mccormick3,defgex_mccormick4;

* McCormicks envelopes
defsigma_mccormick1(grp)$active_sigmag(grp).. -Temp*sigmag(grp) =g= vg(grp)*drGg.up(grp)+vg.up(grp)*drGg(grp)-vg.up(grp)*drGg.up(grp);
defsigma_mccormick2(grp)$active_sigmag(grp).. -Temp*sigmag(grp) =g= vg(grp)*drGg.lo(grp)+vg.lo(grp)*drGg(grp)-vg.lo(grp)*drGg.lo(grp);
defsigma_mccormick3(grp)$active_sigmag(grp).. -Temp*sigmag(grp) =l= vg(grp)*drGg.up(grp)+vg.lo(grp)*drGg(grp)-vg.lo(grp)*drGg.up(grp);
defsigma_mccormick4(grp)$active_sigmag(grp).. -Temp*sigmag(grp) =l= vg(grp)*drGg.lo(grp)+vg.up(grp)*drGg(grp)-vg.up(grp)*drGg.lo(grp);
defgex_mccormick1(grp)$active_gexg(grp).. gexg(grp) =g= vg(grp)*dfGg.up(grp)+vg.up(grp)*dfGg(grp)-vg.up(grp)*dfGg.up(grp);
defgex_mccormick2(grp)$active_gexg(grp).. gexg(grp) =g= vg(grp)*dfGg.lo(grp)+vg.lo(grp)*dfGg(grp)-vg.lo(grp)*dfGg.lo(grp);
defgex_mccormick3(grp)$active_gexg(grp).. gexg(grp) =l= vg(grp)*dfGg.up(grp)+vg.lo(grp)*dfGg(grp)-vg.lo(grp)*dfGg.up(grp);
defgex_mccormick4(grp)$active_gexg(grp).. gexg(grp) =l= vg(grp)*dfGg.lo(grp)+vg.up(grp)*dfGg(grp)-vg.up(grp)*dfGg.lo(grp);
model tcbm_McCormick /defsigma_mccormick1,defsigma_mccormick2,defsigma_mccormick3,defsigma_mccormick4,defgex_mccormick1,defgex_mccormick2,defgex_mccormick3,defgex_mccormick4/;