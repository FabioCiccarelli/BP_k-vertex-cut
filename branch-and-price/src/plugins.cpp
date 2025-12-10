/**@file   plugins.c
 * @brief  include SCIP plugins and custom plugins
 * @author Christopher Hojny
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "plugins.h"
#include "prop_problemsymmetry.h"

/* include header files here, such that the user only has to include
 * scipdefplugins.h
 */
#include "scip/branch_allfullstrong.h"
#include "scip/branch_cloud.h"
#include "scip/branch_distribution.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_gomory.h"
#include "scip/branch_inference.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_lookahead.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_multaggr.h"
#include "scip/branch_nodereopt.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"
#include "scip/branch_vanillafullstrong.h"
#include "scip/compr_largestrepr.h"
#include "scip/compr_weakcompr.h"
#include "scip/cons_and.h"
#include "scip/cons_benders.h"
#include "scip/cons_benderslp.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_cardinality.h"
#include "scip/cons_conjunction.h"
#include "scip/cons_countsols.h"
#include "scip/cons_cumulative.h"
#include "scip/cons_disjunction.h"
#include "scip/cons_indicator.h"
#include "scip/cons_integral.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_linking.h"
#include "scip/cons_logicor.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_or.h"
#include "scip/cons_orbisack.h"
#include "scip/cons_orbitope.h"
#include "scip/cons_pseudoboolean.h"
#include "scip/cons_setppc.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons_superindicator.h"
#include "scip/cons_symresack.h"
#include "scip/cons_varbound.h"
#include "scip/cons_xor.h"
#include "scip/cons_components.h"
#include "scip/disp_default.h"
#include "scip/dialog_default.h"
#include "scip/event_estim.h"
#include "scip/event_solvingphase.h"
#include "scip/expr_abs.h"
#include "scip/expr_entropy.h"
#include "scip/expr_exp.h"
#include "scip/expr_log.h"
#include "scip/expr_pow.h"
#include "scip/expr_product.h"
#include "scip/expr_sum.h"
#include "scip/expr_trig.h"
#include "scip/expr_value.h"
#include "scip/expr_var.h"
#include "scip/heur_actconsdiving.h"
#include "scip/heur_adaptivediving.h"
#include "scip/heur_bound.h"
#include "scip/heur_clique.h"
#include "scip/heur_coefdiving.h"
#include "scip/heur_completesol.h"
#include "scip/heur_conflictdiving.h"
#include "scip/heur_crossover.h"
#include "scip/heur_dins.h"
#include "scip/heur_distributiondiving.h"
#include "scip/heur_dps.h"
#include "scip/heur_dualval.h"
#include "scip/heur_farkasdiving.h"
#include "scip/heur_feaspump.h"
#include "scip/heur_fixandinfer.h"
#include "scip/heur_fracdiving.h"
#include "scip/heur_gins.h"
#include "scip/heur_guideddiving.h"
#include "scip/heur_indicator.h"
#include "scip/heur_indicatordiving.h"
#include "scip/heur_intdiving.h"
#include "scip/heur_intshifting.h"
#include "scip/heur_linesearchdiving.h"
#include "scip/heur_localbranching.h"
#include "scip/heur_locks.h"
#include "scip/heur_lpface.h"
#include "scip/heur_alns.h"
#include "scip/heur_multistart.h"
#include "scip/heur_mutation.h"
#include "scip/heur_mpec.h"
#include "scip/heur_nlpdiving.h"
#include "scip/heur_objpscostdiving.h"
#include "scip/heur_octane.h"
#include "scip/heur_ofins.h"
#include "scip/heur_oneopt.h"
#include "scip/heur_padm.h"
#include "scip/heur_pscostdiving.h"
#include "scip/heur_proximity.h"
#include "scip/heur_randrounding.h"
#include "scip/heur_rens.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_repair.h"
#include "scip/heur_rins.h"
#include "scip/heur_rootsoldiving.h"
#include "scip/heur_rounding.h"
#include "scip/heur_scheduler.h"
#include "scip/heur_shiftandpropagate.h"
#include "scip/heur_shifting.h"
#include "scip/heur_simplerounding.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trivial.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heur_trustregion.h"
#include "scip/heur_trysol.h"
#include "scip/heur_twoopt.h"
#include "scip/heur_undercover.h"
#include "scip/heur_vbounds.h"
#include "scip/heur_veclendiving.h"
#include "scip/heur_zeroobj.h"
#include "scip/heur_zirounding.h"
#include "scip/nlhdlr_bilinear.h"
#include "scip/nlhdlr_convex.h"
#include "scip/nlhdlr_default.h"
#include "scip/nlhdlr_perspective.h"
#include "scip/nlhdlr_quadratic.h"
#include "scip/nlhdlr_quotient.h"
#include "scip/nlhdlr_signomial.h"
#include "scip/nlhdlr_soc.h"
#include "scip/nodesel_bfs.h"
#include "scip/nodesel_breadthfirst.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_uct.h"
#include "scip/nodesel_restartdfs.h"
#include "scip/presol_boundshift.h"
#include "scip/presol_convertinttobin.h"
#include "scip/presol_domcol.h"
#include "scip/presol_dualagg.h"
#include "scip/presol_dualcomp.h"
#include "scip/presol_dualinfer.h"
#include "scip/presol_gateextraction.h"
#include "scip/presol_implics.h"
#include "scip/presol_inttobinary.h"
#include "scip/presol_milp.h"
#include "scip/presol_redvub.h"
#include "scip/presol_qpkktref.h"
#include "scip/presol_trivial.h"
#include "scip/presol_tworowbnd.h"
#include "scip/presol_sparsify.h"
#include "scip/presol_dualsparsify.h"
#include "scip/presol_stuffing.h"
#include "scip/prop_dualfix.h"
#include "scip/prop_genvbounds.h"
#include "scip/prop_nlobbt.h"
#include "scip/prop_obbt.h"
#include "scip/prop_probing.h"
#include "scip/prop_pseudoobj.h"
#include "scip/prop_redcost.h"
#include "scip/prop_rootredcost.h"
#include "scip/prop_symmetry.h"
#include "scip/prop_vbounds.h"
#include "scip/reader_bnd.h"
#include "scip/reader_ccg.h"
#include "scip/reader_cip.h"
#include "scip/reader_cnf.h"
#include "scip/reader_cor.h"
#include "scip/reader_dec.h"
#include "scip/reader_diff.h"
#include "scip/reader_fix.h"
#include "scip/reader_fzn.h"
#include "scip/reader_gms.h"
#include "scip/reader_lp.h"
#include "scip/reader_mps.h"
#include "scip/reader_mst.h"
#include "scip/reader_nl.h"
#include "scip/reader_opb.h"
#include "scip/reader_osil.h"
#include "scip/reader_pip.h"
#include "scip/reader_ppm.h"
#include "scip/reader_pbm.h"
#include "scip/reader_rlp.h"
#include "scip/reader_smps.h"
#include "scip/reader_sol.h"
#include "scip/reader_sto.h"
#include "scip/reader_tim.h"
#include "scip/reader_wbo.h"
#include "scip/reader_zpl.h"
#include "scip/sepa_eccuts.h"
#include "scip/sepa_cgmip.h"
#include "scip/sepa_clique.h"
#include "scip/sepa_closecuts.h"
#include "scip/sepa_aggregation.h"
#include "scip/sepa_convexproj.h"
#include "scip/sepa_disjunctive.h"
#include "scip/sepa_gauge.h"
#include "scip/sepa_gomory.h"
#include "scip/sepa_impliedbounds.h"
#include "scip/sepa_interminor.h"
#include "scip/sepa_intobj.h"
#include "scip/sepa_lagromory.h"
#include "scip/sepa_mcf.h"
#include "scip/sepa_minor.h"
#include "scip/sepa_mixing.h"
#include "scip/sepa_oddcycle.h"
#include "scip/sepa_rapidlearning.h"
#include "scip/sepa_rlt.h"
#include "scip/sepa_zerohalf.h"
#include "scip/scipshell.h"
#include "scip/symmetry.h"
#include "scip/table_default.h"
#include "scip/concsolver_scip.h"
#include "scip/benders_default.h"
#include "scip/cutsel_hybrid.h"
#include "scip/cutsel_dynamic.h"
#include "scip/cutsel_ensemble.h"

#include "scip/expr_varidx.h"
#include "scip/nlpi_ipopt.h"
#include "scip/nlpi_filtersqp.h"
#include "scip/nlpi_worhp.h"
#include "scip/nlpi_all.h"

/** includes default SCIP plugins into SCIP */
static
SCIP_RETCODE SCIPincludeMyDefaultPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include some default dialogs, since other plugins require that at least the root dialog is available */
   SCIP_CALL( SCIPincludeDialogDefaultBasic(scip) );

   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) ); /* nonlinear constraint handler must be before linear due to constraint upgrading */
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) ); /* linear must be before its specializations due to constraint upgrading */
   SCIP_CALL( SCIPincludeConshdlrAnd(scip) );
   SCIP_CALL( SCIPincludeConshdlrBenders(scip) );
   SCIP_CALL( SCIPincludeConshdlrBenderslp(scip) );
   SCIP_CALL( SCIPincludeConshdlrBounddisjunction(scip) );
   SCIP_CALL( SCIPincludeConshdlrCardinality(scip) );
   SCIP_CALL( SCIPincludeConshdlrConjunction(scip) );
   SCIP_CALL( SCIPincludeConshdlrCountsols(scip) );
   SCIP_CALL( SCIPincludeConshdlrCumulative(scip) );
   SCIP_CALL( SCIPincludeConshdlrDisjunction(scip) );
   SCIP_CALL( SCIPincludeConshdlrIndicator(scip) );
   SCIP_CALL( SCIPincludeConshdlrIntegral(scip) );
   SCIP_CALL( SCIPincludeConshdlrKnapsack(scip) );
   SCIP_CALL( SCIPincludeConshdlrLinking(scip) );
   SCIP_CALL( SCIPincludeConshdlrLogicor(scip) );
   SCIP_CALL( SCIPincludeConshdlrOr(scip) );
   SCIP_CALL( SCIPincludeConshdlrOrbisack(scip) );
   SCIP_CALL( SCIPincludeConshdlrOrbitope(scip) );
   SCIP_CALL( SCIPincludeConshdlrPseudoboolean(scip) );
   SCIP_CALL( SCIPincludeConshdlrSetppc(scip) );
   SCIP_CALL( SCIPincludeConshdlrSOS1(scip) );
   SCIP_CALL( SCIPincludeConshdlrSOS2(scip) );
   SCIP_CALL( SCIPincludeConshdlrSuperindicator(scip) );
   SCIP_CALL( SCIPincludeConshdlrSymresack(scip) );
   SCIP_CALL( SCIPincludeConshdlrVarbound(scip) );
   SCIP_CALL( SCIPincludeConshdlrXor(scip) );
   SCIP_CALL( SCIPincludeConshdlrComponents(scip) );

   /* include readers in order of chances to be necessary */
   SCIP_CALL( SCIPincludeReaderMps(scip) );
   SCIP_CALL( SCIPincludeReaderLp(scip) );
   SCIP_CALL( SCIPincludeReaderSol(scip) );
   SCIP_CALL( SCIPincludeReaderOsil(scip) );
   SCIP_CALL( SCIPincludeReaderZpl(scip) );
#ifdef SCIP_WITH_AMPL
   SCIP_CALL( SCIPincludeReaderNl(scip) );
#endif
   SCIP_CALL( SCIPincludeReaderGms(scip) );
   SCIP_CALL( SCIPincludeReaderOpb(scip) );
   SCIP_CALL( SCIPincludeReaderWbo(scip) );
   SCIP_CALL( SCIPincludeReaderPip(scip) );
   SCIP_CALL( SCIPincludeReaderFzn(scip) );
   SCIP_CALL( SCIPincludeReaderCnf(scip) );
   SCIP_CALL( SCIPincludeReaderCip(scip) );
   SCIP_CALL( SCIPincludeReaderSmps(scip) );
   SCIP_CALL( SCIPincludeReaderSto(scip) );
   SCIP_CALL( SCIPincludeReaderTim(scip) );
   SCIP_CALL( SCIPincludeReaderCor(scip) );
   SCIP_CALL( SCIPincludeReaderRlp(scip) );
   SCIP_CALL( SCIPincludeReaderBnd(scip) );
   SCIP_CALL( SCIPincludeReaderDiff(scip) );
   SCIP_CALL( SCIPincludeReaderDec(scip) );
   SCIP_CALL( SCIPincludeReaderFix(scip) );
   SCIP_CALL( SCIPincludeReaderMst(scip) );
   SCIP_CALL( SCIPincludeReaderPpm(scip) );
   SCIP_CALL( SCIPincludeReaderPbm(scip) );
   SCIP_CALL( SCIPincludeReaderCcg(scip) );

   SCIP_CALL( SCIPincludePresolBoundshift(scip) );
   SCIP_CALL( SCIPincludePresolConvertinttobin(scip) );
   SCIP_CALL( SCIPincludePresolDomcol(scip) );
   SCIP_CALL( SCIPincludePresolDualagg(scip) );
   SCIP_CALL( SCIPincludePresolDualcomp(scip) );
   SCIP_CALL( SCIPincludePresolDualinfer(scip) );
   SCIP_CALL( SCIPincludePresolGateextraction(scip) );
   SCIP_CALL( SCIPincludePresolImplics(scip) );
   SCIP_CALL( SCIPincludePresolInttobinary(scip) );
#ifdef SCIP_WITH_PAPILO
   SCIP_CALL( SCIPincludePresolMILP(scip) );
#endif
   SCIP_CALL( SCIPincludePresolQPKKTref(scip) );
   SCIP_CALL( SCIPincludePresolRedvub(scip) );
   SCIP_CALL( SCIPincludePresolTrivial(scip) );
   SCIP_CALL( SCIPincludePresolTworowbnd(scip) );
   SCIP_CALL( SCIPincludePresolSparsify(scip) );
   SCIP_CALL( SCIPincludePresolDualsparsify(scip) );
   SCIP_CALL( SCIPincludePresolStuffing(scip) );
   SCIP_CALL( SCIPincludeNodeselBfs(scip) );
   SCIP_CALL( SCIPincludeNodeselBreadthfirst(scip) );
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );
   SCIP_CALL( SCIPincludeNodeselEstimate(scip) );
   SCIP_CALL( SCIPincludeNodeselHybridestim(scip) );
   SCIP_CALL( SCIPincludeNodeselRestartdfs(scip) );
   SCIP_CALL( SCIPincludeNodeselUct(scip) );
   SCIP_CALL( SCIPincludeBranchruleAllfullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleCloud(scip) );
   SCIP_CALL( SCIPincludeBranchruleDistribution(scip) );
   SCIP_CALL( SCIPincludeBranchruleFullstrong(scip) );
   SCIP_CALL( SCIPincludeBranchruleGomory(scip) );
   SCIP_CALL( SCIPincludeBranchruleInference(scip) );
   SCIP_CALL( SCIPincludeBranchruleLeastinf(scip) );
   SCIP_CALL( SCIPincludeBranchruleLookahead(scip) );
   SCIP_CALL( SCIPincludeBranchruleMostinf(scip) );
   SCIP_CALL( SCIPincludeBranchruleMultAggr(scip) );
   SCIP_CALL( SCIPincludeBranchruleNodereopt(scip) );
   SCIP_CALL( SCIPincludeBranchrulePscost(scip) );
   SCIP_CALL( SCIPincludeBranchruleRandom(scip) );
   SCIP_CALL( SCIPincludeBranchruleRelpscost(scip) );
   SCIP_CALL( SCIPincludeBranchruleVanillafullstrong(scip) );
   SCIP_CALL( SCIPincludeEventHdlrEstim(scip) );
   SCIP_CALL( SCIPincludeEventHdlrSolvingphase(scip) );
   SCIP_CALL( SCIPincludeComprLargestrepr(scip) );
   SCIP_CALL( SCIPincludeComprWeakcompr(scip) );
   SCIP_CALL( SCIPincludeHeurActconsdiving(scip) );
   SCIP_CALL( SCIPincludeHeurAdaptivediving(scip) );
   SCIP_CALL( SCIPincludeHeurBound(scip) );
   SCIP_CALL( SCIPincludeHeurClique(scip) );
   SCIP_CALL( SCIPincludeHeurCoefdiving(scip) );
   SCIP_CALL( SCIPincludeHeurCompletesol(scip) );
   SCIP_CALL( SCIPincludeHeurConflictdiving(scip) );
   SCIP_CALL( SCIPincludeHeurCrossover(scip) );
   SCIP_CALL( SCIPincludeHeurDins(scip) );
   SCIP_CALL( SCIPincludeHeurDistributiondiving(scip) );
   SCIP_CALL( SCIPincludeHeurDps(scip) );
   SCIP_CALL( SCIPincludeHeurDualval(scip) );
   SCIP_CALL( SCIPincludeHeurFarkasdiving(scip) );
   SCIP_CALL( SCIPincludeHeurFeaspump(scip) );
   SCIP_CALL( SCIPincludeHeurFixandinfer(scip) );
   SCIP_CALL( SCIPincludeHeurFracdiving(scip) );
   SCIP_CALL( SCIPincludeHeurGins(scip) );
   SCIP_CALL( SCIPincludeHeurGuideddiving(scip) );
   SCIP_CALL( SCIPincludeHeurZeroobj(scip) );
   SCIP_CALL( SCIPincludeHeurIndicator(scip) );
   SCIP_CALL( SCIPincludeHeurIndicatordiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntdiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntshifting(scip) );
   SCIP_CALL( SCIPincludeHeurLinesearchdiving(scip) );
   SCIP_CALL( SCIPincludeHeurLocalbranching(scip) );
   SCIP_CALL( SCIPincludeHeurLocks(scip) );
   SCIP_CALL( SCIPincludeHeurLpface(scip) );
   SCIP_CALL( SCIPincludeHeurAlns(scip) );
   SCIP_CALL( SCIPincludeHeurNlpdiving(scip) );
   SCIP_CALL( SCIPincludeHeurMutation(scip) );
   SCIP_CALL( SCIPincludeHeurMultistart(scip) );
   SCIP_CALL( SCIPincludeHeurMpec(scip) );
   SCIP_CALL( SCIPincludeHeurObjpscostdiving(scip) );
   SCIP_CALL( SCIPincludeHeurOctane(scip) );
   SCIP_CALL( SCIPincludeHeurOfins(scip) );
   SCIP_CALL( SCIPincludeHeurOneopt(scip) );
   SCIP_CALL( SCIPincludeHeurPADM(scip) );
   SCIP_CALL( SCIPincludeHeurProximity(scip) );
   SCIP_CALL( SCIPincludeHeurPscostdiving(scip) );
   SCIP_CALL( SCIPincludeHeurRandrounding(scip) );
   SCIP_CALL( SCIPincludeHeurRens(scip) );
   SCIP_CALL( SCIPincludeHeurReoptsols(scip) );
   SCIP_CALL( SCIPincludeHeurRepair(scip) );
   SCIP_CALL( SCIPincludeHeurRins(scip) );
   SCIP_CALL( SCIPincludeHeurRootsoldiving(scip) );
   SCIP_CALL( SCIPincludeHeurRounding(scip) );
   SCIP_CALL( SCIPincludeHeurScheduler(scip) );
   SCIP_CALL( SCIPincludeHeurShiftandpropagate(scip) );
   SCIP_CALL( SCIPincludeHeurShifting(scip) );
   SCIP_CALL( SCIPincludeHeurSimplerounding(scip) );
   SCIP_CALL( SCIPincludeHeurSubNlp(scip) );
   SCIP_CALL( SCIPincludeHeurTrivial(scip) );
   SCIP_CALL( SCIPincludeHeurTrivialnegation(scip) );
   SCIP_CALL( SCIPincludeHeurTrustregion(scip) );
   SCIP_CALL( SCIPincludeHeurTrySol(scip) );
   SCIP_CALL( SCIPincludeHeurTwoopt(scip) );
   SCIP_CALL( SCIPincludeHeurUndercover(scip) );
   SCIP_CALL( SCIPincludeHeurVbounds(scip) );
   SCIP_CALL( SCIPincludeHeurVeclendiving(scip) );
   SCIP_CALL( SCIPincludeHeurZirounding(scip) );
   SCIP_CALL( SCIPincludePropDualfix(scip) );
   SCIP_CALL( SCIPincludePropGenvbounds(scip) );
   SCIP_CALL( SCIPincludePropObbt(scip) );
   SCIP_CALL( SCIPincludePropNlobbt(scip) );
   SCIP_CALL( SCIPincludePropProbing(scip) );
   SCIP_CALL( SCIPincludePropPseudoobj(scip) );
   SCIP_CALL( SCIPincludePropRedcost(scip) );
   SCIP_CALL( SCIPincludePropRootredcost(scip) );
   SCIP_CALL( SCIPincludePropVbounds(scip) );
   SCIP_CALL( SCIPincludeSepaCGMIP(scip) );
   SCIP_CALL( SCIPincludeSepaClique(scip) );
   SCIP_CALL( SCIPincludeSepaClosecuts(scip) );
   SCIP_CALL( SCIPincludeSepaAggregation(scip) );
   SCIP_CALL( SCIPincludeSepaConvexproj(scip) );
   SCIP_CALL( SCIPincludeSepaDisjunctive(scip) );
   SCIP_CALL( SCIPincludeSepaEccuts(scip) );
   SCIP_CALL( SCIPincludeSepaGauge(scip) );
   SCIP_CALL( SCIPincludeSepaGomory(scip) );
   SCIP_CALL( SCIPincludeSepaImpliedbounds(scip) );
   SCIP_CALL( SCIPincludeSepaInterminor(scip) );
   SCIP_CALL( SCIPincludeSepaIntobj(scip) );
   SCIP_CALL( SCIPincludeSepaLagromory(scip) );
   SCIP_CALL( SCIPincludeSepaMcf(scip) );
   SCIP_CALL( SCIPincludeSepaMinor(scip) );
   SCIP_CALL( SCIPincludeSepaMixing(scip) );
   SCIP_CALL( SCIPincludeSepaOddcycle(scip) );
   SCIP_CALL( SCIPincludeSepaRapidlearning(scip) );
   SCIP_CALL( SCIPincludeSepaRlt(scip) );
   SCIP_CALL( SCIPincludeSepaZerohalf(scip) );
   SCIP_CALL( SCIPincludeDispDefault(scip) );
   SCIP_CALL( SCIPincludeTableDefault(scip) );
   SCIP_CALL( SCIPincludeConcurrentScipSolvers(scip) );
   SCIP_CALL( SCIPincludeBendersDefault(scip) );
   SCIP_CALL( SCIPincludeCutselEnsemble(scip) );
   SCIP_CALL( SCIPincludeCutselHybrid(scip) );
   SCIP_CALL( SCIPincludeCutselDynamic(scip) );
   SCIP_CALL( SCIPincludeExprhdlrAbs(scip) );
   SCIP_CALL( SCIPincludeExprhdlrCos(scip) );
   SCIP_CALL( SCIPincludeExprhdlrEntropy(scip) );
   SCIP_CALL( SCIPincludeExprhdlrExp(scip) );
   SCIP_CALL( SCIPincludeExprhdlrLog(scip) );
   SCIP_CALL( SCIPincludeExprhdlrPow(scip) );
   SCIP_CALL( SCIPincludeExprhdlrProduct(scip) );
   SCIP_CALL( SCIPincludeExprhdlrSignpower(scip) );
   SCIP_CALL( SCIPincludeExprhdlrSin(scip) );
   SCIP_CALL( SCIPincludeExprhdlrSum(scip) );
   SCIP_CALL( SCIPincludeExprhdlrValue(scip) );
   SCIP_CALL( SCIPincludeExprhdlrVar(scip) );
   SCIP_CALL( SCIPincludeExprhdlrVaridx(scip) );
   SCIP_CALL( SCIPincludeNlhdlrDefault(scip) );
   SCIP_CALL( SCIPincludeNlhdlrConvex(scip) );
   SCIP_CALL( SCIPincludeNlhdlrConcave(scip) );
   SCIP_CALL( SCIPincludeNlhdlrBilinear(scip) );
   SCIP_CALL( SCIPincludeNlhdlrPerspective(scip) );
   SCIP_CALL( SCIPincludeNlhdlrQuadratic(scip) );
   SCIP_CALL( SCIPincludeNlhdlrQuotient(scip) );
   SCIP_CALL( SCIPincludeNlhdlrSignomial(scip) );
   SCIP_CALL( SCIPincludeNlhdlrSoc(scip) );
   SCIP_CALL( SCIPincludeNlpSolverIpopt(scip) );
   SCIP_CALL( SCIPincludeNlpSolverFilterSQP(scip) );
   SCIP_CALL( SCIPincludeNlpSolverWorhp(scip, TRUE) );
   SCIP_CALL( SCIPincludeNlpSolverWorhp(scip, FALSE) );
   SCIP_CALL( SCIPincludeNlpSolverAll(scip) );

   SCIP_CALL( SCIPincludePropProblemsymmetry(scip) );

#ifdef TPI_TNY
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, "TinyCThread", "Small, portable implementation of the C11 threads API (tinycthread.github.io)") );
#endif

   return SCIP_OKAY;
}

/** includes default SCIP plugins and custom plugins into SCIP */
SCIP_RETCODE SCIPincludePlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludeMyDefaultPlugins(scip) );

   return SCIP_OKAY;
}
