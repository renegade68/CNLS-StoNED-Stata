"""Research CRS CNLS-DDF: normalization elimination plus explicit distances.

Clarabel only. Full or exhaustive cutting planes, no working-set pruning.
The original unregularized objective and physical-coordinate certificate are retained.
"""
import time
import numpy as np
from scipy import sparse
import clarabel

def fit_crs(X,Y,gx,gy,*,base,tolerance,normtolerance,max_cuts,full,pivot,
            block_size,adaptive,ipm_max_iter,trace,cut_batch,inner_cap,inner_policy):
    start=time.perf_counter()
    p=base.build_problem(X,Y,gx,gy,pivot)
    n,k=p.n,p.k;q=k-1;B=p.B[:,1:];offset=p.offset
    Z=sparse.block_diag([b.reshape(1,q) for b in B],format='csc')
    Q=sparse.diags(np.r_[np.ones(n),np.zeros(n*q)],format='csc')
    E=sparse.hstack((sparse.eye(n),-Z),format='csc')
    pairs=([(i,h) for i in range(n) for h in range(n) if i!=h] if full else
           [(i,h) for i in range(n) for h in (i-1,i+1) if 0<=h<n])
    active=set(pairs);final_eps=.005*tolerance;eps=final_eps
    history=[];solver_history=[];logs=[];iterations=expansions=refinements=holds=peak=0
    converged=False;reason='maxcut'
    for cp in range(1,max_cuts+1):
        if not np.isfinite(eps) or eps<=0:raise FloatingPointError('interior tolerance underflow')
        own,alt=np.asarray(pairs,dtype=int).T;count=len(pairs)
        rows=np.repeat(np.arange(count),2*q)
        cols=np.c_[own[:,None]*q+np.arange(q),alt[:,None]*q+np.arange(q)].ravel()
        vals=np.c_[B[own],-B[own]].ravel()
        af=sparse.coo_matrix((vals,(rows,cols)),shape=(count,n*q)).tocsc()
        bound=sparse.kron(sparse.eye(n),p.a[p.free].reshape(1,q),format='csc')
        G=sparse.vstack((af,-sparse.eye(n*q),bound),format='csc');G.eliminate_zeros()
        upper=np.r_[np.zeros(count+n*q),np.ones(n)]
        A=sparse.vstack((E,sparse.hstack((sparse.csc_matrix((len(upper),n)),G))),format='csc')
        rhs=np.r_[offset,upper]
        settings=clarabel.DefaultSettings();settings.verbose=bool(trace);settings.max_threads=1;settings.max_iter=ipm_max_iter
        settings.tol_gap_abs=settings.tol_gap_rel=settings.tol_feas=eps
        # crs2: stabilize the internal KKT solve; the unregularized QP is unchanged.
        settings.static_regularization_constant=1e-6
        settings.iterative_refinement_abstol=settings.iterative_refinement_reltol=1e-14
        settings.iterative_refinement_max_iter=30
        model=clarabel.DefaultSolver(Q,np.zeros(n+n*q),A,rhs,[clarabel.ZeroConeT(n),clarabel.NonnegativeConeT(len(upper))],settings)
        steps=[]
        if trace:
            model.print_to_buffer()
            def callback(info):
                steps.append({key:float(getattr(info,key)) for key in ('iterations','res_primal','res_dual','gap_abs','gap_rel','step_length')});return False
            model.set_termination_callback(callback)
        result=model.solve();status=int(result.status);iterations+=result.iterations
        xx=np.asarray(result.x);all_dual=np.asarray(result.z);raw_distance=xx[:n];theta=xx[n:];eqdual=all_dual[:n];dual=all_dual[n:]
        base.finite(xx,all_dual)
        slopes=np.zeros((n,k));slopes[:,p.free]=theta.reshape(n,q)
        slopes[:,p.pivot]=(1-slopes[:,p.free]@p.a[p.free])/p.a[p.pivot]
        coef=np.c_[np.zeros(n),slopes/p.scales];d=np.sum(p.W*coef,axis=1)
        ax=G@theta;slack=ax-upper;lam=dual[:count];mono=dual[count:count+n*q].reshape(n,q);bd=dual[-n:]
        v=p.W[:,1:]/p.scales;grad=v*d[:,None]
        np.add.at(grad,own,lam[:,None]*v[own]);np.add.at(grad,alt,-lam[:,None]*v[own])
        eta=bd-grad[:,p.pivot]/p.a[p.pivot];nu=np.zeros((n,k));nu[:,p.free]=mono;nu[:,p.pivot]=p.a[p.pivot]*bd
        stat=grad+eta[:,None]*p.a-nu
        primal=max(float(np.max(slack,initial=0)),base.infnorm(raw_distance-d))
        dual_feas=max(float(np.max(-dual,initial=0)),float(np.max(-nu,initial=0)))
        stationarity=max(base.infnorm(stat),base.infnorm(raw_distance+eqdual))
        comp=max(base.infnorm(lam*slack[:count]),base.infnorm(nu*slopes))
        t=.05*tolerance
        ep=t*(1+max(base.infnorm(ax),base.infnorm(np.minimum(ax,upper)),base.infnorm(d)))
        ed=t*(1+max(base.infnorm(dual),base.infnorm(nu)))
        es=t*(1+max(base.infnorm(grad),base.infnorm(eta[:,None]*p.a),base.infnorm(nu),base.infnorm(d)))
        ec=t*(1+max(base.infnorm(slack[:count])*base.infnorm(lam),base.infnorm(nu)*base.infnorm(slopes)))
        violation,nviol,candidates=base.scan_all(p,coef,d,active,tolerance,block_size,cut_batch)
        diag=dict(max_viol=violation,max_mono=float(np.max(-coef[:,1:],initial=0)),max_mono_scaled=float(np.max(-slopes,initial=0)),max_normalization=base.infnorm(coef[:,1:]@p.g-1),primal=primal,eps_primal=ep,dual_feasibility=dual_feas,eps_dual_feasibility=ed,stationarity=stationarity,eps_stationarity=es,complementarity=comp,eps_complementarity=ec,rss=float(d@d),mean_distance=float(d.mean()),violating_pairs=nviol,max_equality=base.infnorm(raw_distance-d))
        base.finite(list(diag.values()))
        ok=primal<=ep and dual_feas<=ed and stationarity<=es and comp<=ec and violation<=tolerance and diag['max_mono_scaled']<=tolerance and diag['max_normalization']<=normtolerance
        native=dict(solver_status_name=str(result.status),solver_primal=float(result.r_prim),solver_dual=float(result.r_dual),objective_primal=float(result.obj_val),objective_dual=float(result.obj_val_dual),history=steps)
        base.finite([native[x] for x in ('solver_primal','solver_dual','objective_primal','objective_dual')])
        history.append(dict(iteration=cp,active=count,solver_status=status,max_violation=violation,candidates=len(candidates),inner_iters=result.iterations,inner_tol=eps))
        solver_history.append(dict(**history[-1],**native,primal=primal,stationarity=stationarity,complementarity=comp,eps_primal=ep,eps_stationarity=es,eps_complementarity=ec))
        if trace:logs.append('CRS_CLARABEL_SOLVE '+str(cp)+'\n'+model.get_print_buffer())
        peak=max(peak,sum(base.sparse_bytes(m) for m in (p.Z,p.P,Z,Q,E,G,A)))
        strict=result.status==clarabel.SolverStatus.Solved
        if strict and eps<=final_eps and ok:converged=True;reason='certified';break
        if result.status==clarabel.SolverStatus.MaxIterations:reason='inner_maxiter';break
        if result.status not in (clarabel.SolverStatus.Solved,clarabel.SolverStatus.AlmostSolved):reason='solver_failure';break
        if cp==max_cuts:break
        if candidates:
            pairs.extend(candidates);active.update(candidates);expansions+=1
            proposed=max(final_eps,min(inner_cap,.01*violation)) if adaptive else final_eps
            if inner_policy=='monotone' and cp>1:holds+=int(proposed>eps);proposed=min(eps,proposed)
            eps=proposed
        elif eps>final_eps:eps=final_eps;refinements+=1
        elif strict:eps*=.1;refinements+=1
        else:reason='solver_failure';break
    diag.update(converged=int(converged),stop_reason=reason,cp_iters=cp,n_active=len(pairs),total_inner_iters=iterations,ipm_iterations=iterations,ipm_maxiter=ipm_max_iter,cut_batch=cut_batch,inner_cap=inner_cap,inner_monotone=int(inner_policy=='monotone'),tolerance_holds=holds,osqp_iterations=0,setups=cp,same_qp_resumes=0,expanded_primal_starts=0,constraint_expansions=expansions,strict_refinements=refinements,solver_status=status,inner_tolerance=history[-1]['inner_tol'],matrix_payload_bytes=peak,ipm_matrix_bytes=peak,dense_normal_matrix_bytes=8*(n*k)**2,scan_workspace_bytes=8*min(block_size,n)*n,parameters=n*k,free_slopes=n*q,equalities=n,seconds=time.perf_counter()-start,pivot=p.pivot+1,N=n,m=p.X.shape[1],s=p.Y.shape[1],tolerance=tolerance,normtolerance=normtolerance,objective_ridge=0.,initializer_used=0,initializer_status=0,initializer_iterations=0,initializer_seconds=0.,initializer_matrix_bytes=0,initializer_tolerance=0.,status_resets=0,rts_crs=1)
    diag.update(kkt_static_regularization=1e-6,kkt_refinement_abstol=1e-14,kkt_refinement_reltol=1e-14,kkt_refinement_maxiter=30,kkt_profile='crs2_reg6_refine')
    return base._working_set_metadata(dict(distance=d,residual=-d,coefficients=coef,theta=theta,dual=dual,equality_dual=eqdual,raw_distance=raw_distance,pairs=np.asarray(pairs),lam=lam,nu=nu,eta=eta,scales=p.scales,diagnostics=diag,history=history,solver_history=solver_history,initialization=None,trace='\n'.join(logs)),'none')
