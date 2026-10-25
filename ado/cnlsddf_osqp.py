"""Unregularized sparse VRS CNLS-DDF QP. Python API and file bridge for Stata.

Pairs are (own observation, alternative plane). Every Afriat row is
B_i theta_i - B_i theta_h <= 0; all final checks use original distance units.
OSQP's internal proximal sigma is not a ridge in the estimator's objective.
"""
from __future__ import annotations
import argparse
import contextlib
import csv
from dataclasses import dataclass
import io
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy import sparse
import scipy
import osqp

VERSION = '0.1.0-crs2'


def finite(*arrays):
    if not all(np.all(np.isfinite(a)) for a in arrays):
        raise FloatingPointError('nonfinite data, iterate, or diagnostic')


def infnorm(x):
    finite(x)
    return float(np.max(np.abs(x), initial=0.0))


def sparse_bytes(a):
    return a.data.nbytes + a.indices.nbytes + a.indptr.nbytes


@dataclass
class Problem:
    X: np.ndarray
    Y: np.ndarray
    g: np.ndarray
    scales: np.ndarray
    a: np.ndarray
    pivot: int
    free: np.ndarray
    B: np.ndarray
    offset: np.ndarray
    Z: sparse.csc_matrix
    P: sparse.csc_matrix
    linear: np.ndarray
    W: np.ndarray

    @property
    def n(self): return self.X.shape[0]

    @property
    def k(self): return self.g.size

    def coefficients(self, theta):
        t = np.asarray(theta).reshape(self.n, self.k)
        slopes = np.zeros((self.n, self.k))
        slopes[:, self.free] = t[:, 1:]
        slopes[:, self.pivot] = (1 - slopes[:, self.free] @ self.a[self.free]) / self.a[self.pivot]
        coef = np.column_stack((t[:, 0], slopes / self.scales))
        distance = np.sum(self.W * coef, axis=1)
        finite(coef, slopes, distance)
        return coef, slopes, distance


def build_problem(X, Y, gx, gy, pivot=None):
    X, Y = np.asarray(X, dtype=float), np.asarray(Y, dtype=float)
    gx, gy = np.asarray(gx, dtype=float), np.asarray(gy, dtype=float)
    if X.ndim != 2 or Y.ndim != 2 or X.shape[0] != Y.shape[0] or X.shape[0] < 2 or min(X.shape[1], Y.shape[1]) < 1:
        raise ValueError('matching matrices with N>=2 and nonempty input/output blocks required')
    if gx.shape != (X.shape[1],) or gy.shape != (Y.shape[1],):
        raise ValueError('direction dimensions do not match quantities')
    finite(X, Y, gx, gy)
    g = np.r_[gx, gy]
    if np.any(X < 0) or np.any(Y < 0) or np.any(g < 0) or not np.any(g > 0):
        raise ValueError('nonnegative quantities and a nonzero nonnegative direction required')
    V = np.column_stack((X, -Y))
    scales = np.max(np.abs(V), axis=0)
    scales[scales == 0] = 1
    a, Vs = g / scales, V / scales
    if pivot is None: pivot = int(np.argmax(a))
    if not isinstance(pivot, (int, np.integer)) or not 0 <= pivot < g.size or a[pivot] <= 0:
        raise ValueError('pivot must select a positive scaled direction component')
    free = np.delete(np.arange(g.size), pivot)
    B = np.column_stack((np.ones(len(X)), Vs[:, free] - np.outer(Vs[:, pivot], a[free] / a[pivot])))
    offset = Vs[:, pivot] / a[pivot]
    finite(a, B, offset)
    n, k = B.shape
    Z = sparse.coo_matrix((B.ravel(), (np.repeat(np.arange(n), k), np.arange(n*k))), shape=(n, n*k)).tocsc()
    P = (Z.T @ Z).tocsc()  # PSD, exactly the unregularized least-squares objective.
    linear = np.asarray(Z.T @ offset).ravel()
    finite(P.data, linear)
    return Problem(X, Y, g, scales, a, pivot, free, B, offset, Z, P, linear, np.column_stack((np.ones(n), V)))


def constraints(p, pairs):
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    count = len(pairs); n, k = p.n, p.k
    if count:
        if np.any(pairs < 0) or np.any(pairs >= n) or np.any(pairs[:, 0] == pairs[:, 1]):
            raise ValueError('invalid Afriat pair')
        own, alt = pairs.T
        rows = np.repeat(np.arange(count), 2*k)
        cols = np.column_stack((own[:, None]*k + np.arange(k), alt[:, None]*k + np.arange(k))).ravel()
        vals = np.column_stack((p.B[own], -p.B[own])).ravel()
        afriat = sparse.coo_matrix((vals, (rows, cols)), shape=(count, n*k)).tocsc()
    else: afriat = sparse.csc_matrix((0, n*k))
    cols = (np.arange(n)[:, None]*k + np.arange(1, k)).ravel()
    mono = sparse.coo_matrix((-np.ones(n*(k-1)), (np.arange(n*(k-1)), cols)), shape=(n*(k-1), n*k)).tocsc()
    bound = sparse.coo_matrix((np.tile(p.a[p.free], n), (np.repeat(np.arange(n), k-1), cols)), shape=(n, n*k)).tocsc()
    A = sparse.vstack((afriat, mono, bound), format='csc')
    A.eliminate_zeros()
    upper = np.r_[np.zeros(count+n*(k-1)), np.ones(n)]
    return A, upper


def scan_all(p, coef, distance, active, tol, block_size=128, cut_batch=1):
    """Scan EVERY ordered pair; keep at most cut_batch new pairs per observation.

    The maximum and count never depend on the candidate-storage limit. Peak
    temporary pair storage is block_size*N, not N*N. No early truncation.
    """
    maximum, nviol, candidates = 0.0, 0, []
    alternatives = {}
    for i, h in active: alternatives.setdefault(i, []).append(h)
    for start in range(0, p.n, block_size):
        end = min(start + block_size, p.n)
        raw = distance[start:end, None] - p.W[start:end] @ coef.T
        finite(raw)
        raw[np.arange(end-start), np.arange(start, end)] = 0
        maximum = max(maximum, float(np.max(raw, initial=0)))
        nviol += int(np.count_nonzero(raw > tol))
        for j, i in enumerate(range(start, end)):
            raw[j, i] = -np.inf
            if i in alternatives: raw[j, alternatives[i]] = -np.inf
            for _ in range(min(cut_batch, p.n-1)):
                h = int(np.argmax(raw[j]))
                if raw[j, h] <= tol: break
                candidates.append((i, h))
                raw[j, h] = -np.inf
    return maximum, nviol, candidates


def certificate(p, A, upper, theta, dual, tol, normtol, active, block_size=128, cut_batch=1):
    finite(theta, dual)
    coef, slopes, d = p.coefficients(theta)
    ax, px, ay = A @ theta, p.P @ theta, A.T @ dual
    finite(ax, px, ay)
    # Independent KKT checks, unscaled by OSQP's internal equilibration.
    t = .05 * tol
    projected = np.minimum(ax, upper)
    primal = float(np.max(ax-upper, initial=0))
    dual_feas = float(np.max(-dual, initial=0))
    stationarity = infnorm(px + p.linear + ay)
    comp = infnorm((ax-upper)*dual)
    eps_primal = t*(1+max(infnorm(ax), infnorm(projected)))
    eps_dual = t*(1+infnorm(dual))
    eps_stat = t*(1+max(infnorm(px), infnorm(p.linear), infnorm(ay)))
    eps_comp = t*(1+infnorm(ax-upper)*infnorm(dual))
    violation, nviol, candidates = scan_all(p, coef, d, active, tol, block_size, cut_batch)
    diag = dict(max_viol=violation, max_mono=float(np.max(-coef[:, 1:], initial=0)),
        max_mono_scaled=float(np.max(-slopes, initial=0)),
        max_normalization=infnorm(coef[:, 1:] @ p.g - 1),
        primal=primal, eps_primal=eps_primal, dual_feasibility=dual_feas, eps_dual_feasibility=eps_dual,
        stationarity=stationarity, eps_stationarity=eps_stat, complementarity=comp, eps_complementarity=eps_comp,
        rss=float(d @ d), mean_distance=float(np.mean(d)), violating_pairs=nviol)
    finite(list(diag.values()))
    ok = (primal <= eps_primal and dual_feas <= eps_dual and stationarity <= eps_stat and comp <= eps_comp
          and violation <= tol and diag['max_mono_scaled'] <= tol and diag['max_normalization'] <= normtol)
    return ok, diag, candidates, coef, d


def distance_coordinates(p, A, upper):
    """Exact affine change theta=T*z+c, where z_i=(d_i, free slopes).

    Z*T selects distances and Z*c+offset=0, hence the new Hessian is
    diag(1,0,...,0) per observation and the linear term is exactly zero.
    Multipliers are unchanged: no constraint scaling or pruning is applied.
    """
    n, k = p.n, p.k
    T = sparse.eye(n*k, format='csc') + sparse.coo_matrix(
        (-p.B[:, 1:].ravel(), (np.repeat(np.arange(n)*k, k-1),
         (np.arange(n)[:, None]*k+np.arange(1, k)).ravel())), shape=(n*k, n*k)).tocsc()
    shift = np.zeros(n*k); shift[::k] = -p.offset
    Q = sparse.diags(np.tile(np.r_[1., np.zeros(k-1)], n), format='csc')
    G = (A @ T).tocsc(); G.eliminate_zeros()
    bound = upper-A @ shift
    finite(T.data, shift, Q.data, G.data, bound)
    return T, shift, Q, G, bound


def interior_start(p, A, upper, active, tolerance, normtolerance, max_iter, block_size, trace):
    """Solve the SAME full QP for a starting point; never authorize final success.

    Only a strict Solved result passing the original-coordinate certificate is
    used. A finite unsuccessful attempt falls back to the original cold start.
    Dependency/setup errors propagate rather than silently skipping the method.
    """
    import clarabel
    start = time.perf_counter()
    T, shift, Q, G, bound = distance_coordinates(p, A, upper)
    settings = clarabel.DefaultSettings()
    settings.verbose = bool(trace)
    settings.max_threads = 1
    settings.max_iter = min(max_iter, 100)
    # A more accurate interior start prevents ADMM from moving away from a
    # nearly optimal but weakly identified face. Final OSQP/KKT criteria below
    # are unchanged; this is a stricter initializer, not relaxed certification.
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = .00005*tolerance
    solver = clarabel.DefaultSolver(Q, np.zeros(p.n*p.k), G, bound,
        [clarabel.NonnegativeConeT(len(upper))], settings)
    steps = []
    if trace:
        solver.print_to_buffer()
        def snapshot(info):
            record = {key:float(getattr(info, key)) for key in
                ('iterations','res_primal','res_dual','gap_abs','gap_rel','mu','step_length')}
            # Initial diagnostic values may be infinite; they are not final iterates.
            steps.append({key:value if np.isfinite(value) else None for key,value in record.items()})
            return False
        solver.set_termination_callback(snapshot)
    result = solver.solve()
    theta, dual = np.asarray(T @ np.asarray(result.x)+shift), np.asarray(result.z)
    finite(theta, dual)
    certified, diag, *_ = certificate(p, A, upper, theta, dual, tolerance, normtolerance, active, block_size)
    used = result.status == clarabel.SolverStatus.Solved and certified
    details = dict(version=clarabel.__version__, status=str(result.status),
        status_code=int(result.status), certificate=diag, iterations=int(result.iterations),
        used=bool(used), history=steps, log=solver.get_print_buffer() if trace else '')
    summary = dict(initializer_used=int(used), initializer_status=int(result.status),
        initializer_iterations=int(result.iterations), initializer_seconds=time.perf_counter()-start,
        initializer_matrix_bytes=sum(sparse_bytes(m) for m in (T,Q,G)),
        initializer_tolerance=settings.tol_feas)
    return theta, dual, summary, details


def clarabel_qp(p, A, upper, eps, max_iter, trace):
    """Fresh interior solve of the exact restricted QP, with original-row duals.

    A new model is constructed each call. No warm starts or factor reuse are
    claimed. Clarabel's internal KKT regularization is not an objective ridge.
    """
    import clarabel
    T, shift, Q, G, bound = distance_coordinates(p, A, upper)
    settings = clarabel.DefaultSettings()
    settings.verbose = bool(trace)
    settings.max_threads = 1
    settings.max_iter = max_iter
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = eps
    qp = clarabel.DefaultSolver(Q, np.zeros(p.n*p.k), G, bound,
        [clarabel.NonnegativeConeT(len(upper))], settings)
    steps = []
    if trace:
        qp.print_to_buffer()
        def snapshot(info):
            values = {key:float(getattr(info, key)) for key in
                ('iterations','res_primal','res_dual','gap_abs','gap_rel','mu','step_length')}
            steps.append({key:value if np.isfinite(value) else None for key,value in values.items()})
            return False
        qp.set_termination_callback(snapshot)
    result = qp.solve()
    theta, dual = np.asarray(T @ np.asarray(result.x)+shift), np.asarray(result.z)
    finite(theta, dual)
    native = dict(status=int(result.status), status_name=str(result.status),
        iterations=int(result.iterations), solver_primal=float(result.r_prim),
        solver_dual=float(result.r_dual), objective_primal=float(result.obj_val),
        objective_dual=float(result.obj_val_dual), history=steps,
        matrix_bytes=sum(sparse_bytes(m) for m in (T,Q,G)),
        log=qp.get_print_buffer() if trace else '')
    finite([native[key] for key in ('solver_primal','solver_dual','objective_primal','objective_dual')])
    return theta, dual, native


def fit_clarabel(p, *, tolerance, normtolerance, max_cuts, full, block_size,
                 adaptive, ipm_max_iter, trace, start, cut_batch, inner_cap, inner_policy):
    """Cutting planes (or full reference solve) with strict interior certification."""
    import clarabel
    pairs = ([(i,h) for i in range(p.n) for h in range(p.n) if i != h] if full else
             [(i,h) for i in range(p.n) for h in (i-1,i+1) if 0 <= h < p.n])
    active = set(pairs)
    final_eps = .005*tolerance
    eps = final_eps
    history, solver_history, logs = [], [], []
    iterations = expansions = refinements = payload = ipm_payload = 0
    tolerance_holds = 0
    converged, reason = False, 'maxcut'
    for cp in range(1,max_cuts+1):
        A, upper = constraints(p,pairs)
        if not np.isfinite(eps) or eps <= 0:
            raise FloatingPointError('interior tolerance underflow')
        theta, dual, native = clarabel_qp(p,A,upper,eps,ipm_max_iter,trace)
        iterations += native['iterations']
        payload = max(payload,sum(sparse_bytes(m) for m in (p.Z,p.P,A)))
        ipm_payload = max(ipm_payload,native['matrix_bytes'])
        ok, diag, candidates, coef, d = certificate(p,A,upper,theta,dual,
            tolerance,normtolerance,active,block_size,cut_batch)
        history.append(dict(iteration=cp,active=len(pairs),solver_status=native['status'],
            max_violation=diag['max_viol'],candidates=len(candidates),
            inner_iters=native['iterations'],inner_tol=eps))
        solver_history.append(dict(**history[-1],solver_status_name=native['status_name'],
            **{key:native[key] for key in ('solver_primal','solver_dual','objective_primal','objective_dual','history')},
            **{key:diag[key] for key in ('primal','stationarity','complementarity',
                'eps_primal','eps_stationarity','eps_complementarity')}))
        if trace: logs.append('CLARABEL_SOLVE '+str(cp)+'\n'+native['log'])
        strict = native['status'] == int(clarabel.SolverStatus.Solved)
        if strict and eps <= final_eps and ok:
            converged, reason = True, 'certified'
            break
        if native['status'] == int(clarabel.SolverStatus.MaxIterations):
            reason = 'inner_maxiter'
            break
        if native['status'] not in (int(clarabel.SolverStatus.Solved),int(clarabel.SolverStatus.AlmostSolved)):
            reason = 'solver_failure'
            break
        if cp == max_cuts: break
        if candidates:
            pairs.extend(candidates); active.update(candidates); expansions += 1
            proposed = max(final_eps,min(inner_cap,.01*diag['max_viol'])) if adaptive else final_eps
            # The first QP deliberately starts strict. Allow its initial release,
            # then retain every tightening from QP 2 onward, even after new cuts.
            # State is local to this fit; no IPM iterate or factor is reused.
            if inner_policy == 'monotone' and cp > 1:
                tolerance_holds += int(proposed > eps)
                proposed = min(eps, proposed)
            eps = proposed
        elif eps > final_eps:
            eps = final_eps; refinements += 1
        elif strict:
            eps *= .1; refinements += 1
        else:
            # AlmostSolved is never final success, even when the external
            # certificate passes. Fresh solves do not continue an IPM iterate.
            reason = 'solver_failure'
            break
    diag.update(converged=int(converged),stop_reason=reason,cp_iters=cp,n_active=len(pairs),
        total_inner_iters=iterations,ipm_iterations=iterations,ipm_maxiter=ipm_max_iter,
        cut_batch=cut_batch,inner_cap=inner_cap,
        inner_monotone=int(inner_policy == 'monotone'),tolerance_holds=tolerance_holds,
        osqp_iterations=0,setups=cp,same_qp_resumes=0,expanded_primal_starts=0,
        constraint_expansions=expansions,strict_refinements=refinements,
        solver_status=native['status'],inner_tolerance=history[-1]['inner_tol'],
        matrix_payload_bytes=payload,ipm_matrix_bytes=ipm_payload,
        dense_normal_matrix_bytes=8*(p.n*p.k)**2,scan_workspace_bytes=8*min(block_size,p.n)*p.n,
        parameters=p.n*p.k,seconds=time.perf_counter()-start,pivot=p.pivot+1,
        N=p.n,m=p.X.shape[1],s=p.Y.shape[1],tolerance=tolerance,normtolerance=normtolerance,
        objective_ridge=0.,initializer_used=0,initializer_status=0,initializer_iterations=0,
        initializer_seconds=0.,initializer_matrix_bytes=0,initializer_tolerance=0.,status_resets=0)
    return dict(distance=d,residual=-d,coefficients=coef,diagnostics=diag,history=history,
        theta=theta,dual=dual,pairs=np.asarray(pairs),solver_history=solver_history,
        initialization=None,trace='\n'.join(logs))


_WORKING_SET_POLICIES = {'none': (), 'once': (10,), 'bounded': (5,10,15,20)}


class _WorkingSet:
    def __init__(self, policy, initial):
        if policy not in _WORKING_SET_POLICIES: raise ValueError('unknown working-set policy')
        self.schedule = _WORKING_SET_POLICIES[policy]
        self.protected = set(initial)
        self.removed = set()
        self.pinned = set()
        self.events = self.removals = self.readmissions = 0

    def select(self, cp, pairs, A, upper, theta, dual, distance, eps, tol, strict):
        if cp not in self.schedule or not strict or eps <= .005*tol: return set()
        # Afriat slacks and multipliers in ORIGINAL distance-row coordinates.
        n = len(pairs)
        slack = np.asarray(upper[:n] - A[:n] @ theta)
        own = np.asarray([distance[i] for i,h in pairs])
        margin = 10*max(tol,eps)*(1+np.maximum(np.abs(own),np.abs(own+slack)))
        dual_limit = max(tol,eps)*(1+float(np.max(np.abs(distance),initial=0)))
        eligible = [j for j,pair in enumerate(pairs) if pair not in self.protected
            and pair not in self.removed and slack[j] > margin[j] and abs(dual[j]) <= dual_limit]
        # Keep at least half the working pairs; deterministic slack/index order.
        ranked = sorted(eligible,key=lambda j:(-slack[j],j))[:n//2]
        chosen = set(ranked)
        if chosen:
            self.events += 1;self.removals += len(chosen)
            self.removed.update(pairs[j] for j in chosen)
        return chosen

    def added(self, candidates):
        returning = set(candidates) & self.removed
        self.pinned.update(returning)
        self.readmissions += len(returning)
        return len(returning)


def _fit_clarabel_working_set(p, *, tolerance, normtolerance, max_cuts, full, block_size,
                 adaptive, ipm_max_iter, trace, start, cut_batch, inner_cap, inner_policy, policy='bounded'):
    """Cutting planes (or full reference solve) with strict interior certification."""
    import clarabel
    if full and policy != 'none': raise ValueError('pruning requires cutting planes')
    pairs = ([(i,h) for i in range(p.n) for h in range(p.n) if i != h] if full else
             [(i,h) for i in range(p.n) for h in (i-1,i+1) if 0 <= h < p.n])
    active = set(pairs)
    manager = _WorkingSet(policy, pairs)
    working_history = []
    recovery_history = []
    final_eps = .005*tolerance
    eps = final_eps
    history, solver_history, logs = [], [], []
    iterations = expansions = refinements = payload = ipm_payload = 0
    tolerance_holds = 0
    converged, reason = False, 'maxcut'
    for cp in range(1,max_cuts+1):
        A, upper = constraints(p,pairs)
        if not np.isfinite(eps) or eps <= 0:
            raise FloatingPointError('interior tolerance underflow')
        theta, dual, native = clarabel_qp(p,A,upper,eps,ipm_max_iter,trace)
        iterations += native['iterations']
        payload = max(payload,sum(sparse_bytes(m) for m in (p.Z,p.P,A)))
        ipm_payload = max(ipm_payload,native['matrix_bytes'])
        ok, diag, candidates, coef, d = certificate(p,A,upper,theta,dual,
            tolerance,normtolerance,active,block_size,cut_batch)
        history.append(dict(iteration=cp,active=len(pairs),solver_status=native['status'],
            max_violation=diag['max_viol'],candidates=len(candidates),
            inner_iters=native['iterations'],inner_tol=eps))
        solver_history.append(dict(**history[-1],solver_status_name=native['status_name'],
            **{key:native[key] for key in ('solver_primal','solver_dual','objective_primal','objective_dual','history')},
            **{key:diag[key] for key in ('primal','stationarity','complementarity',
                'eps_primal','eps_stationarity','eps_complementarity')}))
        if trace: logs.append('CLARABEL_SOLVE '+str(cp)+'\n'+native['log'])
        strict = native['status'] == int(clarabel.SolverStatus.Solved)
        if strict and eps <= final_eps and ok:
            converged, reason = True, 'certified'
            break
        # Only a terminal native numerical failure triggers restoration. Budget
        # exhaustion, infeasibility, exceptions/nonfinite output, or a Solved QP
        # failing its external certificate do not acquire an extra attempt.
        terminal_almost = (native['status'] == int(clarabel.SolverStatus.AlmostSolved)
                           and not candidates and eps <= final_eps)
        numerical_failure = native['status'] in (
            int(clarabel.SolverStatus.NumericalError),
            int(clarabel.SolverStatus.InsufficientProgress))
        if not recovery_history and cp < max_cuts and (terminal_almost or numerical_failure):
            restored = sorted(manager.removed - active)
            if restored:
                # Keep solved rows in their existing order and append missing
                # previously removed rows deterministically. The next loop call
                # consumes the original maxcut/IPM budget; never reset either.
                recovery_history.append(dict(iteration=cp,trigger=native['status_name'],
                    solver_status=native['status'],**({'restored_pairs':restored} if trace else {}),
                    restored=len(restored),next_active=len(pairs)+len(restored),
                    inner_tol=min(eps,final_eps)))
                pairs.extend(restored);active.update(restored)
                manager.schedule = ()  # permanently disable further pruning
                eps = min(eps,final_eps)
                continue
        if native['status'] == int(clarabel.SolverStatus.MaxIterations):
            reason = 'inner_maxiter'
            break
        if native['status'] not in (int(clarabel.SolverStatus.Solved),int(clarabel.SolverStatus.AlmostSolved)):
            reason = 'solver_failure'
            break
        if cp == max_cuts: break
        if candidates:
            # Only alter the NEXT QP. Current theta/dual/certificate/history still
            # refer to the exact solved rows, including on every early exit.
            drop = manager.select(cp, pairs, A, upper, theta, dual, d, eps, tolerance, strict)
            if trace:
                removed_pairs = [pair for j,pair in enumerate(pairs) if j in drop]
                readmitted_pairs = [pair for pair in candidates if pair in manager.removed]
            readmitted = manager.added(candidates)
            if drop:
                pairs = [pair for j, pair in enumerate(pairs) if j not in drop]
                active = set(pairs)
            if trace:
                working_history.append(dict(iteration=cp,removed=len(drop),readmitted=readmitted,
                    pinned=len(manager.pinned),next_active=len(pairs)+len(candidates),
                    removed_pairs=removed_pairs,added_pairs=list(candidates),readmitted_pairs=readmitted_pairs))
            pairs.extend(candidates); active.update(candidates); expansions += 1
            proposed = max(final_eps,min(inner_cap,.01*diag['max_viol'])) if adaptive else final_eps
            # The first QP deliberately starts strict. Allow its initial release,
            # then retain every tightening from QP 2 onward, even after new cuts.
            # State is local to this fit; no IPM iterate or factor is reused.
            if inner_policy == 'monotone' and cp > 1:
                tolerance_holds += int(proposed > eps)
                proposed = min(eps, proposed)
            # Once recovered, every subsequent fresh QP remains strict.
            eps = min(eps, final_eps, proposed) if recovery_history else proposed
        elif eps > final_eps:
            eps = final_eps; refinements += 1
        elif strict:
            eps *= .1; refinements += 1
        else:
            # AlmostSolved is never final success, even when the external
            # certificate passes. Fresh solves do not continue an IPM iterate.
            reason = 'solver_failure'
            break
    diag.update(converged=int(converged),stop_reason=reason,cp_iters=cp,n_active=len(pairs),
        total_inner_iters=iterations,ipm_iterations=iterations,ipm_maxiter=ipm_max_iter,
        cut_batch=cut_batch,inner_cap=inner_cap,
        recovery_enabled=1,recovery_attempts=len(recovery_history),
        recovery_restored=sum(h['restored'] for h in recovery_history),
        working_set_prunes=manager.events,working_set_removed=manager.removals,
        working_set_readmitted=manager.readmissions,working_set_pinned=len(manager.pinned),
        inner_monotone=int(inner_policy == 'monotone'),tolerance_holds=tolerance_holds,
        osqp_iterations=0,setups=cp,same_qp_resumes=0,expanded_primal_starts=0,
        constraint_expansions=expansions,strict_refinements=refinements,
        solver_status=native['status'],inner_tolerance=history[-1]['inner_tol'],
        matrix_payload_bytes=payload,ipm_matrix_bytes=ipm_payload,
        dense_normal_matrix_bytes=8*(p.n*p.k)**2,scan_workspace_bytes=8*min(block_size,p.n)*p.n,
        parameters=p.n*p.k,seconds=time.perf_counter()-start,pivot=p.pivot+1,
        N=p.n,m=p.X.shape[1],s=p.Y.shape[1],tolerance=tolerance,normtolerance=normtolerance,
        objective_ridge=0.,initializer_used=0,initializer_status=0,initializer_iterations=0,
        initializer_seconds=0.,initializer_matrix_bytes=0,initializer_tolerance=0.,status_resets=0)
    return dict(distance=d,residual=-d,coefficients=coef,diagnostics=diag,history=history,
        theta=theta,dual=dual,pairs=np.asarray(pairs),solver_history=solver_history,
        initialization=None,trace='\n'.join(logs),working_set_history=working_history,
        recovery_history=recovery_history)


def _working_set_metadata(result, policy):
    # Keep the default solver loops unchanged; expose zero counters consistently
    # for OSQP/full/no-pruning estimates and clear state on every invocation.
    result['working_set'] = policy
    result['diagnostics']['working_set_enabled'] = int(policy != 'none')
    for name in ('working_set_prunes','working_set_removed','working_set_readmitted',
                 'working_set_pinned','recovery_enabled','recovery_attempts','recovery_restored'):
        result['diagnostics'].setdefault(name, 0)
    result.setdefault('working_set_history', [])
    result.setdefault('recovery_history', [])
    return result


def fit_ddf(X, Y, gx, gy, *, tolerance=1e-6, normtolerance=1e-8, max_iter=200000,
            max_cuts=100, full=False, pivot=None, block_size=128, adaptive=True,
            full_start='auto', trace=False, solver='osqp', ipm_max_iter=100,
            cut_batch=1, inner_cap=1e-3, inner_policy='adaptive', working_set='none', rts='vrs'):
    for value in (tolerance, normtolerance):
        if not np.isfinite(value) or value <= 0: raise ValueError('positive finite tolerances required')
    for value in (max_iter, max_cuts, block_size, ipm_max_iter, cut_batch):
        if not isinstance(value, (int, np.integer)) or value < 1: raise ValueError('positive integer budgets required')
    if full_start not in ('auto','cold','interior') or (full_start=='interior' and not full):
        raise ValueError('full_start must be auto/cold, or interior with full=True')
    if solver not in ('osqp','clarabel'):
        raise ValueError('solver must be osqp or clarabel')
    if solver == 'clarabel' and full_start != 'auto':
        raise ValueError('full_start applies only to the OSQP engine')
    if not np.isfinite(inner_cap) or not 0 < inner_cap <= 1e-3:
        raise ValueError('inner_cap must be finite and in (0, 1e-3]')
    if inner_policy not in ('adaptive','monotone'):
        raise ValueError('inner_policy must be adaptive or monotone')
    if (cut_batch != 1 or inner_cap != 1e-3 or inner_policy != 'adaptive') and (solver != 'clarabel' or full):
        raise ValueError('nondefault cut_batch/inner_cap/inner_policy require Clarabel cutting planes')
    if working_set not in _WORKING_SET_POLICIES:
        raise ValueError('working_set must be none, once or bounded')
    if working_set != 'none' and (solver != 'clarabel' or full):
        raise ValueError('working_set pruning requires Clarabel cutting planes')
    if rts not in ('vrs','crs'):
        raise ValueError('rts must be vrs or crs')
    if rts == 'crs':
        if solver != 'clarabel' or working_set != 'none':
            raise ValueError('CRS research candidate requires solver=clarabel and working_set=none')
        from cnlsddf_crs import fit_crs
        return fit_crs(X,Y,gx,gy,base=sys.modules[__name__],tolerance=tolerance,
            normtolerance=normtolerance,max_cuts=max_cuts,full=full,pivot=pivot,
            block_size=block_size,adaptive=adaptive,ipm_max_iter=ipm_max_iter,
            trace=trace,cut_batch=cut_batch,inner_cap=inner_cap,inner_policy=inner_policy)
    use_interior = full and full_start != 'cold'
    start = time.perf_counter()
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        p = build_problem(X, Y, gx, gy, pivot)
        if solver == 'clarabel':
            driver = fit_clarabel if working_set == 'none' else _fit_clarabel_working_set
            options = {} if working_set == 'none' else dict(policy=working_set)
            result = driver(p,tolerance=tolerance,normtolerance=normtolerance,
                max_cuts=max_cuts,full=full,block_size=block_size,adaptive=adaptive,
                ipm_max_iter=ipm_max_iter,trace=trace,start=start,cut_batch=cut_batch,inner_cap=inner_cap,inner_policy=inner_policy,**options)
            return _working_set_metadata(result, working_set)
        pairs = ([(i, h) for i in range(p.n) for h in range(p.n) if i != h] if full
                 else [(i, h) for i in range(p.n) for h in (i-1, i+1) if 0 <= h < p.n])
        active = set(pairs)
        theta = np.zeros(p.n*p.k)
        theta[::p.k] = -np.mean(p.offset)
        history, solver, A, upper = [], None, None, None
        # A tenfold numerical margin below the independent KKT threshold helps
        # resolve fitted distances in nearly flat restricted QPs.
        final_eps = .005*tolerance
        eps = final_eps
        iterations = setups = resumes = primal_starts = 0
        status_resets = 0
        init_summary = dict(initializer_used=0, initializer_status=0, initializer_iterations=0,
                            initializer_seconds=0., initializer_matrix_bytes=0, initializer_tolerance=0.)
        initialization = None
        solver_history = []
        log = io.StringIO()
        payload_peak = 0
        reason, converged = 'maxcut', False
        for cp in range(1, max_cuts+1):
            if solver is None:
                A, upper = constraints(p, pairs)
                initial_dual = np.zeros(len(upper))
                if cp == 1 and use_interior:
                    initial_theta, candidate_dual, init_summary, initialization = interior_start(
                        p,A,upper,active,tolerance,normtolerance,max_iter,block_size,trace)
                    if init_summary['initializer_used']:
                        theta, initial_dual = initial_theta, candidate_dual
                solver = osqp.OSQP()
                with contextlib.redirect_stdout(log) if trace else contextlib.nullcontext():
                    solver.setup(P=sparse.triu(p.P, format='csc'), q=p.linear, A=A,
                        l=np.full(len(upper), -np.inf), u=upper, verbose=bool(trace),
                        eps_abs=eps, eps_rel=eps, max_iter=max_iter, polishing=True,
                        warm_starting=True, adaptive_rho_interval=50, scaled_termination=False)
                # Changed constraint set: primal warm start only; all duals reset.
                solver.warm_start(x=theta, y=initial_dual)
                setups += 1
                primal_starts += int(cp > 1)
                payload_peak = max(payload_peak, sum(sparse_bytes(m) for m in (p.Z, p.P, A)))
            else:
                solver.update_settings(eps_abs=eps, eps_rel=eps)
                # OSQP 1.0.4's bundled core does not reset status in solve() or
                # update_settings(). An identical q update resets info without
                # changing the QP, factoring a new matrix, or dropping warm state.
                solver.update(q=p.linear)
                status_resets += 1
                resumes += 1
            with contextlib.redirect_stdout(log) if trace else contextlib.nullcontext():
                if trace: print('DDF_SOLVE',cp,'eps',eps,'active',len(pairs))
                res = solver.solve(raise_error=False)
            iterations += int(res.info.iter)
            theta = res.x
            ok, diag, candidates, coef, d = certificate(p, A, upper, theta, res.y,
                tolerance, normtolerance, active, block_size)
            history.append(dict(iteration=cp, active=len(pairs), solver_status=int(res.info.status_val),
                max_violation=diag['max_viol'], candidates=len(candidates), inner_iters=int(res.info.iter), inner_tol=eps))
            solver_history.append(dict(**history[-1], solver_primal=float(res.info.prim_res),
                solver_dual=float(res.info.dual_res), rho_updates=int(res.info.rho_updates),
                rho_estimate=float(res.info.rho_estimate), polish_status=int(res.info.status_polish),
                **{key:diag[key] for key in ('primal','stationarity','complementarity',
                    'eps_primal','eps_stationarity','eps_complementarity')}))
            if res.info.status_val == 1 and eps <= final_eps and ok:
                converged, reason = True, 'certified'
                break
            if res.info.status_val not in (1, 2, 7):
                reason = 'solver_failure'
                break
            if cp == max_cuts:
                reason = 'inner_maxiter' if res.info.status_val != 1 else 'maxcut'
                break
            if candidates:
                pairs.extend(candidates); active.update(candidates)
                solver = None
                eps = max(final_eps, min(1e-3, .01*diag['max_viol'])) if adaptive else final_eps
            elif eps > final_eps: eps = final_eps
            elif res.info.status_val == 1: eps *= .1
        diag.update(converged=int(converged), stop_reason=reason, cp_iters=cp, n_active=len(pairs),
            total_inner_iters=iterations, setups=setups, same_qp_resumes=resumes, expanded_primal_starts=primal_starts,
            solver_status=int(res.info.status_val), inner_tolerance=eps,
            matrix_payload_bytes=payload_peak, dense_normal_matrix_bytes=8*(p.n*p.k)**2,
            scan_workspace_bytes=8*min(block_size, p.n)*p.n, parameters=p.n*p.k,
            seconds=time.perf_counter()-start, pivot=p.pivot+1, N=p.n, m=p.X.shape[1], s=p.Y.shape[1],
            tolerance=tolerance, normtolerance=normtolerance, objective_ridge=0.0)
        diag.update(status_resets=status_resets, **init_summary)
        return _working_set_metadata(dict(distance=d, residual=-d, coefficients=coef, diagnostics=diag, history=history,
                    theta=theta, dual=res.y, pairs=np.asarray(pairs), solver_history=solver_history,
                    initialization=initialization, trace=log.getvalue()), working_set)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rts', choices=('vrs','crs'), default='vrs')
    ap.add_argument('--input', required=True)
    ap.add_argument('--inputs', nargs='+', required=True)
    ap.add_argument('--outputs', nargs='+', required=True)
    ap.add_argument('--gx', nargs='+', type=float, required=True)
    ap.add_argument('--gy', nargs='+', type=float, required=True)
    ap.add_argument('--output', required=True, help='output filename stem (.csv/.json/.meta/.history)')
    ap.add_argument('--tolerance', type=float, default=1e-6)
    ap.add_argument('--normtolerance', type=float, default=1e-8)
    ap.add_argument('--maxiter', type=int, default=200000)
    ap.add_argument('--maxcut', type=int, default=100)
    ap.add_argument('--pivot', type=int, default=0, help='1-based; zero chooses the largest scaled direction')
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--full-start', choices=('auto','cold','interior'), default='auto')
    ap.add_argument('--trace', action='store_true', help='include solver iteration logs in the JSON report')
    ap.add_argument('--solver', choices=('osqp','clarabel'), default='osqp')
    ap.add_argument('--ipmmaxiter', type=int, default=100, help='iterations per direct Clarabel call')
    ap.add_argument('--cutbatch', type=int, default=1, help='new Afriat pairs per observation (Clarabel CP)')
    ap.add_argument('--innercap', type=float, default=1e-3, help='intermediate tolerance ceiling (Clarabel CP)')
    ap.add_argument('--inner-policy', choices=('adaptive','monotone'), default='adaptive',
                    help='retain tightening after QP 2 with monotone (Clarabel CP)')
    ap.add_argument('--working-set', choices=('none','once','bounded'), default='none',
                    help='optional slack-row pruning with one-shot numerical recovery (Clarabel CP)')
    a = ap.parse_args()
    with open(a.input, newline='') as f: data = list(csv.DictReader(f))
    X = np.array([[float(row[v]) for v in a.inputs] for row in data])
    Y = np.array([[float(row[v]) for v in a.outputs] for row in data])
    result = fit_ddf(X,Y,a.gx,a.gy,tolerance=a.tolerance,normtolerance=a.normtolerance,
                     max_iter=a.maxiter,max_cuts=a.maxcut,full=a.full,pivot=None if a.pivot==0 else a.pivot-1,
                     full_start=a.full_start,trace=a.trace,solver=a.solver,ipm_max_iter=a.ipmmaxiter,
                     cut_batch=a.cutbatch,inner_cap=a.innercap,inner_policy=a.inner_policy,working_set=a.working_set,rts=a.rts)
    diag = result['diagnostics']
    np.savetxt(a.output+'.csv',np.column_stack((result['distance'],result['residual'],result['coefficients'])),delimiter=',',fmt='%.17g')
    report = dict(version=VERSION, rts=a.rts, solver=a.solver, working_set=a.working_set,
                  working_set_history=result['working_set_history'],recovery_history=result['recovery_history'], inner_policy=a.inner_policy, numpy=np.__version__, scipy=scipy.__version__, osqp=osqp.__version__,
                  diagnostics=diag, history=result['history'],solver_history=result['solver_history'],
                  initialization=result['initialization'],trace=result['trace'])
    Path(a.output+'.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    # Numeric-only metadata for the Stata bridge; no generated executable code.
    numeric = {k:v for k,v in diag.items() if k not in ('stop_reason','kkt_profile')}
    with open(a.output+'.meta','w',newline='') as f:
        writer=csv.writer(f);writer.writerow(numeric.keys());writer.writerow(numeric.values())
    with open(a.output+'.history','w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(result['history'][0]));writer.writeheader();writer.writerows(result['history'])
    Path(a.output+'.reason').write_text(diag['stop_reason']+'\n')


if __name__ == '__main__':
    try: main()
    except Exception as exc:
        print(type(exc).__name__+': '+str(exc),file=sys.stderr)
        sys.exit(1)
