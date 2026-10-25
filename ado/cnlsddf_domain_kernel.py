"""CRS geometry diagnostics only. Never fit, clip quantities, or replace signed distances.

A caller must independently certify the estimator and bind the data/planes. This
pure function checks geometric inputs, not optimizer termination or provenance.
Physical endpoints describe only the supplied polyhedral technology and ray.
"""
import numpy as np

def diagnose(X,Y,gx,gy,coefficients,distance,scales,quantity_tolerance=1e-10,
             normalization_tolerance=1e-8,envelope_tolerance=1e-7,monotonicity_tolerance=1e-7):
    X,Y,gx,gy,c,d,s=(np.asarray(v,dtype=float) for v in (X,Y,gx,gy,coefficients,distance,scales))
    def need(ok,message):
        if not ok:raise ValueError(message)
    need(X.ndim==Y.ndim==2 and len(X)==len(Y)>0 and X.shape[1]>0 and Y.shape[1]>0,'quantity dimensions')
    n,m=X.shape;k=m+Y.shape[1]
    need(gx.shape==(m,) and gy.shape==(Y.shape[1],) and d.shape==(n,) and s.shape==(k,),'direction/distance/scale dimensions')
    need(c.ndim==2 and len(c)>0 and c.shape[1]==k+1,'coefficient dimensions')
    need(all(np.all(np.isfinite(v)) for v in (X,Y,gx,gy,c,d,s)),'nonfinite input')
    need(np.all(X>=0) and np.all(Y>=0) and np.all(s>0),'nonnegative quantities and positive scales required')
    g=np.r_[gx,gy];need(np.all(g>=0) and np.any(g>0),'nonnegative nonzero direction required')
    need(all(np.isfinite(t) and t>0 for t in (quantity_tolerance,normalization_tolerance,envelope_tolerance,monotonicity_tolerance)),'positive finite tolerances required')
    need(np.all(c[:,0]==0),'CRS zero intercept required')
    w=c[:,1:];q=w@g;speed=float(np.max(g/s))
    need(np.isfinite(speed) and speed>0 and np.all(q>0) and np.max(np.abs(q-1))<=normalization_tolerance,'direction normalization')
    need(np.min(w*s)>=-monotonicity_tolerance,'monotonicity')
    V=np.c_[X,-Y];F=np.empty(n);root=np.empty(n)
    for start in range(0,n,128):
        values=V[start:start+128]@w.T
        F[start:start+128]=np.min(values,axis=1)
        # q is numerically close to one, not assumed to be exactly one.
        root[start:start+128]=np.min(values/q,axis=1)
    need(np.all(np.isfinite(F)) and np.all(np.isfinite(root)),'nonfinite envelope/root')
    need(np.max(np.abs(d-F))<=envelope_tolerance,'distance does not match supplied envelope')
    lower=np.max(-Y[:,gy>0]/gy[gy>0],axis=1) if np.any(gy>0) else np.full(n,-np.inf)
    upper=np.min(X[:,gx>0]/gx[gx>0],axis=1) if np.any(gx>0) else np.full(n,np.inf)
    hi=np.minimum(root,upper);feasible=hi>=lower
    raw=np.c_[X-d[:,None]*gx,Y+d[:,None]*gy];scaled=raw/s
    need(np.all(np.isfinite(raw)),'nonfinite raw target')
    low_margin=(root-lower)*speed;cap_margin=(root-upper)*speed
    near=(np.abs(low_margin)<=quantity_tolerance)|(np.abs(cap_margin)<=quantity_tolerance)
    status=np.where(~feasible,'empty',np.where(root>upper,'input_cap','zero_surface')).astype('<U20');status[near]='near_boundary'
    endpoint_kind=np.where(~feasible,'none',np.where(root>upper,'input_cap',np.where(root==upper,'zero_and_input_cap','zero_surface')))
    physical_step=np.where(feasible,hi,np.nan)
    physical=np.c_[X-physical_step[:,None]*gx,Y+physical_step[:,None]*gy]
    physical[~feasible]=np.nan
    plane_margin=np.full(n,np.nan);min_quantity=np.full(n,np.nan)
    for start in range(0,n,128):
        stop=min(n,start+128);ix=np.arange(start,stop)[feasible[start:stop]]
        if not len(ix):continue
        pv=np.c_[physical[ix,:m],-physical[ix,m:]]
        plane_margin[ix]=np.min((pv@w.T)/q,axis=1)*speed
        min_quantity[ix]=np.min(physical[ix]/s,axis=1)
    return dict(distance=d.copy(),envelope=F,ray_zero_step=root,lower=lower,upper=upper,
                raw_target=raw,raw_negative_any=np.any(raw<0,axis=1),
                raw_negative_material=np.any(scaled < -quantity_tolerance,axis=1),
                raw_negative_input=np.any(scaled[:,:m] < -quantity_tolerance,axis=1),
                raw_negative_output=np.any(scaled[:,m:] < -quantity_tolerance,axis=1),
                minimum_scaled_raw=np.min(scaled,axis=1),exact_ray_nonempty=feasible,
                nonnegative_step_available=feasible & (hi>=0),observed_in_technology_exact=F>=0,status=status,
                endpoint_kind=endpoint_kind,physical_step=physical_step,physical_target=physical,
                lower_margin_scaled=low_margin,input_cap_margin_scaled=cap_margin,
                physical_min_scaled_quantity=min_quantity,physical_plane_margin_scaled=plane_margin,
                direction_speed=np.array(speed),normalization_gap=np.array(np.max(np.abs(q-1))))
