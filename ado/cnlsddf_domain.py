"""CSV bridge for the separately certified CRS geometry kernel; no estimation."""
import argparse
from pathlib import Path
import numpy as np
from cnlsddf_domain_kernel import diagnose

FIELDS = ['rawdist','envelope','zero','lower','upper','lbinf','ubinf','step',
          'feasible','improve','observed','status','endpoint','negany','negmat',
          'negx','negy','minraw','marglo','marghi','minqty','minplane']
STATUS = dict(empty=0,zero_surface=1,input_cap=2,near_boundary=3)
ENDPOINT = dict(none=0,zero_surface=1,input_cap=2,zero_and_input_cap=3)
def convert(z):
    keys=['distance','envelope','ray_zero_step','lower','upper',None,None,'physical_step',
          'exact_ray_nonempty','nonnegative_step_available','observed_in_technology_exact',None,None,
          'raw_negative_any','raw_negative_material','raw_negative_input','raw_negative_output',
          'minimum_scaled_raw','lower_margin_scaled','input_cap_margin_scaled',
          'physical_min_scaled_quantity','physical_plane_margin_scaled']
    n=len(z['distance']);a=np.empty((n,len(keys)))
    for j,k in enumerate(keys):
        if k:a[:,j]=z[k]
    a[:,5]=np.isneginf(z['lower']);a[:,6]=np.isposinf(z['upper'])
    a[:,11]=[STATUS[s] for s in z['status']];a[:,12]=[ENDPOINT[s] for s in z['endpoint_kind']]
    return np.c_[a,z['raw_target'],z['physical_target']]
def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--m',type=int,required=True);p.add_argument('--s',type=int,required=True)
    p.add_argument('--gx',nargs='+',type=float,required=True);p.add_argument('--gy',nargs='+',type=float,required=True)
    p.add_argument('--rss',type=float,required=True);p.add_argument('--mean',type=float,required=True)
    args=p.parse_args();m,s=args.m,args.s;k=m+s
    a=np.loadtxt(args.input,delimiter=',',skiprows=1,ndmin=2)
    if m<1 or s<1 or a.shape[1]!=2*k+3 or not np.all(np.isfinite(a)):raise ValueError('Invalid source payload')
    X,Y=a[:,:m],a[:,m:k];d,r=a[:,k],a[:,k+1];c=a[:,k+2:]
    if not np.array_equal(d,-r):raise ValueError('Residual sign mismatch')
    rss=float(d@d);mean=float(d.mean())
    if not np.isfinite(args.rss) or not np.isfinite(args.mean) or abs(rss-args.rss)>1e-12*max(1,rss) or abs(mean-args.mean)>1e-12*max(1,np.max(np.abs(d))):raise ValueError('Source scalar mismatch')
    scales=np.max(np.abs(a[:,:k]),axis=0);scales[scales==0]=1
    z=diagnose(X,Y,args.gx,args.gy,c,d,scales)
    good=z['exact_ray_nonempty']
    if np.min(z['physical_min_scaled_quantity'][good],initial=0)<-1e-10 or np.min(z['physical_plane_margin_scaled'][good],initial=0)<-1e-10:raise ValueError('Endpoint feasibility check failed')
    out=convert(z);names=FIELDS+[f'raw_x{i+1}' for i in range(m)]+[f'raw_y{i+1}' for i in range(s)]+[f'target_x{i+1}' for i in range(m)]+[f'target_y{i+1}' for i in range(s)]
    with Path(args.output+'.csv').open('x') as f:
        f.write(','.join(names)+'\n')
        for row in out:f.write(','.join(format(float(x),'.17g') if np.isfinite(x) else '.' for x in row)+'\n')
    # Completion marker is written last. The Stata importer validates header, shape and missingness.
    Path(args.output+'.ok').write_text('cnlsddf_domain_v1\n')
if __name__=='__main__':main()
