"""Stata 16 bridge to unchanged JSS pyStoNED 0.7.3 QLE; no model reimplementation.

Invoked only by implementation(pystoned) method(qle). Source hashes and the
NumPy/SciPy versions are checked before use. No network calls. Writes only the call-specific result file supplied by Stata.
"""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import numpy as np
import pystoned
from pystoned import StoNED
from pystoned.constant import CET_ADDI, FUN_COST, FUN_PROD, RED_QLE


def verify_reference():
    spec=json.loads(Path(__file__).with_name('stoned_py073_reference.json').read_text())
    for package,version in spec['runtime'].items():
        if importlib.metadata.version(package)!=version:
            raise RuntimeError(f'QLE reference requires {package}=={version}')
    root=Path(pystoned.__file__).parent
    for name,expected in spec['package_files'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=expected:
            raise RuntimeError('Modified reference source: '+name)


class ResidualModel:
    def __init__(self,eps,cost):
        self.eps=eps
        self.x=np.ones((len(eps),1))
        self.y=np.ones(len(eps))
        self.fun=FUN_COST if cost else FUN_PROD
        self.cet=CET_ADDI
        # This adapter is not a new first-stage convergence certificate.
        self.optimization_status=1
    def get_residual(self): return self.eps.copy()


def main(args):
    input_path,output_path,cost=args
    verify_reference()
    eps=np.genfromtxt(input_path,delimiter=',',skip_header=1,ndmin=1)
    obj=StoNED.StoNED(ResidualModel(eps,int(cost)))
    original=StoNED.opt.minimize
    records=[]
    def observe(*a,**kw):
        result=original(*a,**kw);records.append(result);return result
    try:
        StoNED.opt.minimize=observe
        obj.get_unconditional_expected_inefficiency(RED_QLE)
    finally:
        StoNED.opt.minimize=original
    if len(records)!=1: raise RuntimeError('Expected one public QLE optimizer call')
    values=[obj.sigma_u,obj.sigma_v,int(records[0].success),records[0].x[0]]
    verify_reference()
    # Exclusive creation prevents stale output from being mistaken for success.
    with open(output_path,'x') as f:
        f.write('STONED_QLE_RAW_V1\n')
        f.write('\n'.join(format(float(v),'.17g') for v in values)+'\n')
        message=str(records[0].message).replace('\n',' ').replace('"','').replace("'",'')
        f.write('\n'.join((str(records[0].status), str(getattr(records[0],'nfev','')),
                            message, str(obj.sigma_u), str(obj.sigma_v),
                            str(records[0].x[0])))+'\n')
        f.write('END_STONED_QLE_RAW_V1\n')

if __name__=='__main__':main(sys.argv[1:])
