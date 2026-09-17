"""Independent linear-program check of the displayed warming-response slopes."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linprog

root=Path(__file__).resolve().parent
idx=['warm_days','warm_nights','cool_days','cool_nights']
a=pd.read_csv(root/'outputs/tables/fixed_baseline_annual_extreme_indices.csv').groupby('year')[idx].mean().reset_index().merge(pd.read_csv(root/'outputs/tables/regional_temperature_anomaly.csv'),on='year')
old=pd.read_csv(root/'outputs/tables/warming_link_network_quantile_response.csv');rows=[]
for name in idx:
    x=a.regional_temperature_anomaly_c.to_numpy();y=a[name].to_numpy();X=np.column_stack([np.ones(len(x)),x]);n=len(y)
    for t in [.1,.5,.9]:
        fit=linprog(np.r_[0.,0.,np.full(n,t),np.full(n,1-t)],A_eq=np.column_stack([X,np.eye(n),-np.eye(n)]),b_eq=y,bounds=[(None,None)]*2+[(0,None)]*(2*n),method='highs')
        v=old.loc[(old.index_name==name)&(old.tau==t),'slope_per_c'].iloc[0]
        rows.append({'index_name':name,'tau':t,'saved_slope':v,'linear_program_slope':fit.x[1],'absolute_difference':abs(v-fit.x[1]),'solver_success':bool(fit.success)})
result=pd.DataFrame(rows)
result.to_csv(root/'outputs/audit_cleanup/network_warming_solver_check.csv',index=False)
print(result.to_string(index=False))
assert result.solver_success.all()
assert result.absolute_difference.max()<.01, 'Review numerical convergence before using the figure.'
