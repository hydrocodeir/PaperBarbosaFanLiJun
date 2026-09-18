"""Figures separating index construction, threshold ties and station selection."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .publication_figures import INDEXES,NAMES
from .publication_regimes import REGIMES,LABELS

SCENARIOS=['full_w11_t7','full_w5_t7','fixed_w11_t7_raw','fixed_w11_t7_corrected','fixed_w5_t7_corrected','fixed_w5_t8_corrected']
LABELS_INDEX=['Full record · 11 d · T7','Full record · 5 d · T7','Early reference · 11 d · T7 · raw','Early reference · 11 d · T7 · corrected','Early reference · 5 d · T7 · corrected','Early reference · 5 d · T8 · corrected']
COLORS=['#555D63','#899EA7','#BF7635','#267C91','#6974A6','#965D88']
MARKERS=['o','s','^','D','v','P']

def plot_index_period_changes(out,records,save):
    d=pd.read_csv(out/'tables/index_definition_period_summary.csv')
    fig,ax=plt.subplots(figsize=(7.1,4.2),layout='constrained')
    for k,(scenario,label,color,marker) in enumerate(zip(SCENARIOS,LABELS_INDEX,COLORS,MARKERS)):
        local=d.loc[d.scenario==scenario].set_index('index_name').loc[INDEXES]
        ax.scatter(local.change_days,np.arange(4)+(k-2.5)*.10,s=28,color=color,marker=marker,label=label)
    ax.axvline(0,color='#7A8285',lw=.7)
    ax.set(yticks=range(4),yticklabels=NAMES,xlabel='Late minus early annual count (days per year)',title='Fixed 108-station network · 2008–2024 minus 1991–2007',ylim=(-.6,3.6))
    ax.invert_yaxis()
    fig.legend(*ax.get_legend_handles_labels(),loc='outside lower center',ncol=2,fontsize=6.4,frameon=False)
    save(fig,out,'fig04_fixed_baseline_changes',records)

def plot_zero_regimes(out,records,save):
    summary=pd.read_csv(out/'tables/zero_threshold_summary.csv')
    summary=summary.loc[(summary.threshold_mode=='refitted')&(summary.component=='joint_change')]
    fig,axes=plt.subplots(1,2,figsize=(7.1,4.5),sharex=True,sharey=True,layout='constrained')
    specs=[('strict','all','Strict · all stations','#287E9B','o'),('strict','positive','Strict · positive cutoffs','#555D63','D'),('dry_inclusive_only','all','Dry ties included · all','#B47C35','^')]
    for k,(definition,ax) in enumerate(zip(['annual','warm_season'],axes)):
        for j,(rule,stratum,label,color,marker) in enumerate(specs):
            local=summary.loc[(summary.definition==definition)&(summary.rule==rule)&(summary.stratum==stratum)].set_index('climate_regime').loc[REGIMES]
            y=np.arange(6)+(j-1)*.19
            ax.hlines(y,local.ci_low_pp,local.ci_high_pp,color=color,lw=.9)
            ax.scatter(local.estimate_pp,y,s=18,marker=marker,color=color,label=label)
        local=summary.loc[(summary.definition==definition)&(summary.rule=='strict')&(summary.stratum=='all')].set_index('climate_regime').loc[REGIMES]
        for y,(_,r) in enumerate(local.iterrows()):ax.text(.99,y,f'{int(r.n_stations-r.n_zero)}/{int(r.n_stations)}',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=6)
        ax.axvline(0,color='gray',lw=.6);ax.set_xlim(-5,76)
        ax.set(xticks=[0,20,40,60],yticks=range(6),yticklabels=LABELS,title=f"({'ab'[k]}) {'Annual' if k==0 else 'June–September'}")
    axes[0].invert_yaxis()
    fig.supxlabel('Joint-frequency change (percentage points); exploratory 95% intervals',fontsize=8)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='outside upper center',ncol=3,frameon=False,fontsize=6.1)
    save(fig,out,'fig09_compound_climate_regimes',records)

def plot_definition_supplement(out,records,save):
    trends=pd.read_csv(out/'tables/index_definition_trends.csv')
    fig,axes=plt.subplots(2,2,figsize=(7.1,5.6),layout='constrained',sharex=True)
    for k,(index,ax) in enumerate(zip(INDEXES,axes.flat)):
        d=trends.loc[(trends.index_name==index)&(trends.metric=='Delta1')].set_index('scenario').loc[SCENARIOS]
        for y,((_,r),color,marker) in enumerate(zip(d.iterrows(),COLORS,MARKERS)):
            ax.plot([r.ci_low,r.ci_high],[y,y],color=color,lw=1.3);ax.scatter(r.estimate,y,s=22,color=color,marker=marker)
        ax.axvline(0,color='gray',lw=.6,ls='--')
        ax.set(title=f'({chr(97+k)}) {NAMES[k]}',yticks=range(6),yticklabels=['Full 11 / T7','Full 5 / T7','Early 11 / raw','Early 11 / corrected','Early 5 / corrected','Early 5 / T8 corrected'])
        ax.invert_yaxis()
    fig.supxlabel('Upper-minus-lower slope contrast (days per decade); conditional 95% intervals',fontsize=8)
    save(fig,out,'figS14_index_definition_asymmetry',records)
    effects=pd.read_csv(out/'tables/zero_threshold_paired_effects.csv')
    effects=effects.loc[(effects.definition=='warm_season')&(effects.threshold_mode=='refitted')]
    fig,axes=plt.subplots(1,2,figsize=(7.1,4),layout='constrained',sharey=True)
    for k,(contrast,title,color,ax) in enumerate(zip(['dry_tie_change','hot_tie_increment'],['Include dry ties only','Then include hot ties'],['#B47C35','#756898'],axes)):
        d=effects.loc[effects.contrast==contrast].set_index('climate_regime').loc[['All',*REGIMES]]
        ax.hlines(range(7),d.ci_low_pp,d.ci_high_pp,color=color,lw=1.2);ax.scatter(d.estimate_pp,range(7),color=color,s=23)
        ax.axvline(0,color='gray',lw=.7,ls='--');ax.set(title=f'({chr(97+k)}) {title}',yticks=range(7),yticklabels=['Network',*LABELS],xlabel='Paired change (percentage points)')
    axes[0].invert_yaxis()
    save(fig,out,'figS15_zero_threshold_tie_effects',records)
