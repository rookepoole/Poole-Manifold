#!/usr/bin/env python3
"""Independent Python reference checks for the Poole Manifold PM0 base rule.

Created as AI-assisted verification code for the author's 3D B5-7/S5-9
research. No outside or hardware data are used. Reproducible seed 20261018.
Python 3 + numpy + matplotlib (plots optional).
"""
from __future__ import annotations
import csv
import json
import platform
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
ASSETS = ROOT / 'assets'
SEED = 20261018
O3 = tuple((i,j,k) for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1) if (i,j,k)!=(0,0,0))
O2 = tuple((i,j) for i in (-1,0,1) for j in (-1,0,1) if (i,j)!=(0,0))


def step3d(x: np.ndarray, ablated: bool=False) -> np.ndarray:
    """26 offsets; periodic in all 3 dimensions; synchronous, exact integer counts."""
    n = np.zeros(x.shape, dtype=np.int16)
    for offs in O3:
        n += np.roll(x, offs, axis=(0,1,2))
    return ((n >= 5) & (n <= (7 if ablated else 7 + 2*x))).astype(np.uint8)


def step3d_direct(x: np.ndarray) -> np.ndarray:
    """Slow, direct indexing baseline for small arrays; no np.roll calls."""
    size = x.shape[0]
    y = np.zeros_like(x)
    for i in range(size):
        for j in range(size):
            for k in range(size):
                n = sum(int(x[(i+di)%size,(j+dj)%size,(k+dk)%size]) for di,dj,dk in O3)
                y[i,j,k] = int(5 <= n <= 7 + 2*int(x[i,j,k]))
    return y


def mask_rectangle(L: int, a: int, b: int, start=(5,5)) -> np.ndarray:
    d = np.zeros((L,L),dtype=np.uint8)
    i,j=start
    if not (2<=a<=L-i-2 and 2<=b<=L-j-2):
        raise ValueError('Choose a,b >=2 with at least a two-site outer margin')
    d[i:i+a, j:j+b]=1
    return d


def sheet(d: np.ndarray) -> np.ndarray:
    L=d.shape[0]
    assert d.shape==(L,L) and L>=7
    x=np.zeros((L,L,L),dtype=np.uint8)
    x[L//2]=(1-d)
    return x


def predicted_from_mask(d:np.ndarray) -> np.ndarray:
    """The q_D/r_D planar first-step selector, not the 3D update."""
    L=d.shape[0]
    q=np.zeros((L,L),dtype=np.int16)
    for o in O2:
        q+=np.roll(d,o,axis=(0,1))
    r=q+d
    plane=((d==0)&(q<=3))|((d==1)&(q>=1)&(q<=3))
    halo=(r>=2)&(r<=4)
    y=np.zeros((L,L,L),dtype=np.uint8)
    z=L//2
    y[z]=plane
    y[z-1]=halo
    y[z+1]=halo
    return y


def event_counts_2d(d: np.ndarray) -> dict:
    """Compute planar selector births + channels, independent of 3D summation."""
    q=np.zeros_like(d,dtype=np.int16)
    for o in O2:
        q+=np.roll(d,o,axis=(0,1))
    r=q+d
    in_birth=(d==1)&(q>=1)&(q<=3)
    halo_birth=(r>=2)&(r<=4)
    deaths=(d==0)&(q>=4)
    kin=(8-q)[in_birth]
    knorm=(9-r)[halo_birth]
    k=np.concatenate((kin, knorm, knorm))
    return dict(births=int(len(k)),deaths=int(deaths.sum()),
                channels=[int((k==n).sum()) for n in (5,6,7)],
                inplane=int(in_birth.sum()),normal_each=int(halo_birth.sum()))


def render_figures(rects:list[dict]):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    plt.rcParams.update({'font.size':11,'axes.titlesize':12,'axes.labelsize':11,
                         'figure.facecolor':'white', 'savefig.facecolor':'white',
                         'axes.spines.top':False, 'axes.spines.right':False,
                         'font.family':'DejaVu Sans'})

    # Figure 01: 3D cropped voxel input, accepted births, final result.
    d=mask_rectangle(25,3,4)
    x=sheet(d)
    y=step3d(x)
    b=y & (1-x)
    z=12
    a1=x[z-1:z+2, 2:14, 3:16].transpose(1,2,0).astype(bool)
    a2=b[z-1:z+2, 2:14, 3:16].transpose(1,2,0).astype(bool)
    a3=y[z-1:z+2, 2:14, 3:16].transpose(1,2,0).astype(bool)
    fig=plt.figure(figsize=(14,4.5),layout='constrained')
    labels=['A  /  Initial punctured sheet','B  /  New births (40)','C  /  First updated state']
    for i,arr in enumerate((a1,a2,a3)):
        ax=fig.add_subplot(1,3,i+1,projection='3d')
        colors=np.empty(arr.shape,dtype=object)
        colors[:]='#617A8C' if i !=1 else '#EE9B4A'
        if i==2:
            born=a2
            colors[born]='#EE9B4A'
        ax.voxels(arr,facecolors=colors,edgecolor='#F8F9FB',linewidth=.20,shade=True)
        ax.set(xlim=(0,12),ylim=(0,13),zlim=(0,3))
        ax.set_box_aspect((12,13,3),zoom=.86)
        ax.view_init(elev=27,azim=-64)
        ax.set_axis_off()
        ax.set_title(labels[i],fontweight='bold',pad=8)
    fig.suptitle('PM0 B5–7/S5–9  •  3D one-step response to a 3 × 4 sheet puncture',fontweight='bold',fontsize=15)
    fig.savefig(ASSETS/'fig01_3d_rectangular_defect.png',dpi=190,bbox_inches='tight');plt.close(fig)

    # Figure 02: planar mask and the two normal birth layers.
    L=25;z=12
    d=mask_rectangle(L,3,4)
    x=sheet(d);y=step3d(x)
    plane_initial=x[z]; new_inplane=((1-x)*y)[z]
    normal_birth=((1-x)*y)[z+1]
    fig,axes=plt.subplots(1,3,figsize=(13.6,4.7),layout='constrained')
    cmap1=ListedColormap(['#F5F7FA','#283F52'])
    cmap2=ListedColormap(['#F5F7FA','#EE9B4A'])
    titles=['A / Original sheet (puncture)', 'B / Four corner rebirths', 'C / Normal-layer birth halo']
    for ax,mat,cm,title in zip(axes,[plane_initial,new_inplane,normal_birth],[cmap1,cmap2,cmap2],titles):
        ax.imshow(mat[1:20,1:20],interpolation='nearest',vmin=0,vmax=1,cmap=cm,origin='lower')
        ax.set_title(title,fontweight='bold');ax.set_xticks([]);ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(True);s.set_edgecolor('#D6DDE4')
    fig.suptitle('First response: 4 in-plane births + 18 births in each adjacent layer = 40',fontsize=14,fontweight='bold')
    fig.savefig(ASSETS/'fig02_birth_geometry.png',dpi=200,bbox_inches='tight');plt.close(fig)

    # Figure 03: exact birth count vs perimeter, across all rectangles.
    perims=np.array([2*(r['a']+r['b']) for r in rects]);births=np.array([r['births'] for r in rects]);
    xgrid=np.arange(perims.min(),perims.max()+1)
    fig,ax=plt.subplots(figsize=(9.5,5.2),layout='constrained')
    ax.plot(xgrid,2*xgrid+12,label=r'Exact law $B=2P+12$',lw=2.4,c='#294E70',zorder=1)
    ax.scatter(perims,births,s=14,alpha=.42,c='#EE9B4A',label='841 planar q/r selector counts',zorder=2)
    ax.set(xlabel='Rectangle lattice perimeter  P = 2(a + b)',ylabel='First-step birth count',title='Rectangular punctures: exact perimeter-to-event law')
    ax.grid(True,alpha=.15);ax.legend(frameon=False)
    fig.savefig(ASSETS/'fig03_perimeter_law.png',dpi=200,bbox_inches='tight');plt.close(fig)

    # Figure 04: channel counts for square punctures.
    hs=np.arange(2,31)
    fig,ax=plt.subplots(figsize=(9.5,5.2),layout='constrained')
    ax.plot(hs,np.full_like(hs,12),label='Birth channel N=5 (12)',lw=2.2)
    ax.plot(hs,8*hs-16,label='Birth channel N=6 (8h−16)',lw=2.2)
    ax.plot(hs,np.full_like(hs,16),label='Birth channel N=7 (16)',lw=2.2)
    ax.plot(hs,8*hs+12,label='Total (8h+12)',lw=2.8,ls='--')
    ax.set(xlabel='Square hole side length  h',ylabel='First-step event count',title='Event-channel decomposition of square punctures')
    ax.grid(True,alpha=.15);ax.legend(frameon=False)
    fig.savefig(ASSETS/'fig04_channel_spectrum.png',dpi=200,bbox_inches='tight');plt.close(fig)


def main():
    DATA.mkdir(parents=True,exist_ok=True)
    ASSETS.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(SEED)
    assert len(O3)==26 and len(O2)==8

    # Direct no-roll baseline for small random entire 3D states.
    direct_cases=5
    for _ in range(direct_cases):
        s=rng.integers(0,2,size=(7,7,7),dtype=np.uint8)
        assert np.array_equal(step3d(s),step3d_direct(s))

    blank=np.zeros((17,17),dtype=np.uint8)
    whole=sheet(blank)
    assert np.array_equal(step3d(whole),whole)
    assert step3d(whole,ablated=True).sum()==0
    single=blank.copy();single[8,8]=1
    assert np.array_equal(step3d(sheet(single)),sheet(single))

    # 841 2D algebraic-selector audits (all a,b = 2,...,30).
    rectangle_rows=[]
    for a in range(2,31):
        for b in range(2,31):
            d=mask_rectangle(40,a,b)
            obs=event_counts_2d(d)
            expected=(4*a+4*b+12,0,(12,4*a+4*b-16,16),4,2*(a+b)+4)
            actual=(obs['births'],obs['deaths'],tuple(obs['channels']),obs['inplane'],obs['normal_each'])
            assert actual==expected,(a,b,actual,expected)
            rectangle_rows.append(dict(a=a,b=b,area=a*b,perimeter=2*(a+b),births=obs['births'],
                    predicted_births=expected[0],deaths=obs['deaths'],N5=obs['channels'][0],
                    N6=obs['channels'][1],N7=obs['channels'][2],inplane_births=obs['inplane'],
                    births_per_normal_layer=obs['normal_each'],match=True))
    with (DATA/'rectangle_841_cases.csv').open('w',newline='',encoding='utf-8') as out:
        w=csv.DictWriter(out,fieldnames=list(rectangle_rows[0]));w.writeheader();w.writerows(rectangle_rows)

    # 720 random masks checked using full independent 3D cyclic neighbor sum.
    mask_rows=[]
    densities=[.03,.10,.20,.35,.50,.70]
    for idx in range(720):
        density=densities[idx % len(densities)]
        d=(rng.random((17,17))<density).astype(np.uint8)
        x=sheet(d)
        y=step3d(x)
        yp=predicted_from_mask(d)
        difference=int(np.count_nonzero(y!=yp))
        assert difference==0,(idx,difference)
        masked=event_counts_2d(d)
        counts3=dict(births=int(np.sum((1-x)*y)),deaths=int(np.sum(x*(1-y))))
        assert counts3['births']==masked['births'] and counts3['deaths']==masked['deaths']
        mask_rows.append({'case':idx+1,'density_target':density,'actual_removed_cells':int(d.sum()),
                          'predicted_births':masked['births'],'actual_births':counts3['births'],
                          'predicted_deaths':masked['deaths'],'actual_deaths':counts3['deaths'],
                          'mismatched_voxels':difference})
    with (DATA/'random_mask_720_cases.csv').open('w',newline='',encoding='utf-8') as out:
        w=csv.DictWriter(out,fieldnames=list(mask_rows[0]));w.writeheader();w.writerows(mask_rows)

    # Independently compare selected 3D rectangles to 2D selector.
    samples=[]
    for a,b in [(2,2),(3,4),(10,10),(12,12)]:
        d=mask_rectangle(25,a,b)
        x=sheet(d)
        y=step3d(x)
        yp=predicted_from_mask(d)
        counts=event_counts_2d(d)
        assert np.array_equal(y,yp)
        assert counts['births']==int(((1-x)*y).sum())
        samples.append({'a':a,'b':b,'births':counts['births'], 'deaths':counts['deaths'],
                        'channels':counts['channels'],'exact_voxel_match':True})

    report={
      'program':'PM0 deterministic 3D Moore-26 B5-7/S5-9',
      'date':'2026-10-18',
      'rng_seed':SEED,
      'python':platform.python_version(),
      'numpy':np.__version__,
      'verification_method':'Python/NumPy cyclic 26-offset full 3D reference; direct-loop 3D cross-check; independent 2D q/r defect selector',
      'direct_3D_random_cases':direct_cases,
      '3D_random_mask_matches':len(mask_rows),
      '3D_random_mask_voxel_mismatches':sum(row['mismatched_voxels'] for row in mask_rows),
      '2D_rectangle_formula_matches':len(rectangle_rows),
      'flat_sheet_fixed':True,
      'single_pore_fixed':True,
      'ablated_sheet_total_active_after_one_step':0,
      'selected_full_3D_rectangle_cases':samples,
      'not_independent_external_replication':True,
      'excluded_claims':'No general physical, quantum hardware, non-Markovian, universality, or topological protection claim is tested.'
    }
    (DATA/'verification_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    render_figures(rectangle_rows)
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
