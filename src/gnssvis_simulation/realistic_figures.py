"""Render simulations with the existing GNSSVis layout and object table."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from gnssvis.plot import Visualizer, COLORS_DISTINCT
from .beam import beam_offsets_deg, beam_gain


def plot_realistic_figures(r, output: Path):
    sc, g = r.scenario, r.geometry
    output = Path(output)
    t = mdates.date2num(sc.timestamps_utc)
    end = mdates.date2num(sc.end_utc)
    f = sc.frequencies_hz / 1e6
    beam = sc.section('beam')
    p = sc.section('pointing')
    x, y = beam_offsets_deg(g.alt_deg, g.az_deg, p['altitude_deg'], p['azimuth_deg'])
    rows = np.flatnonzero(r.half_power_mask.any(axis=1))
    colors = {row: COLORS_DISTINCT[i] for i, row in enumerate(rows)}
    sun_mask = (beam_gain(r.solar['offset_x_deg'], r.solar['offset_y_deg'],
                         beam['reference_frequency_hz'], beam) >= .5) & (r.solar['alt_deg'] > 0)
    p0 = float(sc.section('noise')['mean_power_relative'])
    limits = (float((10*np.log10(r.total_power/p0)).min()),
              float((10*np.log10(r.total_power/p0)).max()))

    def draw(data, title, name, crop=None, doppler=False, solar=True, gps=True):
        with plt.rc_context(matplotlib.rcParamsDefault):
            fig = plt.figure(figsize=(16,10))
            gs = fig.add_gridspec(2,2,width_ratios=[3,1.2],height_ratios=[1,3],
                                 hspace=.35,wspace=.35)
            top=fig.add_subplot(gs[0,0]); wf=fig.add_subplot(gs[1,0],sharex=top)
            table=fig.add_subplot(gs[:,1]); table.axis('off')
            legend=[]
            for row in (rows if gps else []):
                mask=r.half_power_mask[row]
                values=r.doppler_hz[row]/1000 if doppler else y[row]
                top.scatter(t[mask],values[mask],s=6.25,color=colors[row],alpha=.85,linewidths=0)
                legend.append((g.satellites[row].name,colors[row]))
            if solar and not doppler:
                top.scatter(t[sun_mask],r.solar['offset_y_deg'][sun_mask],s=35,
                            marker='*',color='goldenrod',linewidths=0)
                legend.append(('Sun','goldenrod'))
            top.axhline(0,color='gray',ls='--',lw=.5,alpha=.5)
            top.set_ylabel('Doppler L1 [kHz]' if doppler else 'Relative elevation offset [°]')
            top.set_title('Satellite Doppler inside the half-power beam' if doppler else
                          'Object trajectories relative to the fixed pointing',fontsize=11,fontweight='bold')
            if not doppler: top.set_ylim(-25,25)
            top.grid(True,alpha=.3)
            top.tick_params(labelbottom=False)
            keep=np.ones(len(f),dtype=bool) if crop is None else ((f>=crop[0])&(f<=crop[1]))
            vals=10*np.log10(np.maximum(data[keep],1e-20)/p0)
            half=sc.channel_width_hz/2e6
            im=wf.imshow(vals,origin='lower',aspect='auto',interpolation='none',cmap='plasma',
                         extent=[t[0],end,f[keep][0]-half,f[keep][-1]+half],
                         vmin=limits[0],vmax=limits[1])
            for hz,label in ([(1176.45,'GPS L5'),(1227.60,'GPS L2'),(1575.42,'GPS L1')] if gps else []):
                if f[keep][0]-half<=hz<=f[keep][-1]+half:
                    wf.axhline(hz,color='cyan',ls='--',alpha=.7,lw=1)
                    wf.text(1.015,hz,f'{label}\n{hz:.2f} MHz',transform=wf.get_yaxis_transform(),
                            fontsize=8,va='center')
            wf.set_xlim(t[0],end)
            wf.xaxis.set_major_locator(mdates.HourLocator(interval=4))
            wf.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m\n%H:%M'))
            wf.set_xlabel('Time [UTC]');wf.set_ylabel('Frequency [MHz]')
            wf.set_title(title,fontsize=11,fontweight='bold')
            fig.colorbar(im,ax=wf,orientation='horizontal',label='Relative power [dB re mean background]',pad=.16)
            Visualizer.draw_object_table(table,legend,fontsize=9,row_scale=2.1,name_width=24)
            table.set_title('Objects crossing half power at L1',fontsize=11,fontweight='bold')
            fig.suptitle('GNSSVis — synthetic 24-hour observation',fontsize=13)
            fig.text(.125,.025,'Fixed alt/az pointing; elliptical half-power selection at 1575.42 MHz. '
                     'Gaussian wings contribute outside this contour. Illustrative powers.',fontsize=8)
            fig.savefig(output/name,dpi=150)
            plt.close(fig)

    draw(r.total_power,'Waterfall — GPS + Sun + noise','01_waterfall_24h.png')
    draw(r.noise_power+r.satellite_power,'Waterfall — GPS + noise (Sun excluded)',
         '02_componentes.png',solar=False)
    draw(r.noise_power+r.satellite_power,'Waterfall — GPS L1 detail','03_l1_doppler.png',
         crop=(1558,1590),doppler=True,solar=False)
    draw(r.noise_power+r.solar_power,'Waterfall — Sun + noise (GPS excluded)','05_comparacao.png', gps=False)
    from .realistic import bpsk_channel_psd
    with plt.rc_context(matplotlib.rcParamsDefault):
        fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained')
        freq=np.linspace(-25e6,25e6,2001)
        for rate in [1.023e6,10.23e6]:
            psd=bpsk_channel_psd(freq,50000,0,rate)
            axs[0].plot(freq/1e6,10*np.log10(np.maximum(psd,1e-6)),label=f'{rate/1e6:g} Mchip/s')
        axs[0].set(xlabel='Offset from carrier [MHz]',ylabel='Normalized PSD [dB]',ylim=(-45,1))
        ang=np.linspace(-90,90,500)
        for hz,label in [(1176.45e6,'L5'),(1575.42e6,'L1')]:
            axs[1].plot(ang,beam_gain(ang,0,hz,beam),label=label)
        axs[1].set(xlabel='Horizontal offset [°]',ylabel='Power gain')
        for ax in axs: ax.legend();ax.grid(alpha=.3)
        fig.savefig(output/'04_modelos.png',dpi=150);plt.close(fig)
