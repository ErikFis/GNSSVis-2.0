# src/gnssvis/plot.py
"""Visualization for gnssvis.

Generate radio astronomy plots, including waterfalls, integrated spectra,
and waterfalls with a table of objects in the field of view.

Author: Erik H L Guimarães"""

import numpy as np
import matplotlib.pyplot as plt


class Visualizer:
    """Visualize radio astronomy data."""
    
    def plot_waterfall_with_objects(
        self,
        data: np.ndarray,
        times,
        frequencies: np.ndarray,
        sky_objects: dict,
        telescope_pointing=None,
        beamwidth_deg: float = 20.0
    ) -> None:
        """Plot a waterfall and Skyfield objects in a synchronized layout.
        
        Place the satellite table beside the plots without overlapping them."""
        import matplotlib.dates as mdates
        
        # Convert to datetime
        if hasattr(times, 'datetime'):
            time_axis = np.array(times.datetime)
        else:
            time_axis = np.array([t.to_datetime() for t in times])

        # Layout: two rows and columns, with more room for the table
        # [Trajectories | Legend]
        # [Waterfall    | Legend]
        fig = plt.figure(figsize=(16, 10))
        gs = fig.add_gridspec(2, 2, width_ratios=[3, 1.2], height_ratios=[1, 3], 
                             hspace=0.35, wspace=0.35)
        
        ax_traj = fig.add_subplot(gs[0, 0])  # Trajectories (top left)
        ax_wf = fig.add_subplot(gs[1, 0])     # Waterfall (bottom left)
        ax_legend = fig.add_subplot(gs[:, 1]) # Legend (right, spanning both rows)

        # Plot satellite trajectories
        if sky_objects:
            colors_distinct = [
                '#1f77b4', '#ff7f0e', '#d62728', '#9467bd', '#8c564b', '#e377c2', 
                '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#ff9896',
                '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
                '#393b79', '#637939', '#8c6d31', '#843c39', '#7b4173', '#5254a3',
                '#6b6ecf', '#9c9ede', '#cedb9c', '#e7cb94', '#e7ba52', '#bd9e39',
                '#ad494a', '#d6616b', '#e7969c', '#8ca252', '#b5cf6b', '#cedb9c',
                '#31a354', '#756bb1', '#636363', '#969696', '#bdbdbd', '#d9d9d9'
            ]
            telescope_dec = telescope_pointing.dec.deg if telescope_pointing else 0.0
            telescope_ra = telescope_pointing.ra.deg if telescope_pointing else 0.0
            
            # Legend data
            legend_data = []
            
            for i, (name, coord) in enumerate(sky_objects.items()):
                color = colors_distinct[i % len(colors_distinct)]
                name_lower = name.lower()
                if any(x in name_lower for x in ['gps', 'glonass', 'galileo', 'beidou', 'navstar', 'cosmos', 'satellite']):
                    marker = 'o'
                    markersize = 2.5
                else:
                    marker = '*'
                    markersize = 7
                    
                if hasattr(coord, "ra") and np.ndim(coord.ra.deg) > 0 and len(coord.ra.deg) > 1:
                    dec_relative = coord.dec.deg - telescope_dec
                    ra_relative = coord.ra.deg - telescope_ra
                    
                    ra_relative = np.where(ra_relative > 180, ra_relative - 360, ra_relative)
                    ra_relative = np.where(ra_relative < -180, ra_relative + 360, ra_relative)
                    
                    cos_dec = np.cos(np.radians(telescope_dec))
                    separation = np.sqrt(dec_relative**2 + (ra_relative * cos_dec)**2)
                    
                    fov_threshold = beamwidth_deg / 2.0
                    in_fov = separation <= fov_threshold
                    
                    # Plot all points and highlight those within the FOV
                    if len(dec_relative) > 0:
                        # Subsample the trajectory using a step derived from 128 points
                        if len(dec_relative) > 128:
                            step = len(dec_relative) // 128
                            plot_indices = np.arange(0, len(dec_relative), step)
                        else:
                            plot_indices = np.arange(len(dec_relative))
                        
                        # Separate points inside and outside the FOV
                        in_fov_mask = in_fov[plot_indices]
                        
                        # Plot points outside the FOV with low opacity
                        out_fov_indices = plot_indices[~in_fov_mask]
                        if len(out_fov_indices) > 0:
                            ax_traj.scatter(
                                time_axis[out_fov_indices],
                                dec_relative[out_fov_indices],
                                color=color,
                                marker=marker,
                                s=markersize**2 * 0.3,  # Smaller
                                alpha=0,  # Very transparent
                                linewidths=0,
                                edgecolors='none'
                            )
                        
                        # Plot points inside the FOV with high opacity
                        in_fov_indices = plot_indices[in_fov_mask]
                        if len(in_fov_indices) > 0:
                            ax_traj.scatter(
                                time_axis[in_fov_indices],
                                dec_relative[in_fov_indices],
                                color=color,
                                marker=marker,
                                s=markersize**2,
                                alpha=0.85,
                                linewidths=0,
                                edgecolors='none'
                            )
                        
                        # Add to the legend
                        legend_data.append((name, color))
                else:
                    dec_relative = coord.dec.deg - telescope_dec if hasattr(coord, "dec") else np.nan
                    ax_traj.hlines(
                        y=dec_relative,
                        xmin=time_axis[0],
                        xmax=time_axis[-1],
                        color=color,
                        linestyles="--",
                        alpha=0.4
                    )
                    legend_data.append((name, color))
            
            ax_traj.set_ylabel("Relative Declination [°]", fontsize=10)
            ax_traj.set_title("Satellite trajectory (relative to the pointing)", fontsize=11, fontweight='bold')
            ax_traj.set_ylim(-25, 25)
            ax_traj.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
            ax_traj.grid(True, alpha=0.3)
            ax_traj.xaxis_date()
            ax_traj.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        else:
            ax_traj.set_title("No objects visible in the FoV")
            legend_data = []

        # Bottom panel: waterfall
        vmin, vmax = np.nanpercentile(data, [2, 98])
        im = ax_wf.imshow(
            data,
            origin="lower",
            aspect="auto",
            interpolation="none",
            extent=[mdates.date2num(time_axis[0]), mdates.date2num(time_axis[-1]), frequencies[0], frequencies[-1]],
            cmap="plasma",
            vmin=vmin,
            vmax=vmax
        )

        # Plot GNSS emission frequency lines on the waterfall
        gnss_freqs = {
            'gps': [1575.42, 1227.60, 1176.45],
            'galileo': [1575.42, 1176.45, 1207.14, 1278.75],
            'glonass': [1602.0, 1246.0, 1202.0],
            'beidou': [1561.098, 1207.14, 1268.52, 1176.45]
        }
        freq_min, freq_max = min(frequencies[0], frequencies[-1]), max(frequencies[0], frequencies[-1])
        freq_labels = {}
        for name, color in legend_data:
            name_lower = name.lower()
            for const, freqs in gnss_freqs.items():
                if const in name_lower or (const == 'galileo' and 'gsat' in name_lower):
                    for f in freqs:
                        if freq_min <= f <= freq_max:
                            if f not in freq_labels:
                                freq_labels[f] = []
                            short_name = name.split('(')[0].strip()
                            freq_labels[f].append(short_name)

        time_end = mdates.date2num(time_axis[-1])
        time_range = mdates.date2num(time_axis[-1]) - mdates.date2num(time_axis[0])
        x_text = time_end + time_range * 0.02

        for f, names in freq_labels.items():
            ax_wf.axhline(y=f, color='cyan', linestyle='--', alpha=0.7, linewidth=1.0)
            label_str = ", ".join(names)
            ax_wf.text(x_text, f, label_str, color='black', va='center', ha='left', fontsize=7, clip_on=False)
        ax_wf.xaxis_date()
        ax_wf.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M:%S'))
        ax_wf.set_xlabel("Time (FITS)", fontsize=10)
        ax_wf.set_ylabel("Frequency [MHz]", fontsize=10)
        ax_wf.set_title("Waterfall — UIRAPURU", fontsize=11, fontweight='bold')
        plt.colorbar(im, ax=ax_wf, orientation="horizontal", label="Intensity [dB]", pad=0.12)
        
        # Share the X axis
        ax_traj.set_xlim(ax_wf.get_xlim())
        
        # Right panel: legend (satellite table)
        ax_legend.axis('off')
        
        if legend_data:
            # Create a table of satellites and colors
            table_data = [['Satellite', 'Color']]
            colors_for_table = [['white', 'white']]  # Header
            
            for name, color in legend_data:
                # Truncate long names
                display_name = name[:22] if len(name) > 22 else name
                table_data.append([display_name, '●'])  # Colored symbol
                colors_for_table.append(['white', color])
            
            # Create the table
            table = ax_legend.table(
                cellText=table_data,
                cellLoc='left',
                loc='center',
                colWidths=[0.80, 0.20]
            )
            table.auto_set_font_size(False)
            table.set_fontsize(11)
            table.scale(1, 2.8)
            
            # Set cell colors
            for i, (row_data, row_colors) in enumerate(zip(table_data, colors_for_table)):
                for j, (cell_text, cell_color) in enumerate(zip(row_data, row_colors)):
                    cell = table[(i, j)]
                    if i == 0:  # Header
                        cell.set_facecolor('white')
                        cell.set_text_props(weight='bold', color='black', fontsize=11)
                        cell.set_edgecolor('gray')
                        cell.set_linewidth(0.5)
                    else:
                        if j == 1:  # Color column
                            cell.set_facecolor(cell_color)
                            cell.set_text_props(color=cell_color, weight='bold', fontsize=16)
                            cell.set_edgecolor('gray')
                            cell.set_linewidth(0.5)
                        else:  # Name column
                            cell.set_facecolor('white')
                            cell.set_text_props(color='black', fontsize=10)
                            cell.set_edgecolor('gray')
                            cell.set_linewidth(0.5)
            
#            ax_legend.text(0.5, -0.05, f'Total: {len(legend_data)} satellites', 
#                          ha='center', va='top', fontsize=11, fontweight='bold',
#                          transform=ax_legend.transAxes)
#        else:
#            ax_legend.text(0.5, 0.5, 'No satellites\nno FOV', 
#                          ha='center', va='center', fontsize=12,
#                          transform=ax_legend.transAxes)
        
#        plt.suptitle(f'GNSS analysis - {len(legend_data)} satellites detected', 
#                     fontsize=13, fontweight='bold', y=0.98)
        plt.show()
