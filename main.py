from typing import Any
from pathlib import Path
import numpy as np
import ripser
from matplotlib.axes import Axes
from matplotlib.collections import QuadMesh
import matplotlib.pyplot as plt

def to_phase(arr : np.ndarray) -> np.ndarray[Any, complex]:
    """
    Takes data of the form ROIs x time and returns a
    phase array of shape (time, ) whose complex angle
    corresponds to the "phase" of the data
    
    ## Arguments
    `arr` : `np.ndarray`
        Shape (n_rois, time)

    ## Returns

    `phase` : `np.ndarray`
        Shape (time,)

    """

    circ = np.exp(1j*np.linspace(-np.pi, np.pi, arr.shape[0]))

    # Because the EPGs in the bridge are actually split, this relationship
    # should be [0, 2, 4, 6, 8, ... 1, 3, 5, ...]
    circ = np.array(circ[::2].tolist() + circ[1::2].tolist())

    return circ @ arr

def plot_imaging_heatmap(time : np.ndarray, data : np.ndarray, ax : Axes) -> QuadMesh:
    """
    Plots the EPG imaging data in the standard convention: time axis on the x
    coordinate, EPG rois ordered on the vertical component, shading proportional
    to the ΔF/F

    ## Arguments

    - `time`
        The time axis of the data, x axis, of shape `(time,)`

    - `data`
        The fluorescence data itself, of shape `(n_rois, time)`

    - `ax`
        The matplotlib `Axes` on which to plot.

    ## Returns

    - `pcolor : QuadMesh`
        The output of the `pcolormesh` call that produces the heatmap
    """

    # modify the array `time` so that each pixel is centered on its appropriate timepoint
    time = time.astype(float)
    time[1:] -= np.diff(time)/2
    time[0] -= np.mean(np.diff(time))/2
    time = np.insert(time, 0, time[0] - np.mean(np.diff(time)))

    xx, yy = np.meshgrid(time, np.linspace(-np.pi, np.pi, data.shape[0] + 1))

    pc = ax.pcolormesh(
        xx, yy, data,
        cmap = 'Blues', vmin = 0, vmax = 2.0
    )

    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    return pc

def plot_phase_and_heading(image_timestamps : np.ndarray,
                           phase : np.ndarray,
                           vr_timestamps : np.ndarray,
                           vr_heading : np.ndarray,
                           ax : Axes,
                           subtract_offset : bool = True,
                           ) -> None:
    
    # subtract out the offset between the phase
    # and the VR heading so that they overlay
    # (the offset is arbitrary from fly to fly)
    if subtract_offset:
        # upsample to the vr timestamp positions
        phase_interp = np.exp(1j*np.interp(
            vr_timestamps,
            image_timestamps,
            np.unwrap(np.angle(phase)),
        ))
        offset = np.mean(phase_interp / np.exp(1j*vr_heading))
        phase /= offset

    ax.plot(
        vr_timestamps,
        vr_heading,
        'o',
        markersize = 1,
        color = '#000000',
        label = 'VR heading',
    )

    ax.plot(
        image_timestamps,
        np.angle(phase),
        'o',
        markersize = 1,
        label = 'EPG phase'
    )

    ax.legend(
        loc = 'upper right',
        fontsize = 8,
        markerscale = 2,
    )
    pass

def plot_demo_file(file : Path, downsample_by : int = 1, subtract_offset : bool = True):
    """
    Demonstration of how to plot a few variables
    of the demo files. Produces `matplotlib` plots
    of the data stored in the `.npz` at `file`.

    ## Arguments

    - `file` : `Path`
        Path to the file to analyze
    """
    file_dict = np.load(file, allow_pickle = False)

    # ΔF/F data in ROI x time format
    dfof = file_dict['dfof']
    # Timestamps of each frame in epoch time
    image_timestamps = file_dict['image_timestamps']
    # Orientation in VR space of the imaged fly
    vr_heading = file_dict['vr_heading']
    # Timestamps of the behavior data
    vr_timestamps = file_dict['vr_timestamps']
    # VR position in x, y space as a complex number (x + iy)
    # not used here but you can if you want!
    vr_position = file_dict['vr_position'] #noqa
    # Timestamp of the moment the visual environment is activated
    vr_on_time = file_dict['bar_on_time']

    # Separate views of the data only where the VR is off and
    # the VR is on
    dfof_pre_vr, dfof_during_vr = (
        np.take(dfof, np.where(image_timestamps < vr_on_time)[0], axis = -1),
        np.take(dfof, np.where(image_timestamps >= vr_on_time)[0], axis = -1),
    )

    f, x = plt.subplots(nrows = 2, figsize = (8, 3), sharex = True)
    
    plot_imaging_heatmap(
        image_timestamps,
        dfof,
        x[0]
    )

    phase = to_phase(dfof)

    plot_phase_and_heading(
        image_timestamps,
        phase,
        vr_timestamps,
        vr_heading,
        x[1],
        subtract_offset = subtract_offset,
    )

    xticks = np.linspace(
        image_timestamps[0],
        image_timestamps[-1],
        4
    )

    x[1].set_xticks(
        xticks,
        ((xticks - image_timestamps[0])/1e9).astype(int)
    )
    x[1].set_xlabel('Time into imaging (sec)')

    f_rip, x_rip = plt.subplots(
        nrows = 1, ncols = 2, figsize = (8,4), sharex=True, sharey=True)

    r_pre = ripser.ripser(
        dfof_pre_vr.T[::downsample_by],
        metric = 'correlation',
    )

    r_post = ripser.ripser(
        dfof_during_vr.T[::downsample_by],
        metric = 'correlation',
    )

    for ax_idx, r in enumerate([r_pre, r_post]):
        for c_idx, dim in enumerate(r['dgms']):
            color = 'C' + str(c_idx)
            for (pair_x, pair_y) in dim:
                x_rip[ax_idx].plot(pair_x, pair_y, 'o', markersize = 1, color = color)
        titles = ['Pre VR', 'During VR']
        x_rip[ax_idx].set_title(titles[ax_idx])
    plt.show()

def main():
    plot_demo_file('data/Fly_1/imaging_data.npz',downsample_by=10)
    # print("Hello from nico-ripser!")

if __name__ == "__main__":
    main()
