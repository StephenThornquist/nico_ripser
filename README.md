# Analyses for Niko

Here's a demo script to parse the analyses and run `ripser` on the data.

Simply run `main.py` in an appropriate `Python` environment.

## Installation

I use `uv` for environment management, which is nice because it has lockfiles to
ensure reproducible version control. However, the package itself should build
just fine with just the `pyproject.toml` using `pip`. The two build strategies
are described here:

First clone the github repository onto your machine:

`git clone <>`

Go to the directory containing the repository while using the Python
interpreter of your choice (e.g. in the environment you want to install into) and type

### `uv`

```
uv run main.py
```

### `pip`

```
pip install .
```

This will ensure you have the dependencies you need to run this code.
The dependencies are quite minimal.

Then type

`python main.py`.

This will execute the demo function, which plots some EPG data
alongside its persistence diagrams, both during the pre-VR phase
and the during-VR phase.

## Data

These imaging data are collected from the protocerebral bridge of
female Drosophila melanogaster. They're fairly processed
before they get to the `.npz` files contained herein. The short version
is we:
- Register the imaging data
- Use a circular correlation analysis algorithm (https://github.com/MaimonLab/FourCorr/)
to identify segments of the protocerebral bridge (called glomeruli)
- Sum the fluorescence intensity within each glomerulus during each imaging frame
into a single number $F(t)$.
- Compute a $\Delta F/F_0(t)$ signal by subtracting a baseline estimate $F_0$ from each ROI
at every timepoint and then dividing by that $F_0$ (i.e. $\Delta F/F_0 = (F-F_0)/F_0$). In keeping with the standard
in the literature (though I intend to deviate from this standard at some point),
that $F_0$ value corresponds to the 5th percentile value of $F$ over the duration
of the recording.

The $\Delta F / F_0$ reported here are ordered anatomically from the left half
of the fly's head to the right half of the flies head by index (so the 0th row
of the array corresponds to the leftmost glomerulus, and the 15th row corresponds
to the rightmost glomerulus).

The VR position is collected by analyzing images of a patterned ball on which the
flies walk with a variant of an algorithm called `FicTrac`. We use the VR position
and orientation to compute the visual display, which is just a 12 degree UV-wavelength
vertical bar that rotates as the fly yaws through the virtual world. No features change
as the fly translates through space in these simple experiments.