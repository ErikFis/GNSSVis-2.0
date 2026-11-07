# Interactive notebooks

The notebooks are distributed without saved outputs, embedded images, widget
state, or execution counts, keeping the source package small and portable.

- `gnssvis_simple_demo.ipynb`: waterfall, baseline subtraction, and smoothing.
- `gnssvis_satellite_validation.ipynb`: GNSS satellite identification workflow.

Install the package with `python -m pip install -e .`, then open the notebooks
in your preferred Jupyter environment. Set the FITS path to your own local
observation. Observation data is deliberately excluded. Initialization and
trajectory analysis may download ephemerides and TLEs.

Before committing notebook edits, clear outputs and execution counts. Keep
plots and numerical output outside the source submission. Check the task
budget with `python tools/check_task_size.py`.
