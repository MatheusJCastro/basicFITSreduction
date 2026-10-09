# Basic Data Reduction

*Written by: Matheus J. Castro*

Simple package for calibrating, aligning, and combining astronomical FITS images.

## Installation

```bash
pip install git+https://github.com/MatheusJCastro/basicFITSreduction.git
```

## Usage

Configure and run:

```python
from reduce_raw import ReduceRaw

reduction = ReduceRaw()
reduction.dir_name = "2026-10-06/group2"
reduction.obj_name = "NGC7635"
reduction.combine_filters = ["R", "G", "B"]
reduction.combine = True
reduction.reduce_bias = True
reduction.reduce_flat = True
reduction.run()
```

You can also run the example with `python -m reduce_raw`.

Parts of this project's organization and refactoring were completed with the
assistance of GitHub Copilot.

## License

This project is distributed under the GNU Affero General Public License v3.0.
