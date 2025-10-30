# ice-melt-models

![CI](https://github.com/oscarjtg/ice-melt-models/actions/workflows/ci.yml/badge.svg)

Ice melt models for interactions between ice sheets and polar oceans!

Based on models presented in [Jenkins et al., (2010)](https://doi.org/10.1175/2010JPO4317.1).

For explanations, see the Jupyter notebook located at `examples/melt_model_explanations.ipynb` which explains the motivation for this module and the origin of the melt models.

## Installation instructions

### Method 1 (quick)

Clone this repository, copy the file `icemeltmodel.py` into the directory with your Pythons scripts and import the module
```python
from icemeltmodels import ThreeEquationMeltModel

model = ThreeEquationMeltModel(density_water=1025.)
melt_rate = model.melt_rate(current_speed=0.5, 
                            temperature=2.0,
                            salinity=35.0,
                            pressure=0.0)
```

Also works with `numpy` arrays!

### Method 2 (takes longer, but useful if using this module across multiple projects)

Clone this repo into the `user site' directory of your Python distribution. You can find this directory by running
```
python -m site --user-site
```
Then run 
```
cd <your user site directory>
```
and 
```
git clone https://github.com/oscarjtg/ice-melt-models.git
cd ./ice-melt-models
```

Then (preferably within an active virtual environment) you can install the package by navigating to the installed directory, which contains a `pyproject.yaml` file 
(you should already be in this directory if you ran the `cd` command above) 
and typing into your terminal
```
pip install .
```
Alternatively, from another directory you can run
```
pip install <path-to-ice-melt-models-directory>
```

Then the pip package manager should install the package and its dependencies to your virtual environment (though NB I have not added the dependencies for the example Jupyter notebook).

You can then import the models as before
```python
from icemeltmodels import ThreeEquationMeltModel
```

## TODO:
* Write automated unit tests
* Doxygen documentation