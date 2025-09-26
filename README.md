# ice-melt-models
Ice melt parameterisations for ice sheet-ocean interfaces.

For explanations, see the Jupyter notebook located at examples/melt_model_examples.ipynb .

## Installation instructions

You can run this locally by cloning this repo into the `user site' directory of your Python distribution. You can find this directory by running
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

Then (preferably within an active virtual environment) you can install the package by navigating to the installed directory (which should be done already following the instructions above) and typing
```
pip install .
```
Alternatively from another directory run
```
pip install <path-to-ice-melt-models-directory>
```

Then the pip package manager should install the package and its dependencies to your virtual environment (though NB I have not added the dependencies for the example Jupyter notebook).

You can then access the models by (e.g.)
```
from icemeltmodels import ThreeEquationMeltModel
```

