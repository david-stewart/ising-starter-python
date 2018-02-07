# 2D Ising Model in Python w/ a C-Library for Matrix Logic

Written by Surya Dutta '18 | Original Matlab code written by Jed Thompson '17

Modified to use a C-Library for Matrix Logic by David Stewart

## Regarding the C-Library:
A small c-library `ising_matrix.c` has been written which implements the 2D ising
model for each step.

In order to build the library on your machine, use the following command (if you use
OSX): 

* `clang -shared -Wl,-install_name,ising_matrix.so -o ising_matrix.so -fPIC ising_matrix.c`

If you use Linux:

* `gcc -shared -Wl,-soname,ising_matrix -o ising_matrix.so -fPIC ising_matrix.c`

(For windows... to be updated...)

Once the c-library compiled it is a "shared object" `.so` which python can now call, using
the `ctypes` modules.

The `main.py` and `ising.py` modules have been updated accordingly. The `annealing.py` code
is unaffected.

A few notes:

* Passing any parameters to the C-library will required passing C-language types. 
For example, when calling the function `allocate(int)` in the C-library,
python must pass a C-language int to the function call. Assuming that python
has named the library `c_matrix`, the call int python(with python int `n`) would be:

    `c_matrix.allocate(c_int(n))`

* Python assumes that the return values from calls to the shared
C-library will be integers and will error if they are not.
Therefore the Python code must be told expressly what
the return type is for each C-library function used.
For example, the C-library function with signature `float get_E()`
returns a C-language float. Therefore, before using this function in the Python
code, the following line must first be used:

    `c_matrix.get_E.restype = c_float`

When python gets a C-language float, it can convert then to a 
Python float object with the `float()` function. Therefore, to get the
value of `E` from the C-library into Python:

    `E = float(c_matrix.get_E())`

## New Features
* Command line interface to input simulation parameters
* Better interpretability with pythonic features (e.g. list comprehension)
* Modular codebase for easy changes and experimentation
* Fancy progress bars with time estimations
* Complete error handling with progress save
* ...and more coming!


## Getting Started

### GUI-based approach (Windows/Mac/Linux)

If you prefer to use GUIs as opposed to the command line, this section is for you!

#### Managing your Git Repository (optional, but recommended)

If you would like to use Git Version Control in your team to collaborate on and keep backups of your code, great! If not, no worries - just download the files here and follow the instructions below. Version control is always highly recommended.

There are plenty of great GUIs for Git. My personal favorite is [Github Desktop](https://desktop.github.com/).

If you are new to Git and want to learn more about version control, visit [this website](https://programminghistorian.org/lessons/getting-started-with-github-desktop) for a great primer on version control, Git, and Github Desktop.

The first step is to make a Github account and fork this repository (click on `Fork` in the top right). This will create a copy of this code onto your own account. Now you can follow the instructions for your respective GUI to clone this repository (download the files locally), and start working with the simulation!

#### Installing and using Python

If you don't have Python installed yet, I would highly recommend using the **Anaconda distribution** to install Python 3. You can find the installation instructions [here](https://docs.anaconda.com/anaconda/install/)

Once this is installed on your computer, you will have Python 3 ready to go, as well as important packages like NumPy and SciPy. You can view these packages and install new ones using the [Anaconda Navigator](https://docs.anaconda.com/anaconda/navigator/) (need to install this separately).

In order to edit and run your code, I would recommend [Spyder (Scientific PYthon Development EnviRonment)](https://pythonhosted.org/spyder/) (I know, horrible acronym, but the IDE makes up for it). It should be really easy to edit your code and run it through this environment. **NOTE**: The IPython shell in Spyder does not support nested progress bars, so you will only see one when you run the simulation. In order to see both, you will need to change your run configuration to run in a normal Python shell.

Another optional but cool program you can use is [Jupyter Notebook](http://jupyter.org/) (comes pre-installed with Anaconda). These notebooks support Python code, as well as Markdown and LaTeX, so you can keep all of your code organized and easily testable (hint: use for easier data analysis!). You should be able to open this through Anaconda Navigator.

These are just recommendations - there are plenty of other GUI-based applications for Python development out there (like Enthought Canopy). If you have time, do some playing around and see what you like!

### Command line (Mac/Linux)

1. Fork this repository to your own account

2. Navigate to the folder you would like to use, then use:
  ```bash
  git clone git@github.com:{{your-github-username}}/ising-starter-python ising && cd ising
  ```

3. If you are using conda (recommended), use this to install required packages:
  ```bash
  conda install --yes --file requirements.txt
  ```

  If you have a standalone version of Python 3 installed and are using pip, use this instead (you may need to be a superuser to install Pip packages):
 ```bash
  pip install -r requirements.txt
  ```

4. The code should be ready to run! Use this to run the simulation, and it should save the results automatically to an auto-generated data folder:
  ```bash
  python main.py
  ```

5. If you want to change other parameters of the simulation, you can use:

  ```
  python main.py --help

  2D Ising Model Simulation
  Usage: main.py [OPTIONS]

  Options:
    --t_min FLOAT           Minimum Temperature (inclusive)
    --t_max FLOAT           Maximum Temperature (inclusive)
    --t_step FLOAT          Temperature Step Size
    --n INTEGER             Lattice Size (NxN)
    --num_steps INTEGER     Total Number of Steps
    --num_analysis INTEGER  Number of Steps used in Analysis
    --num_burnin INTEGER    Total Number of Burnin Steps
    --j FLOAT               Interaction Strength
    --b FLOAT               Applied Magnetic Field
    --flip_prop FLOAT       Proportion of Spins to Consider Flipping per Step
    --help                  Show this message and exit.
  ```

  This will list all of the parameters you can change. For example, if you run `python main.py --b=0.5 --flip_prop=0.2`, the simulation will add a magnetic field of 0.5T and increase the flip proportion to 0.2. You can also edit the default parameters directly in the `main.py` file.


### Powershell (Windows)

Coming Soon!

## Understanding the Simulation

There are three important python files in this simulation: `main.py`, `ising.py`, and `annealing.py`

`main.py` is the file you have to run for the simulation. The code in this file takes in the input parameters, runs the Ising model for each temperature step, gets the relevant data, saves it, and gives you a set of nice plots at the end. This is a lot, so we've broken this down into different functions to make it easier to understand/change. Here are the two most important ones:

* `run_simulation`: takes in all the input variables and runs the simulation.

* `calculate_and_save_values`: takes in the energy, magnetization, and spin values from the Ising code, calculates the appropriate statistical values, and saves them to a CSV file. **This is where you should implement code to calculate the other values you are interested in**.

`ising.py` calculates the Ising model at a certain temperature

## Multi-processing

If you feel adventurous and want to use multiprocessing (running the simulation on multiple cores), we got you covered! Use the `main-multiprocessing.py` file to get started. This is identical to the `main.py` file, except there is an additional option for indicating how many processes you want to run.

For example, if you want to run your code on four processes (good place to start since most CPUs have four cores), use the following command:

```
python main-multiprocessing.py --processes=4
```

There isn't a clean way to add progress bars to this, so instead, the program will output the temperatures it is currently computing, as well as the temperatures it is finished with.

**NOTE:** Because of the asynchronous conditions of this program, it may not write your data to the CSV file in the right order. Make sure you sort by temperature before analyzing the data.

## To-Dos

* Add Windows Powershell instructions
* Add more info on Simulation structure
* Optimize (make code faster)

## Acknowledgements

Coming Soon!
