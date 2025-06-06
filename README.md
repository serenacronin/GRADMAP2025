Location to store code for the 2025 GRAD-MAP project.

====================================================================================
# LOG ONTO A UMD ASTRO COMPUTER REMOTELY
====================================================================================

Step 1: Open a terminal on your laptop.

Step 2: Type in the following and hit enter:

    ssh username@ssh.astro.umd.edu

Step 3: Put in your password. Then, type:

    ssh computer_name

Step 4: Navigate to the the directory you want to work in by typing:

    cd /computer_name/username/

To list all files in this directory, type:

    ls

To make a new sub directory in this directory, type:

    mkdir directory_name

To navigate to a new directory, type:

    cd path_to_directory

or, to navigate backwards:

    cd ../path_to_directory

The above assumes you will be working on the computer 'jansky' (which we recommend for long computations).


====================================================================================
# SET UP A CONDA ENVIRONMENT (Mac & Linux)
====================================================================================

Step 0: Make sure you have Anaconda downloaded and ready to go: https://docs.anaconda.com/free/anaconda/install/

Step 1: Create the environment! Use a custom name, for example, "gradmap25":

    conda create --name gradmap25
    
Make sure you press 'y' when prompted.

Step 2: Activate your environment.

    source activate gradmap25

Step 3: Add the packages you need using the following commands:

    conda install python astropy jupyter matplotlib numpy pandas scipy spyder

All the latest versions should be fine, but if you want to download a specific version, this would look like (for example):
    conda install python=3.10.4 astropy=5.1

Make sure you press 'y' when prompted.

Step 4: To deactive your environment, type:

    source deactivate gradmap25

You must activate your environment every time you open a new terminal.


====================================================================================
# SET UP A CONDA ENVIRONMENT (Windows)
====================================================================================

Step 0: Make sure you have Anaconda downloaded and ready to go: https://docs.anaconda.com/free/anaconda/install/

Step 1: Create the environment! Use a custom name, for example, "gradmap25":

    conda create --name gradmap25
    
Make sure you press 'y' when prompted.

Step 2: Activate your environment.

    conda activate gradmap25

Step 3: Add the packages you need using the following commands:

    conda install python astropy jupyter matplotlib numpy pandas scipy

All the latest versions should be fine, but if you want to download a specific version, this would look like (for example):
    conda install python=3.10.4 astropy=5.1
    
Make sure you press 'y' when prompted. The specified versions are just to match what I have tested and run the fitting program on.

Step 4: To deactive your environment, type:

    deactivate gradmap25

You must activate your environment every time you open a new terminal.


====================================================================================
# WORKSPACE IDE
====================================================================================

I recommend using Spyder for this project, but you can also use Jupyter Notebooks. To open them, type into the terminal:

    spyder

and

    jupyter notebook

You will need separate terminals for each.
