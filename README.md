# 2020_Lavickova_model
Code for **A self-regenerating synthetic cell model**, Lavickova B, Laohakunakorn N, and Maerkl SJ (2020) [https://doi.org/10.1101/2020.07.03.185900](https://doi.org/10.1101/2020.07.03.185900)

Tested on Julia 1.4.1, 1.4.2. 

**Usage:**

	git clone https://github.com/nadanai263/2020_Lavickova_model.git
1. If you have a local Julia installation: ```cd``` to the project directory and run

	(v.1.4) pkg> activate .
	(2020_Lavickova_model) pkg> instantiate
More info here: [https://julialang.github.io/Pkg.jl/stable/environments/](https://julialang.github.io/Pkg.jl/stable/environments/).

2. Alternatively you can run the code from a Docker container:

	docker pull nadanai263/2020lavickova

```cd``` to the main project directory and run

	docker run -p 8888:8888 -v "$PWD":/home/jovyan/ 2020lavickova
This mounts the repository's working directory to the container. More information can be found [here](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html).

3. You can also build the same Docker image from the Dockerfile: ```cd``` to the Dockerfile folder and run

	docker build -t lavjul .
where ```lavjul``` is the arbitrary name of the new Docker container. This builds a new container based on the latest jupyter/datascience-notebook image, which contains working distributions of Julia, Python, and R, and a Jupyter notebook installation. To start the container, ```cd``` to the main repo directory and run

	docker run -p 8888:8888 -v "$PWD":/home/jovyan/ lavjul

---

**Files:**

* ./notebooks - interactive Jupyter notebooks with examples of running the model (model_implement.ipynb) as well as generating all plots for the supplemental information (model_plots_*.ipynb). Linked to code in the scripts directory.

* ./scripts - model scripts. models.jl contains the ODEs, callbacks.jl contains the functions required to implement periodic dilution, solve.jl contains the ODE solver code, and run.jl is an example script which runs the model by calling the other scripts. 
