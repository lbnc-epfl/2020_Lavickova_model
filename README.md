# 2020_Lavickova_model
Code for **A self-regenerating synthetic cell model**, Lavickova B, Laohakunakorn N, and Maerkl SJ (2020)
doi: [https://doi.org/10.1101/2020.07.03.185900](https://doi.org/10.1101/2020.07.03.185900)

Tested on Julia 1.4.1, 1.4.2. 

Usage:
	git clone https://github.com/nadanai263/2020_Lavickova_model.git
If you have a local Julia installation: ```cd``` to the project directory and run
	(v.1.4) pkg> activate .
	(2020_Lavickova_model) pkg> instantiate
More info here: [https://julialang.github.io/Pkg.jl/stable/environments/](https://julialang.github.io/Pkg.jl/stable/environments/).

Using Docker: ```cd``` to the Dockerfile folder and run
	docker build -t lavjul .
where ```lavjul``` is the arbitrary name of the new Docker container. 



