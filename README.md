
# FastTemporalScanning.jl
*An algorithm for quickly identifying injection shifts that bring about desired transmission line temperature changes.*

## Motivation
Consider two scenarios one might encounter in transmission system operation:
1. The grid has a number of renewable generation sites whose output is subject to forecast inaccuracy. Some lines are operating near their temperature limits. The transmission system operator would like to know how likely it is for small forecast deviations to bring a line in the network to its temperature limit. If some lines could be brought to their temperature limits via small forecast deviations, what would those deviation patterns look like?
2. A transmission line is operating above its limit temperature, and the operator must act to relieve the overload. Given a set of controllable generators whose actions have various costs (weights), what is the cheapest way for the operator to dispatch generation to drive the line temperature back down to an acceptable level?

Though their physical interpretations differ, the two cases technically have the same input and desired output:

> Given network data and a set of generators or loads presumed to be variable or controllable, find the smallest set of variations over a specified time horizon that brings a transmission line to a certain final temperature. 

The Fast Temporal Scanning algorithm described here can answer this question by solving an appropriate quadratically-constrained quadratic program. The algorithm does this without relying on a black box solver. 

## Mathematical description
The following mathematical program provides a compact description of the motivating question:

$$
\begin{subequations}\label{eq:opt}
\begin{align}
\label{eq:opt-obj}\underset{\vec{r}}{\min} \quad & \vec{r}^\top \mathbf{Q} \vec{r} \\
\nonumber \text{subject to:} & \\
\label{eq:opt-temp} \sum_{k=1}^{n_t} \hat{\theta}_{i^*j^*}(t_k)^2 &= \frac{c_1}{c_4}\left(T_\text{final} - c_7\right) \\
\label{eq:opt-tempvars}\hat{\theta}_{i^*j^*}(t_{k}) &= \theta_{i^*j^*}(t_k)\sqrt{ (e^{c_1\bar{t}})^{n_t-k+1} - (e^{c_1\bar{t}})^{n_t-k} } \\
\label{eq:opt-flow} \sum_j Y_{ij} \theta_{ij,t_k} & = G_{i,t_k} + \vec{r}_{i,t_k}) - D_{i,t_k} \\[-6pt]
\nonumber &\qquad\qquad\qquad\quad~ \forall i \in 1... n_b,~k\in 1... n_t \\[6pt]
\label{eq:opt-conv} \vec{G}_{t_k} &= \vec{G}_{0,t_k} + \alpha_{t_k}\vec{g} ~~ \forall k\in 1\ldots n_t \\
\label{eq:opt-ref} \theta_{ref,t_k} & = 0 \quad\quad\quad\quad\quad~~ \forall k\in 1\ldots n_t
\end{align}
\end{subequations}
$$

The set of generators or loads assumed to be variable or controllable is represented by $\mathbf{r}$. The objective matrix $\mathbf{Q}$ contains weights or covariance coefficients encoding the "cost" of each variation and the statistical relationships between them, if any. The first constraint fixes the final temperature of a chosen line to some fixed value. The next constraint defines the auxiliary angle variables used in the first constraint by connecting them to voltage angle differences across the chosen line. The following constraint implements power balance, while the one after that distributes mismatch between active power generation and load across conventional generators (thereby implementing droop response). The last constraint defines an voltage angle reference node for the system.

Hopefully the notation is mostly intuitive. A few symbols bear further definition:
* Constants $c_1$, $c_4$, and $c_7$ come from the line temperature model. They are defined in terms of line parameters including length and per-unit resistance and conductor parameters including diameter and heat capacity. These values are inputs to the code.
* $t_k$ is the k-th time step.
* $\bar{t}$ is the length of the interval between any two time steps.
* $n_t$ is the number of time steps under consideration.
* $D_{i,t_k}$ represents all non-generation injection at node $i$ and time step $k$. In other words, the net injection at each node is equal to conventional generation $G$, demand and other fixed injections $D$, and changes in variable/controllable generators $\mathbf{r}$. 
* $\alpha_{t_k}$ is the active power mismatch at time step $t_k$.

## Using the code
First, start Julia in the `src` directory and run:

```julia
include("FastTemporalScanning.jl")
using FastTemporalScanning
```

Now you are ready to load some data and pass it to the algorithm.

### Loading data
#### 1. Network data
Much of the network data required by the algorithm is typical of transmission network analysis tools. For that reason, there is a one-line way to load many relevant parts of the input data from a Matpower case of your choice:

```julia
n = network_data("case_name")
```

(Enter `casenames()` to return a list of valid case names.) After you run the line above, `n` will be a dictionary with the following fields:
* `f`: Vector of integers containing the "from" bus for each line.
* `t`: Vector of integers containing the "to" bus for each line.
* `r`: Vector of floats containing the per-unit resistance of each line.
* `x`: Vector of floats containing the per-unit reactance of each line.
* `b`: Vector of floats containing the per-unit susceptance of each bus.
* `G`: Vector of floats containing nominal conventional generation for each bus.
* `k`: Vector of floats containing participation factors for conventional generators.
* `D`: Vector of floats containing nominal demand at each bus.
* `Y`: Admittance matrix (necessary?)
* `Sb`: System base MVA (in VA)

This data may also be extracted from one of Jon Martin's `System` Matlab structures. Suppose you have a file called `data.mat` in a directory called `data` alongside `src`. `data.mat` contains the structure `System` (and potentially other variables as well -- it does not matter). Then you can obtain the same dictionary `n` as above by running

```julia
n = network_data("data.mat")
```

Of course, the third option is to come up with your own means of building `d`. As long as you pay attention to units and all the vectors are the right length, you shouldn't have too much trouble.

#### 2. Thermal data
The heart of the fast temporal scanning algorithm its line temperature model, which maps temperature changes to active power flows on the line at each time step. This thermal model is based on the IEEE 738 standard, and is summarized briefly here. To a close approximation, line temperature may be represented as the solution to the initial value problem

$$
\begin{align}
\frac{dT_l}{dt} &= c_tT_l(t) + c_2.
\end{align}
$$

The solution to this simple IVP is

$$
\begin{align}
T_l(t) &= c_3e^{c_1t} - \frac{c_2}{c_1},
\end{algin}
$$

where $c_3$ is defined in terms of the initial temperature and the other two constants: $c_3 = T_l(t_0) + c_2/c_1$. The constants are derived in our journal article. For the purposes of running the code, the important thing is that they are straightforward functions of the following thermal parameters:

* `rij` [pu] is resistance
* `xij` [pu] is reactance
* `Sb` is system base MVA
* `mCp` [J/m-C] is line heat capacity
* `eta_c` [W/m-C] is conductive heat loss rate coefficient
* `eta_r` [W/m-C^4] is radiative heat loss rate coefficient
* `Tlim` [C] is highest allowable line temperature
* `length` [m] is line length
* `Tamb` [C] is ambient temperature (of air)
* `eta_s` [] is solar heat rate coefficient

The first three parameters come from the network data dictionary `n`. The next four come from conductor material data sheets. Length, of course, varies from line to line. The last two parameters depend on ambient and solar conditions.

As with network data, there are multiple ways to load this data in preparation for a fast temporal scan. If you wish to generate thermal data from line electrical parameters, you can use the dictionary `n` as follows:

```julia
t = thermal_data(n)
```

(This functionality is based on Jon Martin's code for estimating line thermal parameters, itself based on work by Mads Almassalkhi.) After running the above, `t` will be a dictionary with the following fields:
* `D`: Vector of floats containing conductor diameter for each line.
* `mCp`: Vector of floats containing heat capacity of each line.
* `Ilim`: Vector of float containing steady-state current limits for all lines.
* `eta_c`: Vector of floats containing conductive heat rate coefficients for all lines.
* `eta_r`: Vector of floats containing radiative heat rate coefficients for all lines.
* `eta_s`: Vector of floats containing insolation coefficients for all lines.
* `length`: Vector of floats containing line lengths in meters.

The second way to load thermal data is to declare a conductor name for each line. Suppose your network has three lines: the first two are "dove" and the last is "waxwing". You can obtain thermal data for these lines as follows:

```julia
l = ["dove";"dove";"waxwing"]
t = thermal_data(l)
```

(Enter `conductornames()` to return a list of valid conductor nicknames.)

#### 3. Forecast data
Though the network data contains nominal generator and load active power injection values, a major selling point of the Fast Temporal Scanning algorithm is its ability to incorporate changing forecast values across a time horizon. This part of data loading takes a bit more work than the other two. Suppose you wish to consider a network with two buses over three time steps, each 60 seconds in length. Bus 1 has a conventional generator that outputs a constant 1 pu. Bus 2 has no conventional generator. There is a wind node at Bus 1 that is forecast to increase by 0.1 pu every minute until the end of the horizon from an initial output of 0.2 pu. Both nodes have constant loads of 0.5 pu. Here is how you would enter that information:

```julia
f["G"] = [[1.0;0.0];[1.0;0.0];[1.0;0.0]]
f["D"] = [[-0.3;-0.5];[-0.2;-0.5];[-0.1;-0.5]]
f["int_length"] = 60
```

Three things to note. First, the data format is a "vector of vectors". `f["G"]` is a vector of length 3 (the number of time steps), but each each element is itself a vector of length 2 (the number of buses). Second, note that the forecast wind output is encoded in `f["D"]` as negative demand. Why not just add generation in there as well? Because if we did that, there would be no way to implement droop response (where conventional generators divvy up active power mismatch)! Finally, note that `f["int_length"]` contains the interval length in *seconds*.

#### 4. Deviation data
The ultimate goal of the Fast Scanning Algorithm is to orchestrate a set of variable generators or loads in order to influence a transmission line's temperature. Thus, the algorithm needs to know where each variable injection is, what its cost is, and what covariance there is between variable injections. Consider two cases.

**Case 1.** Suppose the variable injections are due to uncertainty in wind farm output, and the objective is to find the smallest set of deviations that bring a certain line in the network to its temperature limit. The network has two buses. There are two wind farms at Bus 1 and one at Bus 2. Covariance between the Bus 1 wind farms is 0.5. The wind farm at Bus 2 is assumed to be independent of the Bus 1 sites. Here is how one would enter the deviation data:

```julia
d["bus"] = [1;1;2]
d["matrix"] =  [1.0 0.5 0.0;
		0.5 1.0 0.0;
		0.0 0.0 1.0]
```

**Case 2**. Now suppose the variable injections correspond to controllable generators, and the objective is to find an efficient set of dispatch instructions that reduces a line's temperature. The system has three buses. Bus 2 has one controllable generator, and Bus 3 has two. The Bus 3 generator is assumed to be twice as expensive to manipulate as the other two. Here is how one would enter the deviation data:

```julia
d["bus"] = [2;2;3]
d["matrix"] =  [1.0 0.0 0.0
		0.0 1.0 0.0
		0.0 0.0 2.0]
```

Note that the matrix is diagonal in this case; generator deviations are assumed to be independent, so off-diagonal terms are zero.

## Reading the code
The project directory is organized like a typical Julia package. All source code resides in `src/`, which stands alone. Data files used to put the source code through its paces are in `data/`.


## Organization
I will go through my notebooks in chronological order, skipping or briefly mentioning work that has been superseded.

I copied the nbs directory. I will now combine all the information into single narrative, deleting the nb copies as I go.

* 2015-02-28. Getting started. Early versions of a few matrix-building functions. Typeset description of how things are kept track of in matrix form.

* 2016-02-17. Understanding conductor assignment. Limit current and base voltage go in. Diameter, aluminum and steel masses, resistance [Ohms/m], and bundle info come out. This function is called by `ConductorModel_setup`. I combined the two files since their separation was more an artifact of MATLAB than logic.
