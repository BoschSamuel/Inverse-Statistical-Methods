## Inverse Statistical Methods

Inverse problems in statistical physics are motivated by the challenges of “big data” in different fields, in particular experiments in biology. In inverse problems, the usual procedure of statistical physics needs to be reversed: Instead of calculating the observables over time on the basis of model parameters, we seek to infer parameters of a model based on observations.

In this project, I will focus on the inverse of the Potts model. The Potts model, a generalization of the Ising model, is a model of interacting spins. The reason of this is theoretical (maximum entropy) and practical: it is the simplest model that can explain two-point statistics in data. Also, it is widely used in the community. For example, the reconstruction of neural and gene regulatory networks, and the determination of the 3D structure of proteins [1].

# Running the code
The main simulation is "main_simulation.cpp"
You can compile it with "g++ -std=c++11 main_simulation.cpp" and run it with "./a.out"

It creates several output files (C.txt, Autocorrelation.txt and Energy_vs_time.txt), which can later be used for analysis.

For checking the accuracy of the main simulation, you can run "Analysis.ipynb" and get the correlation coefficient.

For making initial plots of a simulations, just run "Results.ipynb"


# Project description
Mentors: prof. Paolo De Los Rios Dr. Stefano Zamuner

Inverse statistical methods and pseudolikelihood approximation

What is the aim of the project?

	a) Simulate a Potts model with pairwise and triplet interactions
	b) Code an inference system based on Pott's model from scratch and test it.
	c) The test will be done by inferring the parameter of a Pott's model, starting from an equilibrium ensemble of configurations, simulated using Monte Carlo techniques.
	d) Other tests may involve inferring pairwise couplings of a protein family: Direct Coupling applications (DCA) [2]
	e) In prof. De Los Rios' group, a new technique has been devised to infer higher order interactions from data. For this reason, the code should be devised in such a way that it is easy to extend to this more general case.
	f) If there is enough time, we would test the new techniques to infer the parameters of a simulated Pott's model in which triplet interactions are present.


New skills which will be obtained through the project:

	a) C++ (so far, I only have experience in Python, Matlab and some basics in C)
	b) Implementation of a simple Metropolis algorithm
	c) Block analysis (or other techniques) to check the convergence of the simulation
	d) Theory of inverse Pott's model and maximum likelihood modeling
	e) Code for function minimization/maximization
	f) Code and theory of neural networks with several hidden layers


Further details:

Tremendous efforts have been made to determine the three-dimensional structure of proteins. A linear amino acid chain folds into a convoluted shape [3,4], the folded protein, thus bringing amino acids into close physical proximity that are separated by a long distance along the linear sequence. The three-dimensional structure of a protein determines its physical and chemical properties, and how it interacts with other cellular components: broadly, the shape of a protein determines many aspects of its function. Protein structure determination relies on crystallizing proteins and analyzing the X-ray diffraction pattern of the resulting solid. Given the experimental effort required, the determination of a protein’s structure from its sequence alone has also been key challenge to computational biology for several decades. The computational approach models the forces between amino acids in order to find the low-energy structure a protein in solution will fold into. Depending on the level of detail, this approach requires extensive computational resources.

References:
	[1] Ernst Ising. Beitrag zur theorie des ferromagnetismus. Zeitschrift für Physik, 31(1):253–258, 1925.
	[2] Renfrey Burnard Potts. The mathematical investigation of some cooperative phenomena, 1951.
	[3] M Newman and G Barkema. Monte carlo methods in statistical physics chapter 1-4. Oxford University Press: New York, USA, 1999.
	[4] Faruck Morcos, Andrea Pagnani, Bryan Lunt, Arianna Bertolino, Debora S Marks, Chris Sander, Riccardo Zecchina, José N Onuchic, Terence Hwa, and Martin Weigt. Direct-coupling analysis of residue coevolution captures native contacts across many protein families. Proceedings of the National Academy of Sciences, 108(49):E1293–E1301, 2011.
	[5] Edwin T Jaynes. Information theory and statistical mechanics. Physical review, 106(4):620, 1957.
	[6] Edwin T Jaynes. Information theory and statistical mechanics. ii. Physical review, 108(2):171, 1957.
	[7] T Plefka. Convergence condition of the tap equation for the infinite-ranged ising spin glass model. Journal of Physics A: Mathematical and general, 15(6):1971, 1982.
	[8] Antoine Georges and Jonathan S Yedidia. How to expand around mean-field theory using high-temperature expansions. Journal of Physics A: Mathematical and General, 24(9):2173, 1991.
	[9] Francis Galton. Regression towards mediocrity in hereditary stature. The Journal of the Anthropological Institute of Great Britain and Ireland, 15:246–263, 1886.
	[10] Karl Pearson. Notes on regression and inheritance in the case of two parents proceedings of the royal society of london, 58, 240-242, 1895.
	[11] Sarah Boslaugh. Statistics in a nutshell: A desktop quick reference. " O’Reilly Media, Inc.", 2012.
	[12] J. Taylor. Introduction to Error Analysis, the Study of Uncertainties in Physical Measurements, 2nd Edition. University Science Books, 1997.
