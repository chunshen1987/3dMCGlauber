## 3D Monte-Carlo Glauber Model

This code package simulates high-energy nucleus-nucleus collisions in 3 dimension. The code uses C++11 standard. 

---

### Physics Model

The spatial configuration of the nucleons inside the nucleus is sampled accroding the Woods-Saxon distribution. 

Strings that carry energy will be produced at collision point. The string will evolve accroding to deceleration dynamics, which decelerate the rapidity of the colliding nucleons pair. 

Optionally, valence quarks can be used as the participant degree of the freedom for the collisions. Their incoming rapidities are sampled according to the parton distribution function.

A list of the produced strings (their spatial and momentum information) will be output to a file at the end of the simulations. 

---

### How to use the code

The code can be compiled with CMake. To generate the executable, one can simply type

	mkdir build
	cd build
	cmake ..
	make
	make install

The code package also includes unit test modules for individual class components. To generate the unit test executables, one can build the program as,
	
	mkdir build
	cd build
	cmake .. -Dunittest=ON
	make
	make install

All the unit test executables will be installed to the test/ folder. Under the test/ folder, one can run the bash script run\_all\_UnitTests.sh to execute all the unit test programs.

### Reusing the minimum-bias centrality table

When `cenMin` and `cenMax` select a centrality range narrower than 0–100%,
the program calibrates the cuts from 50,000 inelastic events and saves the
sorted string counts (`Nstrings`) in `minimum_bias_collisions.dat` in the
working directory. Later runs automatically load a matching table and can use
it for different centrality ranges and random seeds at the same collision
energy. Set `minBiasTableFile`
in the input or pass `minBiasTableFile=/path/to/table.dat` to share a table
across run directories or keep separate calibrations.

The program regenerates missing, invalid, or incompatible tables. Choose a
fresh table path after replacing underlying nuclear or parton input tables.
For reproducible production events, keep the same saved table and run seed.
The run that generates the table consumes random numbers during calibration,
so its production events can differ from a run that loads the table with the
same seed. A 0–100% run skips calibration entirely.
