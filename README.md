## 3D Monte-Carlo Glauber Model

This code package simulates high-energy nucleus-nucleus collisions in 3 dimension. 

---

### Physics Model



---

### How to use the code

The code can be compiled with CMake. To generate the executable, one can simply type

	cmake ..
	make
	make install

The code package also includes unit test modules for individual class components. To generate the unit test executables, one can build the program as,
	
	cmake ..-Dtest=ON
	make
	make install

All the unit test executables will be installed to the test/ folder. Under the test/ folder, one can run the bash script run\_all\_UnitTests.sh to execute all the unit test programs.

