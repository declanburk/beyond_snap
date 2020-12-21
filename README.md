# beyond_snap
A package for generating minimum-snap trajectories.
## Installation
A prerequisite is Eigen 3. Linux users can install through apt.
```
sudo apt-get install libeigen3-dev
```
Build the package using cmake.
```
cd build && cmake .. && make
```
Test everything works correctly with the example executable `test_fast_snap`. You can visualise the results with the `read_coefficients.py` script.
```
./test_fast_snap
cp min_snap_coefficients.csv ../utils && cd ../utils
python3 read_coefficients.py
```
## Using the SnapMinimiser class
<!-- The SnapMinimiser class can be used to generate minimum snap trajectories. Create the object and specify the number of waypoints in the trajectory as well as setting a debug flag (set to true for verbose output).
```
int K = 5;
bool debug = false;

SnapMinimiser test_snap(K, debug);
```
Parametrise the trajectory by setting the time and position associated with each waypoint using `set_waypoint()`.
```
// Initialiase to origin at time zero
test_snap.set_waypoint(0, 0.0, 0.0);
// Pass through first waypoint at x=1.5 at t=0.8
test_snap.set_waypoint(1, 0.8, 1.5);
```
The time and position at each waypoint must be specified to generate the trajectory. Call `solve_minimisation()` to calculate the spline coefficients. The coefficients are stored in array of `VectorP<s>` objects of length *k*, so you will need to use `get_P()` to access them.
```
test_snap.solve_minimisation();
// Pointer to array with spline coefficeints. Here I'm using s=5.
VectorP<5> *P = test_snap.get_P();
```
Each element of the array is a `VectorP<s>` object, which contains the coefficients of the corresponding polynomial segment. The coefficients are indexed from order *0* to *2s*. -->