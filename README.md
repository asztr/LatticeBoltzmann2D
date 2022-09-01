# 2D Lattice-Boltzmann

Interactive simulation of a 2D incompressible viscous fluid using the Lattice-Boltzmann D2Q9 method.

### Requirements
```
qmake, libfftw3, libQt
```

### Compilation
```
$ qmake .
$ make
```

### Usage
```
$ ./lbe2d

Left-click: move the fluid.
Right-click: build a rigid block. 
Press "d" to change the visualization mode.
```

### References
<span id="1">[1]</a> R. Bridson, _Fluid Simulation for Computer Graphics_, A.K.Peters (2008).
