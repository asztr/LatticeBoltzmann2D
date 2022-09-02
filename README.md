# 2D Lattice-Boltzmann

Interactive simulation of a 2D incompressible viscous fluid using the Lattice-Boltzmann D2Q9 method.

<div align="center"

<video src="https://user-images.githubusercontent.com/10238412/188178624-deeab462-2ffe-49c9-a36c-99017b388549.mp4" type="video/mp4;"></video>

</div>

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
