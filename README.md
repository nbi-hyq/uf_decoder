# uf_decoder

Union-find decoders for CSS codes (see https://arxiv.org/abs/2407.15988). For the 2d and the (2+1)d toric code, the obtained thresholds are on par with the original union-find decoder from https://arxiv.org/abs/1709.06218. In contrast to the original decoder, our decoders grow clusters node-by-node instead of cluster-by-cluster for syndrome validation, using breadth-first graph traversal of the Tanner graph and not tracking cluster boundaries. (This repository is a preliminary version, and breaking changes may occur in future versions.)

## features
- Fast decoding of topological codes.
- Decoding of general qLDPC codes.
- Include qubit erasure/loss into the decoding.

## compilation / getting started (on Linux systems)
Assuming that a compiler such as ```gcc``` is present on the used system, compilation can be done with the following commands:

```
$ meson setup build --buildtype=release
$ cd build/
$ ninja
```
After compilation, you can use the library ```libSpeedDecoder.so``` from ```Python``` as done in ```py_wrapper/py_decoder.py```.

## examples
For using ```Python```, there are several example scripts in the folder ```py_wrapper/```:
- ```example_toric_code.py``` applies Pauli errors and erasures to a 2d toric code (periodic boundary conditions) and visualizes the decoding.
- ```example_surface_code.py``` applies Pauli errors and erasures to a 2d surface code (with boundaries) and visualizes the decoding.
- ```example_bb_codes.py``` simulates the logical error rate for bivariate bicycle codes from https://arxiv.org/pdf/2308.07915 (batch evaluation)

Note that batch evaluation should be used for simulations with many repetitions for averaging. Otherwise, the library ```libSpeedDecoder.so``` is called very often which may affect the overall simulation time.

Instead of using ```Python```, the ```C```-code can be used independently. To this end, there are similar examples in the folder ```examples/```.

## source code and algorithms
The repository contains all algorithms from https://arxiv.org/abs/2407.15988. The implementation of Algorithms 1, 2, 4 from can be found in ```src/decoder_main.c```. The implementation of Algorithm 3 is in ```src/decoder_ldpc.c``` with ```src/stabilizer_main.c``` containing the used Gaussian elimination with bitwise-logic.
