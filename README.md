# uf_decoder

Union-find decoders for CSS codes. An explanation of all algorithms can be found here: https://arxiv.org/abs/2407.15988.

## features
- Fast decoding of topological codes.
- Decoding of general qLDPC codes.
- Include qubit erasure/loss into the decoding.
- Thresholds of the 2d, (2+1)d toric code are on par with https://arxiv.org/abs/1709.06218.

## getting started step 0) compilation (on Linux systems)
If a compiler such as ```gcc``` is present on the used system, compilation can be done with the following commands:

```
$ pip install meson ninja
$ meson setup build --buildtype=release
$ cd build/
$ ninja
```
After compilation, you can use the decoder from ```Python``` via the library ```libSpeedDecoder.so``` (see ```py_wrapper/py_decoder.py```).

## getting started step 1) examples / getting started
For using ```Python```, there are several example scripts in the folder ```py_wrapper/```:
- ```example_toric_code.py``` applies Pauli errors and erasures to a 2d toric code (periodic boundary conditions) and visualizes the decoding.
- ```example_surface_code.py``` applies Pauli errors and erasures to a 2d surface code (with boundaries) and visualizes the decoding.
- ```example_surface_code_simulate.py``` simulates the logical error rate for a 2d surface/toric (batch evaluation).
- ```example_bb_codes.py``` simulates the logical error rate for bivariate bicycle codes from https://arxiv.org/pdf/2308.07915 (batch evaluation)

When using ```Python``` for simulations with many repetitions for averaging, batch evaluation should be used. Otherwise, the library ```libSpeedDecoder.so``` is called very often which may affect the overall running time. Instead of using ```Python```, the ```C```-code can be used independently. To this end, there are similar examples in the folder ```examples/```.

## source code and algorithms
The repository contains all algorithms from https://arxiv.org/abs/2407.15988. In contrast to the original decoder from https://arxiv.org/abs/1709.06218, our decoders grow clusters node-by-node instead of cluster-by-cluster for syndrome validation, using breadth-first graph traversal of the Tanner graph and not tracking cluster boundaries. In https://arxiv.org/abs/2407.15988, we label different decoding algorithms with numbers 1-4. The implementation of Algorithms 1, 2, and 4 from can be found in ```src/decoder_main.c```. These algorithms can be used for decoding topological codes where every data qubit participates in no more than two parity checks. Algorithm 3 can be used for decoding general CSS qLDPC codes and its implementation can be found in ```src/decoder_ldpc.c``` where also code from ```src/stabilizer_main.c``` is used. Compared to Algorithm 2, Algorithm 3 is less performant in terms of speed and error threshold when applied to topological codes.
