# uf_decoder

This repository contains union-find decoders for CSS codes. For the 2d and the (2+1)d toric code, the obtained thresholds are on par with the original union-find decoder from https://arxiv.org/abs/1709.06218. In contrast to the original decoder, our decoder grows clusters node-by-node instead of cluster-by-cluster for syndrome validation. It avoids tracking cluster boundaries and uses breadth-first graph traversal of the Tanner graph instead. Note that this repository is a preliminary version of the decoder and breaking changes may occur in future versions.

## features
- Fast decoding of topological codes with periodic boundary conditions (Algorithms 1, 2, 4)
- Decoding of general qLDPC codes (Algorithm 3)
- Include qubit erasure/loss into the decoding

## usage (on Linux systems)

The ```C```-code can either be used via one of the examples in the folder ```examples``` or interfaced from ```Python```. For the latter option, you can find examples in the folder ```py_wrapper```. Before usage, the code needs to be compiled. Assuming that a compiler such as ```gcc``` is present on the used system, this can be done with the following commands:

```
$ meson setup build --buildtype=release
$ cd build/
$ ninja
```
For interfacing with ```Python```, you will use the library ```libSpeedDecoder.so``` via ```run_from_py.py```. In the current implementation, data qubits and syndromes share the same data structure with syndromes being at the beginning of the array, followed by data qubits. We may change this in a later version to simplify the interfacing from ```Python``` and improve the time and space complexity.

## source code
The implementation of Algorithms 1, 2, 4 can be found in ```src/decoder_main.c```. The implementation of Algorithm 3 is in ```src/decoder_ldpc.c``` with ```src/stabilizer_main.c``` containing the used Gaussian elimination with bitwise-logic.
