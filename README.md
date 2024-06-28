# uf_decoder

Modified union-find decoder for arbitrary CSS codes. The decoder makes a simplification as it avoids tracking of cluster boundaries and uses instead breadth-first graph traversal of the Tanner graph of the code. The obtained thresholds are on par with https://arxiv.org/abs/1709.06218, yet the decoding is faster due to the made simplifications.

## features
- Fast decoding of topological codes with periodic boundary conditions (Algorithms 1, 2, 3)
- Decoding of general qLDPC codes (Algorithm 4)
- Include qubit erasure/loss into the decoding

## usage

The ```C```-code can either be used via one of the examples in the folder ```examples``` or it can be interfaced from ```Python```. For the latter option, you find examples in the folder ```py_wrapper```. Before usage, the code needs to be compiled with the following commands:

```
$ meson setup build --buildtype=release
$ cd build/
$ ninja
```
For interfacing with ```Python```, you will use the library ```libSpeedDecoder.so``` via ```run_from_py.py```.

## source code
The implementation of Algorithms 1-3 can be found in ```src/decoder_main.c```. The implementation of Algorithm 4 is in ```src/decoder_ldpc.c``` with ```src/stabilizer_main.c``` containing the used Gaussian elimination with bitwise-logic.
