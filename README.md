# G.729AB multichannel codec optimized for TI C6000 DSP
Code was based on ITU-T reference implementation
* Highly optimized
* Added multichannel support and usable interface
* Added self testing
* Passing all test vectors

## How to use
Simply integrate it as part of your code or compile as an independent library (recommended).
Inlude g729itu.h in your code and use.

### Notes
* Use maximum optimization level in your CCS
* Enable caching
* Align context memory
