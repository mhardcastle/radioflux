# radioflux
Measuring radio flux density with the ds9

v 1.0: initial release
v 1.1: major bugfixes
v 1.2: refactor to deal with spectral cubes

## Installation:
```
git clone https://github.com/mhardcastle/radioflux
cd radioflux/
sudo python setup.py install
```
Installing via `sudo pip install git+https://github.com/mhardcastle/radioflux` does not work - ds9crop.ds9 is not copied. Installing via `sudo pip install .` after cloning the repo will work, but the location of ds9crop.ds9 will not be printed.

Tested on Ubuntu 22 + DS9 8.3 and DS9 8.4b2
