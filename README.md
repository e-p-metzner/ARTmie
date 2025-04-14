# ARTmie

---

`ARTmie` is a python library to calculate optical properties of spherical particles.<br>
Those can be simple, single coated or drawn from a particle size distribution.<br>

`ARTmie` uses a lot of improved C++ in the backend for accelaration of various calculations.<br>

The Code was mainly ported from Matlab and Fortran to C++.



## Installation

<h4>from the Python Package Index (PyPI)</h4>

```python
pip install ARTmie
```

<h4>from source</h4>

1) Download the source code from this git by clicking the triangular button next to `Code` and selecting e.g. [zip](https://gitlab.kit.edu/enrico.metzner/artmie/-/archive/main/artmie-main.zip) underneath *Download source code*.<br>
2) Extract to source code to a desired location.<br>
3) Open a terminal with an appropriate python environment (Python 3.0 or higher) and Numpy (1.5 or higher, tested with 2.0).<br>
4) Change to the directory of the ARTmie project (dir above src):<br>
```
cd /path/to/ARTmie
```
5) Install it with the following command:<br>
```python
python -m pip install .
```

<h4>Requirements</h4>

- Python 3.0 or higher
- Numpy 1.5 or higher, preferable 2.0 or higher
- On windows with VS-Code: Microsoft Visual C++ 14.0 or higher



## License

ARTmie is available under the 3-Clause BSD license. See [LICENSE](./LICENSE) for license information.



## Contributor(s)

Enrico P. Metzner



## Roadmap

for version 0.1.1:
 - write documentation
 - create test-script to validate the installation

for version 0.2.0:
 - 1dimensional arrays for input arguments sizepar1 & sizepar2 in Size\_Distribution\_Optics and Size\_Distribution\_Phase\_Function
 - Size\_Distribution\_Phase\_Function can also take an array of angles (theta) instead of predefined equally spaced 721 angles from 0° to 180°
