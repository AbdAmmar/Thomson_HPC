
<br />
<div align="center">
  <a href="https://github.com/ALRAKIK/Thomson">
    <img src="src/logo.png" alt="Logo" width="140" height="140">
  </a>

  <h3 align="center">Thomson Project</h3>

  <p align="center">
    Solve Thomson problem in Clifford tori for one, two and three dimension
    <br />
    <a href="https://github.com/ALRAKIK/Thomson"><strong>Explore the Project »</strong></a>
    <br />
    <br />
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#thomson">About The Project</a>
    </li>
    <li>
      <a href="#prerequisite">Getting Started</a>
      <ul>
        <li><a href="#prerequisite">Prerequisites</a></li>
        <li><a href="#build">Build</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#example">Example</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

# Thomson 
The Thomson problem is a classic problem in physics and mathematics that deals with finding the minimum energy configuration of a collection of point charges (electrons) on the surface of a sphere. The problem was first proposed by the physicist J.J. Thomson in 1904, and it has been studied extensively since then.

The Thomson problem has many applications in physics and chemistry, including the study of atomic and molecular structures. The problem is also of interest in mathematics and computer science, as it involves finding the optimal arrangement of a large number of points in three-dimensional space.

Thomson program has been developed to solve the Thomson problem inside the torus. This program uses Conjucated gradient mathematical technique and computational algorithms to find the optimal configuration of charges on the surface of a torus, which is a doughnut-shaped object.

The program written using FORTRAN 90 language with Python3 GUI to help the user making the input files and running the program.

We define the euclidean distance in the torus by the expression: 

```math
d_{i j}=\sqrt{\sum_{\alpha=1}^n\left[\frac{L_\alpha}{\pi} \sin \left(\frac{\pi}{L_\alpha}\left|x_{\alpha, i}-x_{\alpha, j}\right|\right)\right]^2}
```

![Screenshot](src/GUI.png)

# Prerequisite

* 1 - Fortran compiler (gfortran 9.4 recommended):

```
sudo apt install gfortran
```
* 2 - Fortran library (llapack, openmp):
  
```
sudo apt-get install liblapack-dev libopenblas-dev libomp-dev
```
* 3 - For animation (gnuplot, FFmpeg):

```
sudo apt-get install gnuplot ffmpeg
```
* 4 - Python (Python 3, not working using python 2 because it need some library )

```
sudo apt-get install python3 python3-pip
```
* 5 - Python library (numpy,tkinter,customtkinter) to use the python GUI

```
sudo apt-get install python3-tk
pip3 install customtkinter numpy 
```

 
# Build

* build the program
  
```sh
make Thomson
```
* clean the bin files
```sh
make clean
```
# Usage 

* Using Python code

```sh
python Thomson.py
``` 
* Thomson

```sh
./Thomson "The input file name"
```

* Input file

first line to identify the number of dimension  1,2 or 3.

second line to identify the number of electrons         .

third line for the tolerance                            .

fourth line for the maximum number of iteration         .

then the user can add any keywords for the list (with any order) :

| keywords                      | Description | 
| :---:                         |     :---       |
| `box: num num num`   | change the size of the box on x,y,z  "default $2\pi,2\pi,2\pi$" |
| `random`                      |  start from random geometry (every axis have a random number between  [0,1]  |
| `multiply`                    | multiply the input geometry by the size of the box and distrubute the electron all over the box |
| `show`                        | show all the result (the geometry in every step of iteration, energy of the system , and the norm of the gradient) |
| `density`                     | choose a fixed density (one electron per the unit of length) and ignore the size of the box written before or after  |
| `density rectangle`           | choose a fixed density (one electron per the unit of length) in case of 2D if user want to choose the ratio between the first and the second axis as $\sqrt(3)/2$  |
| `hessian`                     | show the hessian matrix at the convergance |
| `distance`                    | show the euclidean and the geodesic matrix at the convergance |
| `animation`                   | make a video for the optimizations steps or if the user are at a fixed point (minimum or saddle point) gives a picture  |
| `animation origin`            | same as animation but force one electron to be at the origin of the torus [0,0,0] |

geometry (the last keywork) , user have to define the input geometry below it until he reach the last number of electron like this: 

```sh
geomrtry
1   number number number 
2   number number number
3   number number number
.
.
.
.
```


# Example

there are some example of the inputfile in the [examples folder](https://github.com/ALRAKIK/Thomson/tree/main/example) you can try to use them.
  
# Contact

Stefano Evangelisti - stefano.lcpq@gmail.com

Arjan Berger - arjan.berger@irsamc.ups-tlse.fr

Amer Alrakik - alrakikamer@gmail.com

Project Link: [https://github.com/ALRAKIK/Thomson](https://github.com/ALRAKIK/Thomson)

# License

Distributed under the MIT License. See `LICENSE.txt` for more information.

# Acknowledgments

Thanks to ... 
