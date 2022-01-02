# Efficent generation of Mandelbrot set imagesequences using MPI and OpenMP

## Members

|     Name     |  Surname  |    MAT     |
| :----------: | :-------: | :--------: |
|    Luigi     | Riz | **223823** |
| Raffaele |  Pojer   | **224012** |

## Compiling

```
$ ./compile.sh mandelbrot.c
```

## Run

```
$ # manage the execution parameters into mandelbrot.sh
$ qsub mandelbrot.sh
```

## Arguments for the c application

- -?, --help

  Show helpful information

- -f, --zoom-factor=zoomFactor

  Use zoomFactor as magnification factor

- -i, --iterations=ITER

  Use ITER as IterationMax

- -r, --resolution=R

  Produce RxR images

- -t, --threads=threadCount

  Use threadCount threads
  
- -x, --final-x=X

  Specify horizontal coordinate for the zoom
  
- -y, --final-y=Y

  Specify vertical coordinate for the zoom
  
- -z, --zooms=zoomCount

  Produce zoomCount images in output
  
