
## Overview

**Soft Cell Simulator** is an interactive web application developed in **R** using the **Shiny** framework. The application visualizes the iterative algorithm for transforming an admissible convex polygon into a planar soft cell.

The implementation demonstrates how repeated corner cutting and boundary correction gradually smooth the polygon boundary while preserving the distinguished vertices (P^+) and (P^-). Users can explore the effect of the algorithm parameter (\rho), inspect intermediate iterations, and visualize the geometric construction step by step.

## Live Application

The application is publicly available at:

https://019e31a3-c969-cbad-46ce-10a667550a85.share.connect.posit.cloud/

## Features

* Interactive simulation of the soft-cell construction algorithm.
* Adjustable cutting parameter (\rho).
* Animation of successive iterations.
* Visualization of the current polygonal approximation.
* Optional display of:

  * the initial polygon,
  * interior chain vertices,
  * newly generated cut points,
  * boundary-corrected cut points.
* Highlighting of the distinguished vertices (P^+) and (P^-).
* Real-time rendering using **ggplot2**.

## Mathematical Background

The implemented algorithm operates on two boundary chains connecting the distinguished vertices (P^+) and (P^-).

At each iteration:

1. Every interior vertex is replaced by two new points obtained through a corner-cutting rule.
2. Near the distinguished vertices, a boundary correction procedure is applied.
3. The correction uses orthogonal projection onto the vertical lines passing through (P^+) and (P^-).
4. The resulting chains define a refined polygonal approximation.

As the number of iterations increases, the generated polygon converges toward a smooth planar soft cell whose boundary arcs meet tangentially at (P^+) and (P^-).

## Controls

### Parameter (\rho)

Controls the intensity of the local corner-cutting operation.

* Smaller values produce gentler modifications.
* Larger values produce more aggressive smoothing.

### Number of Iterations (k)

Determines how many iterations of the algorithm are applied.

The slider supports animation, allowing users to observe the evolution of the polygonal sequence.

### Visualization Options

* **Show initial polygon** – displays the original polygon as a dashed outline.
* **Show chain vertices** – displays the interior vertices of the current chains.
* **Show cut points** – displays the newly generated cut points and highlights boundary-corrected points.

## Technology Stack

* R
* Shiny
* bslib
* ggplot2

## Running Locally

### Install Dependencies

```r
install.packages(c(
  "shiny",
  "bslib",
  "ggplot2"
))
```

### Launch Application

```r
shiny::runApp("app.R")
```

or open `app.R` in RStudio and click **Run App**.

## Repository Structure

```text
.
├── app.R
├── README.md
```

## Author

Eleonora Detić
Algebra Bernays University, Zagreb

