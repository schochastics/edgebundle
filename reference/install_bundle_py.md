# install python dependencies for hammer bundling

install datashader and scikit-image

## Usage

``` r
install_bundle_py(method = "auto", conda = "auto")
```

## Arguments

- method:

  Installation method (by default, "auto" automatically finds a method
  that will work in the local environment, but note that the
  "virtualenv" method is not available on Windows)

- conda:

  Path to conda executable (or "auto" to find conda using the PATH and
  other conventional install locations)
