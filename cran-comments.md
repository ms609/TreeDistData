## Test environments
* Windows 10 on local machine, R 3.6.0
* Windows 10 via check_win_devel(), R devel
* ubuntu 14.04.5 LTS (on travis-ci), R 3.6.0 and release
* Using check_rhub()

## R CMD check results
There were no ERRORs or WARNINGs.

There were two NOTEs:

* New submission. 
  - This package is a new submission.

* checking installed package size ... NOTE
    installed size is  5.4Mb
    sub-directories of 1Mb or more:
      data   5.3Mb
      
  - This package contains large datasets used to build the vignettes of
    the forthcoming package`TreeDist`.  I have compressed the datasets as far
    as possible.
    
    These datasets should very seldom change, and are primarily used to build 
    vignettes (though there are some other use cases);
    as such, I feel that it is most appropriate to serve them as a standalone
    package.


## Downstream dependencies
There are currently no downstream dependencies for this package.
