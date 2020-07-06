## Test environments
* Windows 10 on local machine, R 4.0.2
* Windows 10 via check_win_devel(), R devel
* ubuntu 14.04.5 LTS (on travis-ci), R 3.6.0 and release
* Using check_for_cran()

## R CMD check results
There were no ERRORs or WARNINGs.

There were two NOTEs:

> New submission. 

This package is a new submission.


> Possibly mis-spelled words in DESCRIPTION:
>   unrooted (10:73)

False positive: the spelling is correct.

> Availability using Additional_repositories specification:
>   ?   ?   https://ms609.github.io/packages



> Size of tarball: 18827486 bytes
      
This package contains large datasets used to build the vignettes of
the package 'TreeDist'.  I have compressed the datasets as far
as possible.
  
These datasets should very seldom change, and are primarily used to build 
vignettes (though there are some other use cases);
as such, I feel that it is most appropriate to serve them as a standalone
package.


## Downstream dependencies
There are currently no downstream dependencies for this package.
