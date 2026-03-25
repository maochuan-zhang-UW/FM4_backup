c nx0 = maximum source-station distance bins for look-up tables
c nd0 = maximum source depth bins for look-up tables
c nindex = maximum number of velocity model look-up tables
c dep1,dep2,dep3 = minimum source depth, maximum, and interval
c del1,del2,del3 = minimum source-station distance, maximum, and interval
c pmin = minimum ray parameter for ray tracing
c nump = number of rays traced
*      parameter(nx0=101,nd0=14,nindex=10)
*      parameter(dep1=0,dep2=39,dep3=3)
*      parameter(del1=0,del2=200,del3=2)
*      parameter(pmin=0,nump=9000)

      parameter(nx0=251,nd0=41,nindex=10)
      parameter(dep1=0.0,dep2=6.0,dep3=0.2)
      parameter(del1=0.0,del2=50.0,del3=0.2)
*      parameter(pmin=0,nump=90000)

*      parameter(nx0=101,nd0=14,nindex=10)
*      parameter(dep1=0,dep2=26,dep3=2)
*      parameter(del1=0,del2=200,del3=2)
      parameter(pmin=0,nump=9000)
