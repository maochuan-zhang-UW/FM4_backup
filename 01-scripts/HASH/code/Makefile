OBJS1=	fmech_subs.o \
        pol_subs.o \
        station_subs.o \
        uncert_subs.o \
	util_subs.o

OBJS2=	fmech_subs.o \
        pol_subs.o \
        station_subs.o \
        uncert_subs.o \
	util_subs.o \
        vel_subs.o

OBJS3=  pol_subs.o \
        station_subs.o \
        uncert_subs.o \
	util_subs.o \
        vel_subs.o \
        fmamp_subs_ww.o

OBJS3b=  pol_subs.o \
        station_subs.o \
        uncert_subs.o \
	util_subs.o \
        vel_subs.o \
        fmech_subs_ww.o \
        fmamp_subs.o

OBJS4=  pol_subs.o \
        station_subs.o \
        uncert_subs.o \
	util_subs.o \
        vel_subs.o \
        fmech_subs_10.o \
        fmamp_subs.o

hash_driver1: Makefile hash_driver1.f $(OBJS1)
	gfortran -std=legacy  hash_driver1.f $(OBJS1) -O -o hash_driver1

hash_driver2: Makefile hash_driver2.f $(OBJS2)
	gfortran -std=legacy  hash_driver2.f $(OBJS2) -O -o hash_driver2

hash_driver3_only_pol: Makefile hash_driver3_only_pol.f $(OBJS3b)
	gfortran -std=legacy  hash_driver3.f $(OBJS3b) -O -o hash_driver3_only_pol

hash_driver3: Makefile hash_driver3.f $(OBJS3)
	gfortran -std=legacy  hash_driver3.f $(OBJS3) -O -o hash_driver3 -fbounds-check

hash_driver3_10: Makefile hash_driver3.f $(OBJS4)
	gfortran -std=legacy  hash_driver3_10.f $(OBJS4) -O -o hash_driver3_10 -fbounds-check

pol_subs.o : Makefile pol_subs.f
	gfortran -std=legacy pol_subs.f -O -c -fbounds-check

fmech_subs.o : Makefile fmech_subs.f
	gfortran -std=legacy fmech_subs.f -O -c -fbounds-check

fmech_subs_ww.o : Makefile fmech_subs_ww.f
	gfortran -std=legacy fmech_subs_ww.f -O -c -fbounds-check

fmech_subs_10.o : Makefile fmech_subs_10.f
	gfortran -std=legacy fmech_subs_10.f -O -c -fbounds-check

fmamp_subs.o : Makefile fmamp_subs.f
	gfortran -std=legacy fmamp_subs.f -O -c -fbounds-check

fmamp_subs_ww.o : Makefile fmamp_subs_ww.f
	gfortran -std=legacy fmamp_subs_ww.f -O -c -fbounds-check

station_subs.o : Makefile station_subs.f
	gfortran -std=legacy station_subs.f -O -c 

uncert_subs.o : Makefile uncert_subs.f
	gfortran -std=legacy uncert_subs.f -O -c -fbounds-check

util_subs.o : Makefile util_subs.f
	gfortran -std=legacy util_subs.f -O -c -fbounds-check

vel_subs.o : Makefile vel_subs.f
	gfortran -std=legacy vel_subs.f -O -c -fbounds-check



