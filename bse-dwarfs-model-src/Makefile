DEFAULT: bse
CMPLR = gfortran
#FFLAGS = -ffixed-line-length-none
FFLAGS = -fbounds-check -ffixed-line-length-none -g -ffpe-trap=invalid,zero,overflow
LFLAGS = const_bse.h zdata.h rand.h
 

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

SRCE2 = \
bse.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f lt2ubv.f mix.f mlwind.f mrenv.f orbel.f ubvtab.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f wd2ubv.f loadTab.f sort.f gt2mag.f imf.f gasdev.f getstar.f \
zdist.f galacticModel.f obstest.f printbpp.f RLOFanalysis.f getm2.f
 
OBJT2 = $(SRCE2:.f=.o)

SRCGR = \
grid.f comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f \
hrdiag.f instar.f kick.f lt2ubv.f mix.f mlwind.f mrenv.f orbel.f ubvtab.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f wd2ubv.f loadTab.f sort.f gt2mag.f imf.f gasdev.f getstar.f \
zdist.f galacticModel.f obstest.f getm2.f
 
OBJGR = $(SRCGR:.f=.o)


bse:    $(OBJT2) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT2) -o bse 

grid: $(OBJGR) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJGR) -o grid
	

OSG = galacticModel.o sampleGalaxy.o ran3.o
sampleGalaxy: $(OSG)
	$(CMPLR) $(FFLAGS) $(OSG) -o sampleGalaxy
