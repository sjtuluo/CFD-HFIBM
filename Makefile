# CIP LINUX:
  FCOMPILER1 = ifort #gfortran ifort #g77
  FCOMPILER2 = gfortran #gfortran ifort #g77 
#  FFALGS  = -g -shared-intel -mcmodel=large -openmp  #-p 
  FFLAGS1  = -qopenmp -O2 -fp-model source -ipo -extend-source 132 -shared-intel -mcmodel=large #-g -shared-intel -check bounds #-axSSE4.1   #-shared-intel -mcmodel=medium  #-p 
  FFLAGS2  = -O2 -fp-model source -ipo -extend-source 132  #-fopenmp #-O3  #-shared-intel -mcmodel=medium  #-p 
  FLINK  = ${FFLAGS1} #-lgfortran #-static #-lgfortran
#  DIR1 = /home/zjx/water_entry/new
  INSTALLDIR = #/home/zy/omp/ocerm/
#

  EXOBJS = ADVTK.o ADVTKD.o ADVU.o ADVV.o ADVVIS.o ADVW.o ARCHIVE.o ATRDE.o \
           BCDATA.o BCITERATION.o BCOND.o BRINV.o DESSA.o DESSST.o DYN.o ELTION.o FIRST.o GETCOR.o\
           GRAD.o IBM.o IBMALDF.o IBMALGC.o IBMALIDC.o IBMALIDP.o IBMARCHIVE.o IBMCLBD.o\
		   IBMCLBP.o IBMINIT.o IBMRGD.o IBMVISUAL.o INIVOR.o OCERM.o PORECAL.o PROFTK.o PROFTKD.o \
		   PROFV.o PROFVIS.o PROFW.o REUV.o SAVEDIVERGE.o SETDOM.o SGSMODEL.o SOLVE3DPOLCG.o SOLVEELFPOLCG.o \
		   STATISTICS.o SUBGRIDH.o SUBGRIDV.o TURGEN.o TVDSCHEMEH.o TVDSCHEMEV.o UPDATEFLOW.o\
           UVFN.o VERTVEL.o WAVEBREAKING.o WAVEGEN.o WENO.o WREAL.o ZEROS.o
#
#  MYINCS =


ocerm:$(EXOBJS)
	$(FCOMPILER1) $(FLINK) $(EXOBJS) -o $(INSTALLDIR)HydroFlow

# SUFFIXES-DEFINITION:
%.o : %.F ./Include/OCERM_INF ./Include/VORGEN_INF
	$(FCOMPILER1) -c   $(FFLAGS1)  $<

update:
	find . -name "*.F" -exec touch {} \; 

clean:
	rm *.o  HydroFlow
         
#
