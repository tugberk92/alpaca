FC = gfortran

#####################

HOME = $(PWD)
SOURCEDIR = $(PWD)/src
INCPATH = $(SOURCEDIR)/inc

OBJ_PATH = $(PWD)/obj/

FFLAGS 	= -fno-automatic -fno-f2c -O2 -g  -I$(INCPATH)

DIRS	 =	$(SOURCEDIR)/int:\
		$(SOURCEDIR)/main:\
		$(SOURCEDIR)/phase:\
		$(SOURCEDIR)/subamps:\
		$(SOURCEDIR)/cscalc:\
		$(SOURCEDIR)/unw:\
		$(SOURCEDIR)/user:\
		$(SOURCEDIR)/var:\
		$(SOURCEDIR)/ion:\

VPATH = $(DIRS)

#############

Ionf = \
ionpars.o \
rhonorm.o \
rho.o \
rhoxycalc.o \
rhoxy.o \
tp.o \
tpcalc.o \
rhoxyint.o \
tpint.o \

Intf = \
vegas.o \
rann.o \

Phasef = \
2bodyn.o \
boost.o \
rambo.o \
axionps.o \
gamdecay.o \

Subampsf = \
axion.o \

Cscalcf = \
formfacgam.o \
cscalc.o \

Userf = \
cuts.o \
histo.o \

Mainf = \
headerout.o \
main.o \
process.o \

Unwf = \
unweight.o \
unwprint.o \
headerlhe.o \

Varf = \
string.o \
varfuncs.o \

#

sCODEi = $(Mainf) $(Intf) $(Phasef) $(Subampsf) $(Cscalcf) $(Userf) $(Unwf) $(Varf) $(Ionf)

sCODE = $(patsubst %,$(OBJ_PATH)%,$(sCODEi))

###########

all : alpaca

alpaca.o:	alpaca.f
	$(FC) $(FFLAGS) -I$(INCPATH) -c  $< -o $@

$(OBJ_PATH)%.o: %.f
	$(FC) $(FFLAGS) -I$(INCPATH) -c  $< -o $@

alpaca : $(OBJ_PATH)alpaca.o $(sCODE) 
	$(FC) $^ -o bin/$@

clean:
	rm -f  bin/alpaca *.o $(OBJ_PATH)*.o
