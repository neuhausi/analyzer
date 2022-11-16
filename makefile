# $Id: makefile,v 1.36 2009/06/29 16:14:52 neuhausi Exp $

## unix

CC 	= gcc

#CFLAGS  = -O
#CFLAGS  =-Wall -g -pg # all warnings with profiler and debuging
#CFLAGS  =-Wall -g # all warnings
CFLAGS  = -Wall

CP	= cp

CURDIR  = .

DATE	= $(shell date +%m_%d_%y)

DESTDIR = bin

CAT	= cat

DIFF	= diff

ECHO	= echo

EDITOR  = xemacs

GZIP	= /bin/gzip -fq9

LIBS    = -lm

LINK    = $(CC) $(CFLAGS) -o $(CURDIR)/$@

MAKE    = make

MKDIR	= mkdir

MV      = mv

RM      = rm -rf

TAR	= tar -cvf

TARF	= analyzer_$(DATE).tar

TART	= analyzer_$(DATE)_t.tar

TEST    = test

## Files

MAINS   = anova.c correct.c fit.c glm.c hyper.c ks.c mixed.c nipals.c nna.c pca.c svd.c ttest.c zzz.c

PROGS   = anova correct fit glm hyper ks mixed mixedtemp nipals nna pca svd ttest zzz

PROGST  = anovat correctt fitt glmt hypert kst mixedt nipalst nnat pcat svdt ttestt

UTILS   = files.c function.c matrix.c methods.c metrics.c solution.c stats.c utils.c var_metrics.c

HEADERS = anova.h correct.h fit.h glm.h hyper.h ks.h mixed.h mixedtemp.h nipals.h nna.h pca.h svd.h ttest.h files.h function.h matrix.h methods.h metrics.h solution.h stats.h utils.h var_metrics.h

## Dependencies

ANOVA   = utils.o files.o stats.o anova.o

CORRECT = utils.o files.o stats.o methods.o correct.o

FIT     = utils.o files.o stats.o solution.o function.o fit.o

GLM     = utils.o files.o stats.o solution.o matrix.o glm.o

HYPER   = utils.o files.o stats.o hyper.o

KS      = utils.o files.o stats.o ks.o

MIXED   = utils.o files.o stats.o solution.o matrix.o mixed.o

MIXEDTEMP   = utils.o files.o stats.o solution.o matrix.o mixedtemp.o

NIPALS  = utils.o files.o stats.o solution.o matrix.o methods.o nipals.o

NNA     = utils.o files.o stats.o metrics.o var_metrics.o methods.o nna.o

PCA     = utils.o files.o stats.o solution.o methods.o pca.o

SVD     = utils.o files.o stats.o solution.o methods.o svd.o

TTEST   = utils.o files.o stats.o ttest.o

ZZZ	= utils.o files.o methods.o solution.o function.o matrix.o metrics.o var_metrics.o stats.o zzz.o

## Targets

all:
	@$(ECHO) ""	
	@$(ECHO) "Making all the ANALYZER programs"
	@$(ECHO) "Using compiler $(CC) with options $(CFLAGS)"
	@$(ECHO) "Making in directory $(CURDIR)"
	@$(ECHO) ""
	@$(MAKE) -i $(PROGS)
	@$(ECHO) ""	
	@$(ECHO) "All programs made OK!"
	@$(ECHO) "Comments to isaac.neuhaus@bms.com"	
	@$(ECHO) ""	

archive:
	$(call createdir, analyzer)
	$(CP) $(MAINS) $(UTILS) $(HEADERS) makefile README analyzer
	$(TAR) $(TARF) analyzer
	#$(GZIP) $(TARF)
	@$(RM) analyzer
	$(call createdir, analyzer)
	$(CP) -R t analyzer
	$(RM) analyzer/t/CVS
	$(TAR) $(TART) analyzer
	#$(GZIP) $(TART)
	@$(RM) analyzer

clean:
	$(RM) *.o core $(PROGS)

edit:
	$(EDITOR) $(MAINS) $(UTILS) $(HEADERS) makefile README &

help:
	$(foreach PRG, $(PROGS), ./$(PRG) -v;)

install: 
	$(call createdir, $(DESTDIR))
	$(CP) $(PROGS) $(DESTDIR)

test:
	@$(ECHO) "Runing tests"
	@$(ECHO) ""
	@$(MAKE) -i $(PROGST)
	@$(ECHO) ""
	@$(ECHO) "Comments to isaac.neuhaus@bms.com"	
	@$(ECHO) ""	

## Programs

anova: $(ANOVA)
	$(LINK) $(ANOVA) $(LIBS)

anovat: 
	@$(ECHO) ""
	@$(ECHO) "Tests for anova"
	@$(ECHO) "  1. one-way (balanced)"
	@$(CURDIR)/anova -d ./t/1f.3l.b.dat -f ./t/1f.3l.b.fac > ./t/anova.out
	@$(ECHO) "  2. one-way (unbalanced)"
	@$(CURDIR)/anova -d ./t/1f.3l.u.dat -f ./t/1f.3l.u.fac >> ./t/anova.out
	@$(ECHO) "  3. one-way (missing data)"
	@$(CURDIR)/anova -d ./t/1f.3l.md.dat -f ./t/1f.3l.b.fac >> ./t/anova.out
	@$(ECHO) "  4. two-way between interaction (balanced)"
	@$(CURDIR)/anova -d ./t/2f.3l.2l.b.dat -f ./t/2f.3l.2l.b.fac >> ./t/anova.out
	@$(ECHO) "  5. two-way between interaction (unbalanced)"
	@$(CURDIR)/anova -d ./t/2f.3l.2l.u.dat -f ./t/2f.3l.2l.u.fac >> ./t/anova.out
	@$(ECHO) "  6. two-way between interaction (missing data)"
	@$(CURDIR)/anova -d ./t/2f.3l.2l.md.dat -f ./t/2f.3l.2l.b.fac >> ./t/anova.out
	@$(ECHO) "  7. two-way between no interaction (balanced)"
	@$(CURDIR)/anova -d ./t/2f.3l.2l.b.dat -f ./t/2f.3l.2l.b.fac -n >> ./t/anova.out
	@$(ECHO) "  8. two-way between no interaction (unbalanced)"
	@$(CURDIR)/anova -d ./t/2f.3l.2l.u.dat -f ./t/2f.3l.2l.u.fac -n >> ./t/anova.out
	@$(ECHO) "  9. two-way between no interaction (missing data)"
	@$(CURDIR)/anova -d ./t/2f.3l.2l.md.dat -f ./t/2f.3l.2l.b.fac -n >> ./t/anova.out
	@$(ECHO) " 10. three-way between interaction (balanced)"
	@$(CURDIR)/anova -d ./t/3f.2l.2l.2l.b.dat -f ./t/3f.2l.2l.2l.b.fac >> ./t/anova.out
	@$(ECHO) " 11. three-way between interaction (unbalanced)"
	@$(CURDIR)/anova -d ./t/3f.2l.2l.2l.u.dat -f ./t/3f.2l.2l.2l.u.fac >> ./t/anova.out
	@$(ECHO) " 12. three-way between interaction (missing data)"
	@$(CURDIR)/anova -d ./t/3f.2l.2l.2l.md.dat -f ./t/3f.2l.2l.2l.b.fac >> ./t/anova.out
	@$(ECHO) " 13. three-way between no interaction (balanced)"
	@$(CURDIR)/anova -d ./t/3f.2l.2l.2l.b.dat -f ./t/3f.2l.2l.2l.b.fac -n >> ./t/anova.out
	@$(ECHO) " 14. three-way between no interaction (unbalanced)"
	@$(CURDIR)/anova -d ./t/3f.2l.2l.2l.u.dat -f ./t/3f.2l.2l.2l.u.fac -n >> ./t/anova.out
	@$(ECHO) " 15. three-way between no interaction (missing data)"
	@$(CURDIR)/anova -d ./t/3f.2l.2l.2l.md.dat -f ./t/3f.2l.2l.2l.b.fac -n >> ./t/anova.out
	@$(ECHO) " 16. one-way within subject"
	@$(CURDIR)/anova -d ./t/2f.6l.2l.w.dat -f ./t/2f.6l.2l.w.fac -t within >> ./t/anova.out
	@$(ECHO) " 17. two-way within subject"
	@$(CURDIR)/anova -d ./t/3f.4l.2l.4l.w.dat -f ./t/3f.4l.2l.4l.w.fac -t within >> ./t/anova.out
	@$(ECHO) " 18. two-way mixed within subject"
	@$(CURDIR)/anova -d ./t/3f.4l.3l.4l.m.dat -f ./t/3f.4l.3l.4l.m.fac -t mixed >> ./t/anova.out
	$(call filecheck, anova, ./t/anova.out, ./t/anova.reg)

correct: $(CORRECT)
	$(LINK) $(CORRECT) $(LIBS)

correctt: 
	@$(ECHO) ""
	@$(ECHO) "Tests for correct"	
	@$(ECHO) "  1. FDR"
	@$(CURDIR)/correct -d ./t/correct.dat -i 2 > ./t/correct.out
	$(call filecheck, correct, ./t/correct.out, ./t/correct.reg)

fit: $(FIT)
	$(LINK) $(FIT) $(LIBS)

fitt: 
	@$(ECHO) ""
	@$(ECHO) "Tests for fit"	
	@$(ECHO) "  1. No fixed parameters"
	@$(CURDIR)/fit -d t/fit.dat -f t/fit.fac -c t/fit.cov -m Dose -a Treatment -t nlls4p -z C> ./t/fit.out
	$(call filecheck, fit, ./t/fit.out, ./t/fit.reg)

glm: $(GLM)
	$(LINK) $(GLM) $(LIBS)

glmt: 
	@$(ECHO) ""
	@$(ECHO) "Tests for glm"	
	@$(ECHO) "  1. Fixed model"
	@$(CURDIR)/glm -d ./t/3f.2l.2l.2l.b.dat -f ./t/3f.2l.2l.2l.b.fac -m Grp1 -m Grp2 -m Grp3 > ./t/glm.out
	@$(ECHO) "  2. Mixed model"
	@$(CURDIR)/glm -d ./t/2f.6l.2l.w.dat -f ./t/2f.6l.2l.w.fac -m Grp -m Sub >> ./t/glm.out
	@$(ECHO) "  3. Missing treatment combination (empty cells)"
	@$(CURDIR)/glm -d t/emptycells.dat -f t/emptycells.fac -i 'R*T' >> ./t/glm.out
	@$(ECHO) "  4. Missing treatment combination (nested)"
	@$(CURDIR)/glm -d t/emptycells.dat -f t/emptycells.fac -i 'T*N' >> ./t/glm.out
	@$(ECHO) "  5. Missing covariates"
	@$(CURDIR)/glm -d t/glm.dat -c t/glm.cov -f t/glm.fac -m CD -m FA -v -n Var1 >> ./t/glm.out
	$(call filecheck, glm, ./t/glm.out, ./t/glm.reg)

hyper: $(HYPER)
	$(LINK) $(HYPER) $(LIBS)

hypert: 
	@$(ECHO) ""
	@$(ECHO) "Tests for hyper"	
	@$(ECHO) "  1. Small data set"
	@$(CURDIR)/hyper -d ./t/hyper.small.dat -l ./t/hyper.lis > ./t/hyper.out
	@$(ECHO) "  1. Big data set"
	@$(CURDIR)/hyper -d ./t/hyper.big.dat -n GeneList2 >> ./t/hyper.out
	$(call filecheck, hyper, ./t/hyper.out, ./t/hyper.reg)

ks: $(KS)
	$(LINK) $(KS) $(LIBS)

kst: 
	@$(ECHO) ""
	@$(ECHO) "Tests for ks"	
	@$(ECHO) "  1. NR example"
	@$(CURDIR)/ks -d ./t/ks.dat -r A0 > ./t/ks.out
	$(call filecheck, ks, ./t/ks.out, ./t/ks.reg)

mixed: $(MIXED)
	$(LINK) $(MIXED) $(LIBS)

mixedt: 
	@$(ECHO) ""
	@$(ECHO) "Tests for mixed"	
	@$(ECHO) "  1. Mixed model"
	@$(CURDIR)/mixed -d ./t/3f.2l.2l.2l.b.dat -f ./t/3f.2l.2l.2l.b.fac -m Grp1 -m Grp2 -m Grp3 > ./t/mixed.out
	$(call filecheck, mixed, ./t/mixed.out, ./t/mixed.reg)

mixedtemp: $(MIXEDTEMP)
	   $(LINK) $(MIXEDTEMP) $(LIBS)

nipals: $(NIPALS)
	$(LINK) $(NIPALS) $(LIBS)

nipalst: 
	@$(ECHO) ""
	@$(ECHO) "Tests for nipals"	
	@$(ECHO) "  1. No misisng data"
	@$(CURDIR)/nipals -d ./t/nipals.dat -a > ./t/nipals.out
	@$(ECHO) "  2. Misisng data"
	@$(CURDIR)/nipals -d ./t/nipals.md.dat -a >> ./t/nipals.out
	$(call filecheck, nipals, ./t/nipals.out, ./t/nipals.reg)

nna: $(NNA)
	$(LINK) $(NNA) $(LIBS)

nnat: 
	@$(ECHO) ""
	@$(ECHO) "Tests for nna"	
	@$(ECHO) "  1. Corregulation (euclidean)"
	@$(CURDIR)/nna -d ./t/nna.dat -n M55593_at > ./t/nna.out
	$(call filecheck, nna, ./t/nna.out, ./t/nna.reg)

pca: $(PCA)
	$(LINK) $(PCA) $(LIBS)

pcat: 
	@$(ECHO) ""
	@$(ECHO) "Tests for pca"	
	@$(ECHO) "  1. Samples and genes"
	@$(CURDIR)/pca -d ./t/pca.dat > ./t/pca.out
	$(call filecheck, pca, ./t/pca.out, ./t/pca.reg)

svd: $(SVD)
	$(LINK) $(SVD) $(LIBS)

svdt: 
	@$(ECHO) ""
	@$(ECHO) "Tests for svd"	
	@$(ECHO) "  1. Correlation"
	@$(CURDIR)/svd -d ./t/svd.dat > ./t/svd.out
	$(call filecheck, svd, ./t/svd.out, ./t/svd.reg)

ttest: $(TTEST)
	$(LINK) $(TTEST) $(LIBS)

ttestt: 
	@$(ECHO) ""
	@$(ECHO) "Tests for ttest"	
	@$(ECHO) "  1. t-test (equal variance)"
	@$(CURDIR)/ttest -d ./t/ttest.dat -f ./t/ttest.fac -t equal > ./t/ttest.out
	@$(ECHO) "  2. t-test (unequal variance)"
	@$(CURDIR)/ttest -d ./t/ttest.dat -f ./t/ttest.fac -t unequal >> ./t/ttest.out
	@$(ECHO) "  3. t-test (paired)"
	@$(CURDIR)/ttest -d ./t/ttest.dat -f ./t/ttest.fac -t paired >> ./t/ttest.out
	$(call filecheck, ttest, ./t/ttest.out, ./t/ttest.reg)

zzz:	$(ZZZ)
	$(LINK) $(ZZZ) $(LIBS)

define filecheck
  @if ($(DIFF) $(2) $(3)); \
    then $(ECHO) "Passed"; $(ECHO) ""; \
    else $(ECHO) "Tests for $(1) failed"; $(ECHO) ""; \
  fi
endef

define createdir
  @if (! $(TEST) -e $(1)); \
    then $(ECHO) "creating dir $(1)"; $(MKDIR) $(1); \
  fi
endef
