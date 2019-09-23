##########################################
###     LB2 Project: SS predictor      ###
##########################################

SHELL := /bin/bash

## Path definition
WD := /Users/davidebuzzao/Bioinformatica_Bologna/Second_Year/LB2/ModuleB/Project/

DATA := $(WD)Data/
DB := $(WD)DB/
SCRIPTS := $(WD)Scripts/
PLOTS := $(WD)Plots/

TS := $(DATA)TrainingSet/
BS := $(DATA)BlindingSet/

DSSP := $(BS)DSSP/
DSSP_RAW := $(DSSP)raw/
DSSP_CLEAN := $(DSSP)clean/
DSSP_SS := $(DSSP)ss/

FASTA1 := $(BS)FASTA/
FASTA2 := $(TS)FASTA/
PDB := $(BS)PDB/

PSIBLAST := $(TS)PSIBLAST/
PSSM :=$(PSIBLAST)/pssm
ALIGNMENTS := $(PSIBLAST)/alignments
PROFILESbin := $(PSIBLAST)/BINprofiles

## DataBase section
DB_UNIPROTfasta := $(DB)uniprot_sprot.fasta
DB_jPred4fasta := $(DB)jPred4.fasta
## Scripts section
SS_COMPO := $(SCRIPTS)SS_Composition.py
R_COMPO := $(SCRIPTS)R_Composition.py
TAXA_COMPO := $(SCRIPTS)Taxonomic_classification.py
SCOP_COMPO := $(SCRIPTS)SCOP_Composition.py
SLW17_COMPO := $(SCRIPTS)SLW17_Composition.py
SORT_CLUSTERS := $(SCRIPTS)Sort_cluster.py
SEQ_PROFILE := $(SCRIPTS)Extract_SeqProf.py
DSSP_INFO := $(SCRIPTS)get_DSSP.py

# User-dependent variables 
COV := 0
SEQ_ID := 30
E_VAL := 0.01
TABRES := $(BS)BlindSet.csv

all: SS_Composition.png R_Composition.png TAXA_cat.txt SCOP_Composition.png SLW17_compo.png TrainingSet BlindSet
TrainingSet: $(TS)OK.profile
BlindSet: $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.fasta $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.clean.list 

######################################
##	Training Set data preparation	##
######################################

$(TS)OK.profile: $(TS)jpred4.domain.list $(TS)OK.psiblast
	for domain in $$(cat $<); \
	do \
	python3 $(SEQ_PROFILE) <(grep -E '^[[:space:]]*[0-9]' $(PSSM)$$domain'.pssm') $$domain $(PROFILESbin); \
	done

$(TS)OK.psiblast: $(TS)jpred4.domain.list
	for id in $$(cat $<); \
	do \
	psiblast -query $(FASTA2)$$id.fasta -db $(DB_UNIPROTfasta) -evalue $(E_VAL) -num_iterations 3 \ 
	-out_ascii_pssm $(PSSM)$$id.pssm -num_descriptions 10000 -num_alignments 10000 -out $(ALIGNMENTS)$$id.alns.blast; \
	done

##############################
##	Blind Set construction	##
##############################
$(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.clean.list: $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.list $(BS)OK.ss
	for id in $$(ls $(DSSP_SS) | sed 's/.ss//g'); \
	do \
	grep $$id $<; \
	done > $@

$(BS)OK.ss: $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.list $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.fasta $(BS)OK.dssp
	for chain in $$(cat $<); \
	do \
	python3 $(DSSP_INFO) $$chain <(grep -A 1 $$chain $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.fasta) $(DSSP) NO; \
	done
	touch $@

$(BS)OK.dssp: $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.list $(BS)OK.pdb
	for chain in $$(cat $< | cut -d "_" -f 1); \
	do \
	mkdssp -i $(PDB)$$chain.pdb -o $(DSSP_RAW)$$chain.dssp; \
	done
	touch $@

$(BS)OK.pdb: $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.list
    cat $< | awk -F "_" '{print toupper($$1)}' | xargs -I '{}' -P 2 -n 1 \
        wget -nv -nc -P $(PDB) https://files.rcsb.org/download/'{}'.pdb
	touch $@

$(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.fasta: $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.list $(BS)BlindSet_c$(COV)_i$(SEQ_ID).clust.sort.fasta
	for id in $$(cat $<); do grep -A 1 $$id $(BS)BlindSet_c$(COV)_i$(SEQ_ID).clust.sort.fasta; done > $@

$(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred-30.list: $(BS)BlindSet_c$(COV)_i$(SEQ_ID).clust.sort.list $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred+30.list 
	comm -23 <(cat $< | sort) <(cat $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred+30.list) > $@

$(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred+30.list: $(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred.blast
	cat $< | awk -F "\t" '{if($$3>=30.0) print $$1}' | sort > $@

$(BS)BlindSet_c$(COV)_i$(SEQ_ID)TOjPred.blast: $(BS)BlindSet_c$(COV)_i$(SEQ_ID).clust.sort.fasta
	blastp -query $< -db $(DB)jPred4.fasta -out $@ -evalue $(E_VAL) -outfmt 6 -num_threads 2

$(DB)jPred4.fasta: $(DB)jPred4.fasta
	makeblastdb -in $< -dbtype prot

$(DB)jPred4.fasta: 
	for fasta in $$(ls $(FASTA1)); do cat $(FASTA1)$$fasta; done > $@


$(BS)BlindSet_c$(COV)_i$(SEQ_ID).clust.sort.list: $(BS)BlindSet_c$(COV)_i$(SEQ_ID).clust.sort.fasta
	grep '^>' $< | sed 's/>//g' > $@

$(BS)BlindSet_c$(COV)_i$(SEQ_ID).clust.sort.fasta: $(BS)BlindSet_c$(COV)_i$(SEQ_ID).sort $(BS)BlindSet.fasta
	for id in $$(cut -d ' ' -f 2 $< | awk -F ":" '{print $$1}'); do grep -A 1 $$id $(BS)BlindSet.fasta; done > $@

$(BS)BlindSet_c$(COV)_i$(SEQ_ID).sort: $(BS)BlindSet.res $(BS)BlindSet_c$(COV)_i$(SEQ_ID).clust
	python3 $(SORT_CLUSTERS) $^ > $@

$(BS)BlindSet_c$(COV)_i$(SEQ_ID).clust: $(BS)BlindSet.fasta
	blastclust -i $^ -o $@ -L 0._c$(COV) -S $(SEQ_ID)


$(BS)BlindSet.fasta: $(TABRES)
	tail -n +2 $< | head -n +5549 | sed 's/"//g' | awk -F "," '{print ">"$$1"_"$$2"\n"$$4}' > $@

$(BS)BlindSet.res: $(TABRES)
	cat $< | sed 's/"//g' | awk -F "," '{print $$1"_"$$2,$$3}' | tail -n +2 | head -n +5549 > $@


##############################
##	Statistical Analysis	##
##############################

## Sliding Window (size = 17) SS-contextualized residue composition
SLW17_compo.png: $(BS)FASTA_cat.fasta $(PLOTS)DSSP_cat.dssp
	python3 $(SLW17_COMPO) $< $(PLOTS)DSSP_cat.dssp 17 $(PLOTS)
	touch 'SLW17_compo.png'

## SCOP Classification
SCOP_Composition.png: $(PLOTS)jpred4.scop.txt legend.txt $(PLOTS)jpred4.scop.txt
	python3 $(SCOP_COMPO) $< legend.txt $(PLOTS)jpred4.scop.txt $(PLOTS)
	
$(PLOTS)jpred4.scop.txt: $(BS)jpred4.list.txt $(BS)dir.cla.scope.2.06-stable.txt
	for domain in $$(cat $<); do grep $$domain dir.cla.scope.2.06-stable.txt; done > $@
	
$(PLOTS)legend.txt: $(BS)dir.des.scope.2.06-stable.txt
	grep -v '^#' $< | awk -F "\t" '{print $$3, $$5}' > $@

## Taxonomic Classification
TAXA_cat.txt: jpred4.list.fasta
	grep '^>' $< | python3 $(TAXA_COMPO) $< | sort -t ':' -rgk2 > $@


## Total residue composition 
R_Composition.png: $(PLOTS)FASTA_cat.fasta
	python3 $(R_COMPO) $< $(PLOTS)

$(PLOTS)FASTA_cat.fasta: 
	for fasta in `ls $(FASTA1)`; do cat $(FASTA1)$$fasta | grep -v '^>'; done > $@


## SS composition Statistics
SS_Composition.png: $(PLOTS)DSSP_cat.dssp
	python3 $(SS_COMPO) $< $(PLOTS)

$(PLOTS)DSSP_cat.dssp: 
	for dssp in `ls $(DSSP)`; do cat $(DSSP)$$dssp | grep -v '^>'; done > $@

clean:
	rm $(PLOTS)DSSP_cat.dssp $(PLOTS)FASTA_cat.fasta
cleanall: 
	rm $(PLOTS)* 