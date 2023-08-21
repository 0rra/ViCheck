# ViCheck

Scripts ran in following order:

(results processing and annotation scripts)

1. ProcessBLASTHits.py
2. GetTax.py
3. MergeBLASTHits.py
4. GetSubId.py
5. SubjectHitAnno.py
   
(plotting and analysis scripts)

6. SummariseHead.py
7. SummariseHitFiles.py
8. SummariseSubjectFiltes.py
9. FinalFiltering.py
    
(tool to flag problematic sequences in other BLAST datasets)

10. CheckMyBLAST.py

# examples of how ProcessBLASTHits.py was run
python ProcessBLASTHits.py refseq_hits.txt RefSeqHits --OriginDatabase RefSeq --AddDb
python ProcessBLASTHits.py partial_hits.txt PartialHits  --OriginDatabase GenBank --AddDb
python ProcessBLASTHits.py complete_hits.txt CompleteHits --OriginDatabase GenBank --AddDb
# examples of how MergeBLASTHits.py was run
python MergeBLASTHits.py RefSeqHits.csv AllHits PartialHits.csv CompleteHits.csv --AddDb
python MergeBLASTHits.py PartialHits.csv GenBankHits CompleteHits.csv --AddDb
# example of how CheckMyBLAST.py can be run
python CheckMyBlast.py complete_test.txt testfilter.txt ProblematicVirus.csv ProblematicNonVirus.csv --Virus --Strict
