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


#ProcessBLASTHits.py
Takes in BLASTn results file and file of retroviruses ids. Processes the reuslts by splitting columns, flagging retroviruses, categorises where the potential contamination subsequence is within the virus query sequence, and uses GetTax.py functon get_tax to add taxonomy information.
The results are added to CSV by default, but can also be added to database.
If the option to AddDb is used without specifying, it will add to "Hits.db" by default but can be added to a database of a user-set name.

Examples of how ProcessBLASTHits.py was run in the HIVES analysis:
python ProcessBLASTHits.py refseq_hits.txt RefSeqHits --OriginDatabase RefSeq --AddDb
python ProcessBLASTHits.py partial_hits.txt PartialHits  --OriginDatabase GenBank --AddDb
python ProcessBLASTHits.py complete_hits.txt CompleteHits --OriginDatabase GenBank --AddDb

# MergeBLASTHits.py

Used to combine blastn results and save results together as one file or one database table. Takes in one main file, and at least one other file to merge with. 
The output name does not need a file extension as it will automatically written to csv. The output name is also used to create the database table for the merged results.
If the option to AddDb is used without specifying, it will add to "Hits.db" by default.


Examples of how MergeBLASTHits.py was run in the HIVES analysis:
python MergeBLASTHits.py RefSeqHits.csv AllHits PartialHits.csv CompleteHits.csv --AddDb
python MergeBLASTHits.py PartialHits.csv GenBankHits CompleteHits.csv --AddDb


# CheckMyBLAST.py
Works by taking in both ProblematicVirus.csv and ProblematicNonVirus.csv, generated by FinalFiltering.py
Has the option to remove Virus sequences,  Non-virus sequences or both.
If option Strict chosen, then it will filter out any contaminted sequences that were flagged to be removed from the dataset
Will automically create a summary file of any matching sequences in the input BLAST file, with why the sequence was flagged and the recommended action.

Example of how CheckMyBLAST.py can be used to filter out contaminated virus sequences from a BLAST results file:
python CheckMyBlast.py complete_test.txt testfilter.txt ProblematicVirus.csv ProblematicNonVirus.csv --Virus --Strict


