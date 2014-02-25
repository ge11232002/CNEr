##############   IO functions ###############




############## Read CNE from SQLite ########
dbName = "/mnt/biggley/home/gtan/work/projects/CNEr/CNErData/cne.sqlite"
tableName = "cne2wBf_danRer7_hg19_40_50"

chr = "chr16"
CNEstart = 45000000
CNEend = 60000000
minLength = 50
windowSize = 300

fetchedCNERanges = readCNERangesFromSQLite(dbName, tableName, chr, CNEstart, CNEend, whichAssembly="2", minLength)

cne2wBf_danRer7_hg19_40_50 = CNEAnnotate(dbName, tableName, whichAssembly="2", chr, CNEstart, CNEend, windowSize, minLength)
cne2wBf_danRer7_hg19_48_50 = CNEAnnotate(dbName, "cne2wBf_danRer7_hg19_48_50", whichAssembly="2", chr, CNEstart, CNEend, windowSize, minLength)
cne2wBf_danRer7_hg19_49_50= CNEAnnotate(dbName, "cne2wBf_danRer7_hg19_49_50", whichAssembly="2", chr, CNEstart, CNEend, windowSize, minLength)


########### Plot in Gviz ###############
genome = "hg19"
strand = "+"
chr = "chr16"
dataMatrix= cne2wBf_danRer7_hg19_40_50

dTrack1 = DataTrack(start=dataMatrix[ ,1], end=dataMatrix[ ,1], 
                    data=dataMatrix[ ,2], chromosome=chr, strand=strand, 
                    genome=genome, type="horiz", 

                    name="hg19_canFam3_len50_id900_300")


axisTrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
plotTracks(list(ideoTrack, axisTrack, dTrack1), collapseTranscripts = TRUE, shape = "arrow", from= 45000000, to=60000000, showId = TRUE, extend.left = 20000)


