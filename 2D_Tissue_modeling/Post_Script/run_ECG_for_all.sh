
run() {
	
	folder=$1

declare -a subfolders=(BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_K1  \
BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_NaV  \
BCL.500.ISO.0.1.CaMKII.CaMKII_db.No_CaMKII_GapG  \
BCL.500.ISO.0.CaMKII.CaMKII_inhb  \
BCL.500.ISO.0.CaMKII.CaMKII_db  \
BCL.500.ISO.0.CaMKII.Default  \
BCL.500.ISO.0.1.CaMKII.CaMKII_inhb  \
BCL.500.ISO.0.1.CaMKII.CaMKII_db  \
BCL.500.ISO.0.1.CaMKII.Default) 


for sub in ${subfolders[@]}
	do
		echo $sub
		./ECG 60 0 0 ../$folder/$sub 0 10000 > log
	done

}


run D_100%
run D_75%
run D_50%
run D_25%




