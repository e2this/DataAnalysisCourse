import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
import numpy as np
import cmapPy.pandasGEXpress.write_gctx as wg

sig_info = pd.read_csv("../Data/GSE92742_Broad_LINCS_sig_info.txt", sep="\t")

cell_list=sig_info['cell_id'].drop_duplicates()

def intersect(a, b):
    return list(set(a) & set(b))
	
for cell in cell_list:
	sig_info_cell_cp = sig_info[["sig_id","pert_id","pert_idose","pert_itime"]][(sig_info["cell_id"] == cell) & (sig_info["pert_type"] == "trt_cp")]
	sig_info_cell_vehicle = sig_info[["sig_id","pert_id","pert_idose","pert_itime"]][(sig_info["cell_id"] == cell) & (sig_info["pert_type"] == "ctl_vehicle")]
	if ((len(sig_info_cell_cp) != 0) & (len(sig_info_cell_vehicle) != 0)):
		pert_id_cp_list = list(set(sig_info_cell_cp["pert_id"]))
		x=sig_info_cell_cp["pert_id"]
		
		i=0
		while i < len(pert_id_cp_list):
			pert_ids=[]
			idx_pert_id=[]
			while (len(idx_pert_id) < 1000) & (i < len(pert_id_cp_list)) :
				temp = [k for k,item in enumerate(x) if (item in pert_id_cp_list[i])]
				idx_pert_id = idx_pert_id + list(temp)
				i=i+1
			
			pert_ids_names=list(set(sig_info_cell_cp["pert_id"].iloc[idx_pert_id]))
			pert_ids= sig_info_cell_cp["sig_id"].iloc[idx_pert_id]
			#pert_ids= sig_info_cell_cp["sig_id"][sig_info_cell_cp["pert_id"]== pert_id]
			
			dose_cp_list= list(set(sig_info_cell_cp["pert_idose"].iloc[idx_pert_id]))
			time_cp_list= list(set(sig_info_cell_cp["pert_itime"].iloc[idx_pert_id]))


			#x=sig_info_cell_vehicle["pert_idose"]
			#idx_idose=[i for i,item in enumerate(x) if (item in dose_cp_list)]
			#x=sig_info_cell_vehicle["pert_itime"]
			#idx_itime=[i for i,item in enumerate(x) if (item in time_cp_list)]
			#idx_both= intersect(idx_idose,idx_itime)
							
			#if( len(idx_both) != 0):
			#	ctl_ids = sig_info_cell_vehicle["sig_id"][idx_both]
			#elif( len (idx_idose) != 0):
			#	ctl_ids = sig_info_cell_vehicle["sig_id"][idx_idose]
			#else:
			ctl_ids = sig_info_cell_vehicle["sig_id"]
				
				
			tot_ids= ctl_ids.tolist() + pert_ids.tolist()
			print( "cell = %s, pert_id= %s \nN_pert=%d, N_ctl=%d" % (cell,"nothing",len(pert_ids),len(ctl_ids)))
			data = parse("../Data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", cid=tot_ids)
			wg.write(data, "Subset_GCT/GSE92742_Level5_%s_%d" % (cell,i))
		
			for pert_id in pert_ids_names:
				for dose in dose_cp_list:
					for time in time_cp_list:
						pert_ids_cls= sig_info_cell_cp["sig_id"][(sig_info_cell_cp["pert_id"]== pert_id)&(sig_info_cell_cp["pert_idose"]==dose) & (sig_info_cell_cp["pert_itime"]== time)]
						if(len(pert_ids_cls) != 0):
							ctl_ids_cls= sig_info_cell_vehicle["sig_id"][(sig_info_cell_vehicle["pert_id"]== pert_id)&(sig_info_cell_vehicle["pert_idose"]==dose) & (sig_info_cell_vehicle["pert_itime"]== time)]
							if(len(ctl_ids_cls) == 0):
								ctl_ids_cls= sig_info_cell_vehicle["sig_id"][(sig_info_cell_vehicle["pert_id"]== pert_id)&(sig_info_cell_vehicle["pert_idose"]==dose)]
								if(len(ctl_ids_cls) == 0):
									ctl_ids_cls= sig_info_cell_vehicle["sig_id"][(sig_info_cell_vehicle["pert_id"]== pert_id)]
									if(len(ctl_ids_cls) ==0):
										ctl_ids_cls= sig_info_cell_vehicle["sig_id"]
		
							sample_ids=data.col_metadata_df.index.tolist()
		
							q=["" for i in range(len(sample_ids))]
							v=["" for i in range(len(sample_ids))]

							for i in range(len(sample_ids)):
								if (sample_ids[i] in ctl_ids_cls.tolist()):
									v[i]='1'
									q[i]="Ctrl"
								elif (sample_ids[i] in pert_ids_cls.tolist()):
									v[i]='2'
									q[i]="Pert"
								else:
									v[i]='0'
									q[i]="None"

							q_unique=[]
							for item in q:
								if(item not in q_unique):
									q_unique=q_unique + [item]

							line2 = " ".join(q_unique)
							line2="# "+line2

							n_gr=len(q_unique)
							N=len(sample_ids)

							line1="%d %d 1" % (N,n_gr)

							line3=" ".join(v)
							
							f=open("CLS/GSE92742_Level5_%s_%s_%s_%s.cls" % (cell,pert_id,dose,time))
							f.write("%s\n%s\n%s" % (line1,line2,line3))
							f.close()
									