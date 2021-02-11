my_tumor_data = "MyTumor.tpm.txt"
data_dir = "./"
my_tumor_mat = data.table::fread(file.path(data_dir, my_tumor_data)) %>% as.data.frame()

